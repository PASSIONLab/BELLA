/*******************************************************************************************
 *	This code was extracted from DALIGNER. Only the code essential for the local alignment
 *	functionality (independent of the DAZZ_DB or other infrastructure) is retained herein.
 *	The code herein has not been modified over the original. See original header comments below.
 *
 *	Extracted by :  Marquita M. Ellis
 *	Date			:  February 2018
 *
 *  Local alignment module.  Routines for finding local alignments given a seed position,
 *    representing such an l.a. with its interval and a set of pass-thru points, so that
 *    a detailed alignment can be efficiently computed on demand.
 *
 *  All routines work on a numeric representation of DNA sequences, i.e. 0 for A, 1 for C,
 *    2 for G, and 3 for T.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include "align.h"

/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 *  originally located in DB.c
 ********************************************************************************************/
char Prog_Name[] = "O(Nd) alignment kernel";

void *Malloc(int64 size, char *mesg)
{ void *p;

  if ((p = malloc(size)) == NULL)
    { if (mesg == NULL)
        EPRINTF(EPLACE,"%s: Out of memory\n",Prog_Name);
      else
        EPRINTF(EPLACE,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

void *Realloc(void *p, int64 size, char *mesg)
{ if (size <= 0)
    size = 1;
  if ((p = realloc(p,size)) == NULL)
    { if (mesg == NULL)
        EPRINTF(EPLACE,"%s: Out of memory\n",Prog_Name);
      else
        EPRINTF(EPLACE,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

/****************************************************************************************\
*                                                                                        *
*  Working Storage Abstraction                                                           *
*                                                                                        *
\****************************************************************************************/

typedef struct            //  Hidden from the user, working space for each thread
  { int     vecmax;
    void   *vector;
    int     celmax;
    void   *cells;
    int     pntmax;
    void   *points;
    int     tramax;
    void   *trace;
  } _Work_Data;

Work_Data *New_Work_Data()
{ _Work_Data *work;

  work = (_Work_Data *) Malloc(sizeof(_Work_Data),"Allocating work data block");
  if (work == NULL)
    EXIT(NULL);
  work->vecmax = 0;
  work->vector = NULL;
  work->pntmax = 0;
  work->points = NULL;
  work->tramax = 0;
  work->trace  = NULL;
  work->celmax = 0;
  work->cells  = NULL;
  return ((Work_Data *) work);
}

static int enlarge_vector(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->vector,max,"Enlarging DP vector");
  if (vec == NULL)
    EXIT(1);
  work->vecmax = max;
  work->vector = vec;
  return (0);
}

static int enlarge_points(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->points,max,"Enlarging point vector");
  if (vec == NULL)
    EXIT(1);
  work->pntmax = max;
  work->points = vec;
  return (0);
}

void Free_Work_Data(Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  if (work->vector != NULL)
    free(work->vector);
  if (work->cells != NULL)
    free(work->cells);
  if (work->trace != NULL)
    free(work->trace);
  if (work->points != NULL)
    free(work->points);
  free(work);
}

/****************************************************************************************\
*                                                                                        *
*  ADAPTIVE PATH FINDING                                                                 *
*                                                                                        *
\****************************************************************************************/

  //  Absolute/Fixed Parameters

#define BVEC  uint64     //  Can be uint32 if PATH_LEN <= 32

#define TRIM_LEN    15   //  Report as the tip, the last wave maximum for which the last
                         //     2*TRIM_LEN edits are prefix-positive at rate ave_corr*f(bias)
                         //     (max value is 20)

#define PATH_LEN    60   //  Follow the last PATH_LEN columns/edges (max value is 63)

  //  Derivative fixed parameters

#define PATH_TOP  0x1000000000000000ll   //  Must be 1 << PATH_LEN
#define PATH_INT  0x0fffffffffffffffll   //  Must be PATH_TOP-1
#define TRIM_MASK 0x7fff                 //  Must be (1 << TRIM_LEN) - 1
#define TRIM_MLAG 200                    //  How far can last trim point be behind best point
#define WAVE_LAG   30                    //  How far can worst point be behind the best point

static double Bias_Factor[10] = { .690, .690, .690, .690, .780,
                                  .850, .900, .933, .966, 1.000 };
//  Adjustable paramters

typedef struct
{ double ave_corr;
  int    trace_space;
  int    reach;
  float  freq[4];
  int    ave_path;
  int16 *score;
  int16 *table;
} _Align_Spec;

/* Fill in bit table: TABLE[x] = 1 iff the alignment modeled by x (1 = match, 0 = mismatch)
     has a non-negative score for every suffix of the alignment under the scoring scheme
     where match = MATCH and mismatch = -1.  MATCH is set so that an alignment with TRIM_PCT
     matches has zero score ( (1-TRIM_PCT) / TRIM_PCT ).                                     */

#define FRACTION 1000  //  Implicit fractional part of scores, i.e. score = x/FRACTION

typedef struct
  { int    mscore;
    int    dscore;
    int16 *table;
    int16 *score;
  } Table_Bits;

static void set_table(int bit, int prefix, int score, int max, Table_Bits *parms)
{ if (bit >= TRIM_LEN)
    { parms->table[prefix] = (int16) (score-max);
      parms->score[prefix] = (int16) score;
    }
  else
    { if (score > max)
        max = score;
      set_table(bit+1,(prefix<<1),score - parms->dscore,max,parms);
      set_table(bit+1,(prefix<<1) | 1,score + parms->mscore,max,parms);
    }
}

/* Create an alignment specification record including path tip tables & values */

Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq, int reach)
{ _Align_Spec *spec;
  Table_Bits   parms;
  double       match;
  int          bias;

  spec = (_Align_Spec *) Malloc(sizeof(_Align_Spec),"Allocating alignment specification");
  if (spec == NULL)
    EXIT(NULL);

  spec->ave_corr    = ave_corr;
  spec->trace_space = trace_space;
  spec->reach       = reach;
  spec->freq[0]     = freq[0];
  spec->freq[1]     = freq[1];
  spec->freq[2]     = freq[2];
  spec->freq[3]     = freq[3];

  match = freq[0] + freq[3];
  if (match > .5)
    match = 1.-match;
  bias = (int) ((match+.025)*20.-1.);
  if (match < .2)
    { fprintf(stderr,"Warning: Base bias worse than 80/20%% ! (New_Align_Spec)\n");
      fprintf(stderr,"         Capping bias at this ratio.\n");
      bias = 3;
    }

  spec->ave_path = (int) (PATH_LEN * (1. - Bias_Factor[bias] * (1. - ave_corr)));
  parms.mscore   = (int) (FRACTION * Bias_Factor[bias] * (1. - ave_corr));
  parms.dscore   = FRACTION - parms.mscore;

  parms.score = (int16 *) Malloc(sizeof(int16)*(TRIM_MASK+1)*2,"Allocating trim table");
  if (parms.score == NULL)
    { free(spec);
      EXIT(NULL);
    }
  parms.table = parms.score + (TRIM_MASK+1);

  set_table(0,0,0,0,&parms);

  spec->table = parms.table;
  spec->score = parms.score;

  return ((Align_Spec *) spec);
}

void Free_Align_Spec(Align_Spec *espec)
{ _Align_Spec *spec = (_Align_Spec *) espec;
  free(spec->score);
  free(spec);
}

//...

/****************************************************************************************\
*                                                                                        *
*  LOCAL ALIGNMENT FINDER: forward_/reverse_wave and Local_Alignment                     *
*                                                                                        *
\****************************************************************************************/
/* At each furthest reaching point, keep a-coordinate of point (V), bitvector
     recording the last TRIM_LEN columns of the implied alignment (T), and the
     # of matches (1-bits) in the bitvector (M).                               */

typedef struct
  { int ptr;
    int diag;
    int diff;
    int mark;
  } Pebble;

static int VectorEl = 6*sizeof(int) + sizeof(BVEC);

static int forward_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int *mind, int maxd, int mida, int minp, int maxp, int aoff, int boff)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int     REACH       = spec->reach;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = maxd;
  low = *mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  /* Compute 0-wave starting from mid-line */

  more  = 1;
  aclip =  INT32_MAX;
  bclip = -INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a  = aseq + hgh;
    for (k = hgh; k >= low; k--)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = (((y+k)+(TRACE_SPACE-aoff))/TRACE_SPACE-1)*TRACE_SPACE+aoff;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = na;
        ha  = avail++;
        na += TRACE_SPACE;

        nb = ((y+(TRACE_SPACE-boff))/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: %d,%d,0,%d\n",avail,-1,k,nb); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = nb;
        hb  = avail++;
        nb += TRACE_SPACE;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip < k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4)
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y += 1;
          }
        c = (y << 1) + k;

        while (y+k >= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }
        while (y >= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = nb;
            hb  = avail++;
            nb += TRACE_SPACE;
          }

        if (c > besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a -= 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (hgh >= aclip)
        { hgh = aclip-1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (low <= bclip)
        { low = bclip+1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip =  INT32_MAX;
      bclip = -INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nFORWARD WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  /* Compute successive waves until no furthest reaching points remain */

  while (more && lasta >= besta - TRIM_MLAG)
    { int     k, n;
      int     ua, ub;
      BVEC    t;
      int     am, ac, ap;
      char   *a;

      low -= 1;
      hgh += 1;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move;
          int64 vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low >= minp)
        { NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = -1;
        }
      else
        low += 1;

      if (hgh <= maxp)
        { NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh]  = am = -1;
        }
      else
        am = V[--hgh];

      dif += 1;

      ac = V[hgh+1] = V[low-1] = -1;
      a  = aseq + hgh;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = hgh; k >= low; k--)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          ap = ac;
          ac = am;
          am = V[d = k-1];

          if (ac < am)
            if (am < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = am+1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac+2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip < k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4)
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y += 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k >= NA[k])
            { if (cells[ha].mark < NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] += TRACE_SPACE;
            }

          while (y >= NB[k])
            { if (cells[hb].mark < NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] += TRACE_SPACE;
            }

          if (c > besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a -= 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta-besty] != 4)
            more = 1;
          if (hgh >= aclip)
            { hgh = aclip-1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (low <= bclip)
            { low = bclip+1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip =  INT32_MAX;
          bclip = -INT32_MAX;
        }

      n = besta - WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] < n)
          hgh -= 1;
        else
          { while (V[low] < n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0 && REACH)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida-k)/2;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        d = cells[h].diff;
        atrace[atlen++] = (uint16) (d-e);
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimd-e);
        atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
        atrace[atlen-2] = (uint16) (atrace[atlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida+k)/2;
    e = 0;
    low = k;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,(mida-k)/2); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark + k;
        d = cells[h].diff;
        btrace[btlen++] = (uint16) (d-e);
        btrace[btlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b-k != trimy)
      { btrace[btlen++] = (uint16) (trimd-e);
        btrace[btlen++] = (uint16) (trimx-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }
    else if (b != trimx)
      { btrace[btlen-1] = (uint16) (btrace[btlen-1] + (trimx-b));
        btrace[btlen-2] = (uint16) (btrace[btlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }

    apath->aepos = trimx;
    apath->bepos = trimy;
    apath->diffs = trimd;
    apath->tlen  = atlen;
    bpath->tlen  = btlen;
  }

  *mind = low;
  return (0);
}

/*** Reverse Wave ***/

static int reverse_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int mind, int maxd, int mida, int minp, int maxp, int aoff, int boff)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int     REACH       = spec->reach;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = maxd;
  low = mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  more  = 1;
  aclip = -INT32_MAX;
  bclip =  INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a = aseq + low;
    for (k = low; k <= hgh; k++)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = (((y+k)+(TRACE_SPACE-aoff)-1)/TRACE_SPACE-1)*TRACE_SPACE+aoff;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y+k;
        ha  = avail++;

        nb = ((y+(TRACE_SPACE-boff)-1)/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: -1,%d,0,%d\n",avail,k,nb+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y;
        hb  = avail++;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip > k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4)
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y -= 1;
          }
        c = (y << 1) + k;

        while (y+k <= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }
        while (y <= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = nb;
            hb  = avail++;
            nb -= TRACE_SPACE;
          }

        if (c < besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a += 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (low <= aclip)
        { low = aclip+1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (hgh >= bclip)
        { hgh = bclip-1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip = -INT32_MAX;
      bclip =  INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nREVERSE WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  while (more && lasta <= besta + TRIM_MLAG)
    { int    k, n;
      int    ua, ub;
      BVEC   t;
      int    am, ac, ap;
      char  *a;

      low -= 1;
      hgh += 1;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move, vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low >= minp)
        { NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[++low];

      if (hgh <= maxp)
        { NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh] = INT32_MAX;
        }
      else
        hgh -= 1;

      dif += 1;

      ac = V[hgh+1] = V[low-1] = INT32_MAX;
      a  = aseq + low;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = low; k <= hgh; k++)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          am = ac;
          ac = ap;
          ap = V[d = k+1];

          if (ac > ap)
            if (ap > am)
              { c = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ap-1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac > am)
              { c  = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac-2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip > k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4)
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y -= 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k <= NA[k])
            { if (cells[ha].mark > NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] -= TRACE_SPACE;
            }
          while (y <= NB[k])
            { if (cells[hb].mark > NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] -= TRACE_SPACE;
            }

          if (c < besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a += 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
            more = 1;
          if (low <= aclip)
            { low = aclip+1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (hgh >= bclip)
            { hgh = bclip-1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip = -INT32_MAX;
          bclip =  INT32_MAX;
        }

      n = besta + WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] > n)
          hgh -= 1;
        else
          { while (V[low] > n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0 && REACH)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark - k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != aoff)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimy;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
        if (apath->tlen == 0)
          { atrace[--atlen] = (uint16) (b-a);
            atrace[--atlen] = (uint16) (d-e);
          }
        else
          { atrace[1] = (uint16) (atrace[1] + (b-a));
            atrace[0] = (uint16) (atrace[0] + (d-e));
          }
        b = a;
        e = d;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (uint16) (b-a);
            d = cells[h].diff;
            atrace[--atlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (uint16) (b-trimy);
            atrace[--atlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen+1] = (uint16) (atrace[atlen+1] + (b-trimy));
            atrace[atlen]   = (uint16) (atrace[atlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark + k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,b-k); fflush(stdout);
#endif
    if ((b-k)%TRACE_SPACE != boff)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimx;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark + k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
        if (bpath->tlen == 0)
          { btrace[--btlen] = (uint16) (b-a);
            btrace[--btlen] = (uint16) (b-a);
          }
        else
          { btrace[1] = (uint16) (btrace[1] + (b-a));
            btrace[0] = (uint16) (btrace[0] + (d-e));
          }
        b = a;
        e = d;
      }

    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark + k;
            btrace[--btlen] = (uint16) (b-a);
            d = cells[h].diff;
            btrace[--btlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b-k != trimy)
          { btrace[--btlen] = (uint16) (b-trimx);
            btrace[--btlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
        else if (b != trimx)
          { btrace[btlen+1] = (uint16) (btrace[btlen+1] + (b-trimx));
            btrace[btlen]   = (uint16) (btrace[btlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
      }

    apath->abpos = trimx;
    apath->bbpos = trimy;
    apath->diffs = apath->diffs + trimd;
    apath->tlen  = apath->tlen  - atlen;
    apath->trace = atrace + atlen;
    bpath->tlen  = bpath->tlen  - btlen;
    bpath->trace = btrace + btlen;
  }

  return (0);
}


/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

Path *Local_Alignment(Alignment *align, Work_Data *ework, Align_Spec *espec,
                      int low, int hgh, int anti, int lbord, int hbord)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath, *bpath;
  int   aoff, boff;
  int   minp, maxp;
  int   selfie;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    if (hgh-low >= 7500)
      wsize = VectorEl*(hgh-low+1);
    else
      wsize = VectorEl*10000;
    if (wsize >= work->vecmax)
      if (enlarge_vector(work,wsize))
        EXIT(NULL);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 4*maxtp*sizeof(uint16) + sizeof(Path);
    if (wsize > work->pntmax)
      if (enlarge_points(work,wsize))
        EXIT(NULL);

    apath = align->path;
    bpath = (Path *) work->points;

    apath->trace = ((uint16 *) (bpath+1)) + maxtp;
    bpath->trace = ((uint16 *) apath->trace) +  2*maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  selfie = (align->aseq == align->bseq);

  if (lbord < 0)
    { if (selfie && low >= 0)
        minp = 1;
      else
        minp = -INT32_MAX;
    }
  else
    minp = low-lbord;
  if (hbord < 0)
    { if (selfie && hgh <= 0)
        maxp = -1;
      else
        maxp = INT32_MAX;
    }
  else
    maxp = hgh+hbord;

  if (ACOMP(align->flags))
    { aoff = align->alen % spec->trace_space;
      boff = 0;
    }
  else if (COMP(align->flags))
    { aoff = 0;
      boff = align->blen % spec->trace_space;
    }
  else
    { aoff = 0;
      boff = 0;
    }

  if (forward_wave(work,spec,align,bpath,&low,hgh,anti,minp,maxp,aoff,boff))
    EXIT(NULL);

#ifdef DEBUG_PASSES
  printf("F1 (%d,%d) ~ %d => (%d,%d) %d\n",
         (2*anti+(low+hgh))/4,(anti-(low+hgh))/4,hgh-low,
         apath->aepos,apath->bepos,apath->diffs);
#endif

  if (reverse_wave(work,spec,align,bpath,low,low,anti,minp,maxp,aoff,boff))
    EXIT(NULL);

#ifdef DEBUG_PASSES
  printf("R1 (%d,%d) => (%d,%d) %d\n",
         (anti+low)/2,(anti-low)/2,apath->abpos,apath->bbpos,apath->diffs);
#endif

  bpath->diffs = apath->diffs;
  if (ACOMP(align->flags))
    { uint16 *trace = (uint16 *) apath->trace;
      uint16  p;
      int     i, j;

      bpath->aepos = apath->bepos;
      bpath->bepos = apath->aepos;
      bpath->abpos = apath->bbpos;
      bpath->bbpos = apath->abpos;

      apath->abpos = align->alen - bpath->bepos;
      apath->bbpos = align->blen - bpath->aepos;
      apath->aepos = align->alen - bpath->bbpos;
      apath->bepos = align->blen - bpath->abpos;
      i = apath->tlen-2;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          p = trace[i+1];
          trace[i+1] = trace[j+1];
          trace[j+1] = p;
          i -= 2;
          j += 2;
        }
    }
  else if (COMP(align->flags))
    { uint16 *trace = (uint16 *) bpath->trace;
      uint16  p;
      int     i, j;

      bpath->abpos = align->blen - apath->bepos;
      bpath->bbpos = align->alen - apath->aepos;
      bpath->aepos = align->blen - apath->bbpos;
      bpath->bepos = align->alen - apath->abpos;
      i = bpath->tlen-2;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          p = trace[i+1];
          trace[i+1] = trace[j+1];
          trace[j+1] = p;
          i -= 2;
          j += 2;
        }
    }
  else
    { bpath->aepos = apath->bepos;
      bpath->bepos = apath->aepos;
      bpath->abpos = apath->bbpos;
      bpath->bbpos = apath->abpos;
    }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",((COMP(align->flags) || ACOMP(align->flags)) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 1; h < apath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }

  { uint16 *trace = (uint16 *) bpath->trace;
    int     a, h;

    printf("\nB-path (%d,%d)->(%d,%d)",bpath->abpos,bpath->bbpos,bpath->aepos,bpath->bepos);
    printf(" %c [%d,%d]\n",((COMP(align->flags) || ACOMP(align->flags)) ? 'c' : 'n'),
                           align->blen,align->alen);
    a = bpath->bbpos;
    for (h = 1; h < bpath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }
#endif

  return (bpath);
}


/****************************************************************************************\
*                                                                                        *
*  EXTENSION VERSION OF LOCAL ALIGNMENT                                                  *
*                                                                                        *
\****************************************************************************************/

static int VectorEn = 4*sizeof(int) + sizeof(BVEC);

static int forward_extend(_Work_Data *work, _Align_Spec *spec, Alignment *align,
                          int midd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *NA;
  int    *_HA, *_NA;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha;
  int     morea, morey, mored;
  int     moreha;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = midd;
  low = midd;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEn;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _NA = _HA + vlen;
    _T  = ((BVEC *) (_NA + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    NA = _NA-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  /* Compute 0-wave starting from mid-line */

  more  = 1;
  aclip =  INT32_MAX;
  bclip = -INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  morem  = -1;

  { int   k;
    char *a;

    a  = aseq + hgh;
    for (k = hgh; k >= low; k--)
      { int     y, c, d;
        int     ha, na;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k)/TRACE_SPACE)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = na;
        ha  = avail++;
        na += TRACE_SPACE;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip < k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4)
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y += 1;
          }
        c = (y << 1) + k;

        while (y+k >= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }

        if (c > besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        NA[k] = na;

        a -= 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (hgh >= aclip)
        { hgh = aclip-1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
            }
        }
      if (low <= bclip)
        { low = bclip+1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
            }
        }
      aclip =  INT32_MAX;
      bclip = -INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nFORWARD WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  /* Compute successive waves until no furthest reaching points remain */

  while (more && lasta >= besta - TRIM_MLAG)
    { int     k, n;
      int     ua;
      BVEC    t;
      int     am, ac, ap;
      char   *a;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move;
          int64 vd, md, had, nad, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEn))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEn;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _NA = _HA + vlen;
              _T  = ((BVEC *) (_NA + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          NA = _NA-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          V[low]  = -1;
        }
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          V[hgh]  = am = -1;
        }
      else
        am = V[hgh];
      dif += 1;

      ac = V[hgh+1] = V[low-1] = -1;
      a  = aseq + hgh;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = -1;
      for (k = hgh; k >= low; k--)
        { int     y, m;
          int     ha;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          ap = ac;
          ac = am;
          am = V[d = k-1];

          if (ac < am)
            if (am < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = am+1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
              }
          else
            if (ac < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ac+2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip < k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4)
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y += 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k >= NA[k])
            { if (cells[ha].mark < NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] += TRACE_SPACE;
            }

          if (c > besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;

          a -= 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta-besty] != 4)
            more = 1;
          if (hgh >= aclip)
            { hgh = aclip-1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                }
            }
          if (low <= bclip)
            { low = bclip+1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                }
            }
          aclip =  INT32_MAX;
          bclip = -INT32_MAX;
        }

      n = besta - WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] < n)
          hgh -= 1;
        else
          { while (V[low] < n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    int     atlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
      }
    else
      trimx = trima-trimy;

    atlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida-k)/2;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        d = cells[h].diff;
        atrace[atlen++] = (uint16) (d-e);
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimd-e);
        atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
        atrace[atlen-2] = (uint16) (atrace[atlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }

    apath->aepos = trimx;
    apath->bepos = trimy;
    apath->diffs = trimd;
    apath->tlen  = atlen;
  }

  return (0);
}

static int reverse_extend(_Work_Data *work, _Align_Spec *spec, Alignment *align,
                          int midd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *NA;
  int    *_HA, *_NA;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha;
  int     morea, morey, mored;
  int     moreha;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = midd;
  low = midd;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEn;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _NA = _HA + vlen;
    _T  = ((BVEC *) (_NA + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    NA = _NA-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  more  = 1;
  aclip = -INT32_MAX;
  bclip =  INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  morem  = -1;

  { int   k;
    char *a;

    a = aseq + low;
    for (k = low; k <= hgh; k++)
      { int     y, c, d;
        int     ha, na;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k+TRACE_SPACE-1)/TRACE_SPACE-1)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y+k;
        ha  = avail++;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip > k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4)
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y -= 1;
          }
        c = (y << 1) + k;

        while (y+k <= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }

        if (c < besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        NA[k] = na;

        a += 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (low <= aclip)
        { low = aclip+1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
            }
        }
      if (hgh >= bclip)
        { hgh = bclip-1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
            }
        }
      aclip = -INT32_MAX;
      bclip =  INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nREVERSE WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  while (more && lasta <= besta + TRIM_MLAG)
    { int    k, n;
      int    ua;
      BVEC   t;
      int    am, ac, ap;
      char  *a;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move, vd, md, had, nad, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEn))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEn;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _NA = _HA + vlen;
              _T  = ((BVEC *) (_NA + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          NA = _NA-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[low];
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          V[hgh] = INT32_MAX;
        }
      dif += 1;

      ac = V[hgh+1] = V[low-1] = INT32_MAX;
      a  = aseq + low;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = -1;
      for (k = low; k <= hgh; k++)
        { int     y, m;
          int     ha;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          am = ac;
          ac = ap;
          ap = V[d = k+1];

          if (ac > ap)
            if (ap > am)
              { c = am-1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ap-1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
              }
          else
            if (ac > am)
              { c  = am-1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ac-2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip > k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4)
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y -= 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k <= NA[k])
            { if (cells[ha].mark > NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] -= TRACE_SPACE;
            }

          if (c < besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;

          a += 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
            more = 1;
          if (low <= aclip)
            { low = aclip+1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                }
            }
          if (hgh >= bclip)
            { hgh = bclip-1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                }
            }
          aclip = -INT32_MAX;
          bclip =  INT32_MAX;
        }

      n = besta + WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] > n)
          hgh -= 1;
        else
          { while (V[low] > n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    int     atlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
      }
    else
      trimx = trima-trimy;

    atlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr;
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark - k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != 0)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimy;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
        atrace[--atlen] = (uint16) (b-a);
        atrace[--atlen] = (uint16) (d-e);
        b = a;
        e = d;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (uint16) (b-a);
            d = cells[h].diff;
            atrace[--atlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (uint16) (b-trimy);
            atrace[--atlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen+1] = (uint16) (atrace[atlen+1] + (b-trimy));
            atrace[atlen]   = (uint16) (atrace[atlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
      }

    apath->abpos = trimx;
    apath->bbpos = trimy;
    apath->diffs = trimd;
    apath->tlen  = - atlen;
    apath->trace = atrace + atlen;
  }

  return (0);
}

/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

int Find_Extension(Alignment *align, Work_Data *ework, Align_Spec *espec,
                   int diag, int anti, int lbord, int hbord, int prefix)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath;
  int   minp, maxp;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    wsize = VectorEn*10000;
    if (wsize >= work->vecmax)
      if (enlarge_vector(work,wsize))
        EXIT(1);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 2*maxtp*sizeof(uint16);
    if (wsize > work->pntmax)
      if (enlarge_points(work,wsize))
        EXIT(1);

    apath = align->path;
    apath->trace = ((uint16 *) work->points) + maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  if (lbord < 0)
    minp = -INT32_MAX;
  else
    minp = diag-lbord;
  if (hbord < 0)
    maxp = INT32_MAX;
  else
    maxp = diag+hbord;

  if (prefix)
    { if (reverse_extend(work,spec,align,diag,anti,minp,maxp))
        EXIT(1);
      apath->aepos = (anti+diag)/2;
      apath->bepos = (anti-diag)/2;
#ifdef DEBUG_PASSES
      printf("E1 (%d,%d) => (%d,%d) %d\n",
             (anti+diag)/2,(anti-diag)/2,apath->abpos,apath->bbpos,apath->diffs);
#endif
    }
  else
    { if (forward_extend(work,spec,align,diag,anti,minp,maxp))
        EXIT(1);
      apath->abpos = (anti+diag)/2;
      apath->bbpos = (anti-diag)/2;
#ifdef DEBUG_PASSES
      printf("F1 (%d,%d) => (%d,%d) %d\n",
             (anti+diag)/2,(anti-diag)/2,apath->aepos,apath->bepos,apath->diffs);
#endif
     }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",(COMP(align->flags) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 1; h < apath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }
#endif

  return (0);
}
