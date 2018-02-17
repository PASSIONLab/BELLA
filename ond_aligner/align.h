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
#ifndef ALIGN_H_
#define ALIGN_H_

/*
 * originally located in DB.h
 */
#define EPRINTF fprintf
#define EPLACE  stderr
#define EXIT(x) exit (1)

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef signed char        int8;
typedef signed short       int16;
typedef signed int         int32;
typedef signed long long   int64;

/*******************************************************************************************
 *
 *  UTILITIES
 *
 * originally in DB.h
 ********************************************************************************************/

//  The following general utilities return NULL if any of their input pointers are NULL, or if they
//    could not perform their function (in which case they also print an error to stderr).

void *Malloc(int64 size, char *mesg);      				//  Guarded versions of malloc, realloc              //  Guarded versions of malloc, realloc
void *Realloc(void *object, int64 size, char *mesg);     //  and strdup, that output "mesg" to

/*** PATH ABSTRACTION:

     Coordinates are *between* characters where 0 is the tick just before the first char,
     1 is the tick between the first and second character, and so on.  Our data structure
     is called a Path refering to its conceptualization in an edit graph.

     A local alignment is specified by the point '(abpos,bbpos)' at which its path in
     the underlying edit graph starts, and the point '(aepos,bepos)' at which it ends.
     In otherwords A[abpos+1..aepos] is aligned to B[bbpos+1..bepos] (assuming X[1] is
     the *first* character of X).

     There are 'diffs' differences in an optimal local alignment between the beginning and
     end points of the alignment (if computed by Compute_Trace), or nearly so (if computed
     by Local_Alignment).

     Optionally, a Path can have additional information about the exact nature of the
     aligned substrings if the field 'trace' is not NULL.  Trace points to either an
     array of integers (if computed by a Compute_Trace routine), or an array of unsigned
     short integers (if computed by Local_Alignment).

     If computed by Local_Alignment 'trace' points at a list of 'tlen' (always even) short
     values:

            d_0, b_0, d_1, b_1, ... d_n-1, b_n-1, d_n, b_n

     to be interpreted as follows.  The alignment from (abpos,bbpos) to (aepos,bepos)
     passes through the n trace points for i in [1,n]:

            (a_i,b_i) where a_i = floor(abpos/TS)*TS + i*TS
                        and b_i = bbpos + (b_0 + b_1 + b_i-1)

     where also let a_0,b_0 = abpos,bbpos and a_(n+1),b_(n+1) = aepos,bepos.  That is, the
     interior (i.e. i != 0 and i != n+1) trace points pass through every TS'th position of
     the aread where TS is the "trace spacing" employed when finding the alignment (see
     New_Align_Spec).  Typically TS is 100.  Then d_i is the number of differences in the
     portion of the alignment between (a_i,b_i) and (a_i+1,b_i+1).  These trace points allow
     the Compute_Trace routines to efficiently compute the exact alignment between the two
     reads by efficiently computing exact alignments between consecutive pairs of trace points.
     Moreover, the diff values give one an idea of the quality of the alignment along every
     segment of TS symbols of the aread.

     If computed by a Compute_Trace routine, 'trace' points at a list of 'tlen' integers
     < i1, i2, ... in > that encodes an exact alignment as follows.  A negative number j
     indicates that a dash should be placed before A[-j] and a positive number k indicates
     that a dash should be placed before B[k], where A and B are the two sequences of the
     overlap.  The indels occur in the trace in the order in which they occur along the
     alignment.  For a good example of how to "decode" a trace into an alignment, see the
     code for the routine Print_Alignment.

***/

typedef struct
  { void     *trace;
    int       tlen;
    int       diffs;
    int       abpos, bbpos;
    int       aepos, bepos;
  } Path;


/*** ALIGNMENT ABSTRACTION:

     An alignment is modeled by an Alignment record, which in addition to a *pointer* to a
     'path', gives pointers to the A and B sequences, their lengths, and indicates whether
     the B-sequence needs to be complemented ('comp' non-zero if so).  The 'trace' pointer
     of the 'path' subrecord can be either NULL, a list of pass-through points, or an exact
     trace depending on what routines have been called on the record.

     One can (1) compute a trace, with Compute_Trace, either from scratch if 'path.trace' = NULL,
     or using the sequence of pass-through points in trace, (2) print an ASCII representation
     of an alignment, or (3) reverse the roles of A and B, and (4) complement a sequence
     (which is a reversible process).

     If the alignment record shows the B sequence as complemented, *** THEN IT IS THE
     RESPONSIBILITY OF THE CALLER *** to make sure that bseq points at a complement of
     the sequence before calling Compute_Trace or Print_Alignment.  Complement_Seq complements
     the sequence a of length n.  The operation does the complementation/reversal in place.
     Calling it a second time on a given fragment restores it to its original state.

     With the introduction of the DAMAPPER, we need to code chains of alignments between a
     pair of sequences.  The alignments of a chain are expected to be found in order either on
     a file or in memory, where the START_FLAG marks the first alignment and the NEXT_FLAG all
     subsequent alignmenst in a chain.  A chain of a single LA is marked with the START_FLAG.
     The BEST_FLAG marks one of the best chains for a pair of sequences.  The convention is
     that either every record has either a START- or NEXT-flag, or none of them do (e.g. as
     produced by daligner), so one can always check the flags of the first alignment to see
     whether or not the chain concept applies to a given collection or not.
***/

#define COMP_FLAG  0x1
#define ACOMP_FLAG 0x2   //  A-sequence is complemented, not B !  Only Local_Alignment notices

#define COMP(x)   ((x) & COMP_FLAG)
#define ACOMP(x)  ((x) & ACOMP_FLAG)

typedef struct
	{ Path   *path;
	  uint32  flags;        /* Pipeline status and complementation flags          */
	  char   *aseq;         /* Pointer to A sequence                              */
	  char   *bseq;         /* Pointer to B sequence                              */
	  int     alen;         /* Length of A sequence                               */
	  int     blen;         /* Length of B sequence                               */
	} Alignment;

  /* Many routines like Local_Alignment, Compute_Trace, and Print_Alignment need working
	 storage that is more efficiently reused with each call, rather than being allocated anew
	 with each call.  Each *thread* can create a Work_Data object with New_Work_Data and this
	 object holds and retains the working storage for routines of this module between calls
	 to the routines.  If enough memory for a Work_Data is not available then NULL is returned.
	 Free_Work_Data frees a Work_Data object and all working storage held by it.
  */

  typedef void Work_Data;

  Work_Data *New_Work_Data();

  void       Free_Work_Data(Work_Data *work);

  /* Local_Alignment seeks local alignments of a quality determined by a number of parameters.
	 These are coded in an Align_Spec object that can be created with New_Align_Spec and
	 freed with Free_Align_Spec when no longer needed.  There are 4 essential parameters:

	 ave_corr:    the average correlation (1 - 2*error_rate) for the sought alignments.  For Pacbio
					data we set this to .70 assuming an average of 15% error in each read.
	 trace_space: the spacing interval for keeping trace points and segment differences (see
					description of 'trace' for Paths above)
	 freq[4]:     a 4-element vector where afreq[0] = frequency of A, f(A), freq[1] = f(C),
					freq[2] = f(G), and freq[3] = f(T).  This vector is part of the header
					of every DAZZ database (see db.h).
	 reach:       a boolean, if set alignment extend to the boundary when reasonable, otherwise
					the terminate only at suffix-positive points.

	 If an alignment cannot reach the boundary of the d.p. matrix with this condition (i.e.
	 overlap), then the last/first 30 columns of the alignment are guaranteed to be
	 suffix/prefix positive at correlation ave_corr * g(freq) where g is an empirically
	 measured function that increases from 1 as the entropy of freq decreases.  If memory is
	 unavailable or the freq distribution is too skewed then NULL is returned.

	 You can get back the original parameters used to create an Align_Spec with the simple
	 utility functions below.
  */

  typedef void Align_Spec;

  Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq, int reach);

  void        Free_Align_Spec(Align_Spec *spec);

/* Local_Alignment finds the longest significant local alignment between the sequences in
 'align' subject to:

   (a) the alignment criterion given by the Align_Spec 'spec',
   (b) it passes through one of the points (anti+k)/2,(anti-k)/2 for k in [low,hgh] within
		 the underlying dynamic programming matrix (i.e. the points on diagonals low to hgh
		 on anti-diagonal anti or anti-1 (depending on whether the diagonal is odd or even)),
   (c) if lbord >= 0, then the alignment is always above diagonal low-lbord, and
   (d) if hbord >= 0, then the alignment is always below diagonal hgh+hbord.

 The path record of 'align' has its 'trace' filled from the point of view of an overlap
 between the aread and the bread.  In addition a Path record from the point of view of the
 bread versus the aread is returned by the function, with this Path's 'trace' filled in
 appropriately.  The space for the returned path and the two 'trace's are in the working
 storage supplied by the Work_Data packet and this space is reused with each call, so if
 one wants to retain the bread-path and the two trace point sequences, then they must be
 copied to user-allocated storage before calling the routine again.  NULL is returned in
 the event of an error.

 Find_Extension is a variant of Local_Alignment that simply finds a local alignment that
 either ends (if prefix is non-zero) or begins (if prefix is zero) at the point
 (anti+diag)/2,(anti-diag)/2).  All other parameters are as before.  It returns a non-zero
 value only when INTERACTIVE is on and it cannot allocate the memory it needs.
 Only the path and trace with respect to the aread is returned.  This routine is experimental
 and may not persist in later versions of the code.
*/

Path *Local_Alignment(Alignment *align, Work_Data *work, Align_Spec *spec,
					int low, int hgh, int anti, int lbord, int hbord);

int   Find_Extension(Alignment *align, Work_Data *work, Align_Spec *spec,    //  experimental !!
                     int diag, int anti, int lbord, int hbord, int prefix);

#endif /* ALIGN_H_ */
