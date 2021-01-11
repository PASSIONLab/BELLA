/*************************************************************************/
/*************************************************************************/
/* A PARALLEL RANDOM NUMBER GENERATION SYSTEM IN SPRNG FORMAT            */
/*                                                                       */ 
/* Author: Ashok Srinivasan,                                             */
/*            NCSA, University of Illinois, Urbana-Champaign             */
/* E-Mail: ashoks@ncsa.uiuc.edu                                          */
/*                                                                       */ 
/* Disclaimer: NCSA expressly disclaims any and all warranties, expressed*/
/* or implied, concerning the enclosed software.  The intent in sharing  */
/* this software is to promote the productive interchange of ideas       */
/* throughout the research community. All software is furnished on an    */
/* "as is" basis. No further updates to this software should be          */
/* expected. Although this may occur, no commitment exists. The authors  */
/* certainly invite your comments as well as the reporting of any bugs.  */
/* NCSA cannot commit that any or all bugs will be fixed.                */
/*                                                                       */
/* Note: Data stored by 'pack_rng' is NOT in a machine independent       */
/*       format. For that, please take a look at some SPRNG examples     */
/*       (lcg/lcg.c, lfg/lfg.c, etc).                                    */
/*************************************************************************/
/*************************************************************************/

/*             This is version 0.2, created 13 Apr 1998                  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "memory.h"
#include "interface.h"

/*** Change this to the type of generator you are implementing ***/
#define GENTYPE "Sample Generator"

#define NPARAMS 1		/*** number of valid parameters ***/
int MAX_STREAMS = (1<<19);  /*** Maximum number of independent streams ***/

struct rngen
{
  char *gentype;
  int stream_number;
  int nstreams;
  int init_seed;
  int parameter;
  int narrays;
  int *array_sizes;
  int **arrays;
  int spawn_offset;
  /*** declare other variables here ***/
  double x;			/* please remove this in your implementation */
};

int NGENS=0;		  /* number of random streams in current process */





/* Initialize random number stream */

#ifdef __STDC__
int *init_rng( int gennum, int total_gen,  int seed, 
	       int param)
#else
int *init_rng(gennum,total_gen,seed,param)
int gennum,param,seed,total_gen;
#endif
{
/*      gives back one stream (node gennum) with updated spawning         */
/*      info; should be called total_gen times, with different value      */
/*      of gennum in [0,total_gen) each call                              */
  struct rngen *genptr;
  int i;
  
  if (total_gen <= 0) /* Is total_gen valid ? */
  {
    total_gen = 1;
    fprintf(stderr,"WARNING - init_rng: Total_gen <= 0. Default value of 1 used for total_gen\n");
  }

  if (gennum >= MAX_STREAMS) /* check if gen_num is valid    */
    fprintf(stderr,"WARNING - init_rng: gennum: %d > maximum number of independent streams: %d\n\tIndependence of streams cannot be guranteed.\n",
	    gennum, MAX_STREAMS); 

  if (gennum < 0 || gennum >= total_gen) /* check if gen_num is valid    */
  {
    fprintf(stderr,"ERROR - init_rng: gennum %d out of range [%d,%d).\n",
	    gennum, 0, total_gen); 
    return (int *) NULL;
  }

  if (param < 0 || param >= NPARAMS)     /* check if parameter is valid */
  {
    fprintf(stderr,"WARNING - init_rng: parameter not valid. Using Default parameter.\n");
    param = 0;
  }
  
  genptr = (struct rngen *) mymalloc(1*sizeof(struct rngen)); 
  if(genptr == NULL)	 /* check if memory allocated for data structure */
    return NULL;
  
  /* Initiallize data structure variables */
  genptr->gentype = GENTYPE;
  genptr->stream_number = gennum;
  genptr->nstreams = total_gen;
  genptr->init_seed = seed & 0x7fffffff;  /* Only 31 LSB of seed considered */
  genptr->parameter = param;
  genptr->spawn_offset = total_gen;
  
  /*** Change the next line depending on your generators data needs ***/
  genptr->narrays = 0;		/* number of arrays needed by your generator */

  if(genptr->narrays > 0)
  {
    genptr->array_sizes = (int *) mymalloc(genptr->narrays*sizeof(int));
    genptr->arrays = (int **) mymalloc(genptr->narrays*sizeof(int *));
    if(genptr->array_sizes == NULL || genptr->arrays == NULL)
      return NULL;
  }
  else
  {
    genptr->array_sizes = NULL;
    genptr->arrays = NULL;
  }
  
  /*** Change the next line depending on your generators data needs ***/
  		/* initiallize ...array_sizes to the sizes of the arrays */


  for(i=0; i<genptr->narrays; i++)
  {
    genptr->arrays[i] = (int *) mymalloc(genptr->array_sizes[i]*sizeof(int)); 
    if(genptr->arrays[i] == NULL)  /* check if memory allocated for data structure */
      return NULL;
  }
  
  /*** Add initialization statements for your data in the arrays and other 
    variables you have defined ***/
  genptr->x = 0.1;  /* Please remove this statement in your implementation */
  
  NGENS++;			/* NGENS = # of streams */
  
  return (int *) genptr;
} 



/* Returns a double precision random number */

#ifdef __STDC__
double get_rn_dbl(int *igenptr)
#else
double get_rn_dbl(igenptr)
int *igenptr;
#endif
{
    struct rngen *genptr = (struct rngen *) igenptr;

    /*** Use the data in the structure genptr to update your state.
      Replace the statements below by those that return a double precision
      number in [0,1).  ***/

    if(genptr->x>0.85)	        /* the statements below are here just for    */
      genptr->x = 0.0;	/* testing purposes. Please _do_  remove */
    else                /* them in your implementation.          */ 
      genptr->x += 0.1;
    
    return genptr->x;			
} 



/* Return a random integer */

#ifdef __STDC__
int get_rn_int(int *igenptr)
#else
int get_rn_int(igenptr)
int *igenptr;
#endif
{
  /* If you have a more efficient way of computing the random integer in
     [0,2^31), then please replace the statement below with your scheme. */

  return (int) (get_rn_dbl(igenptr)*0x80000000);
} 



/* Return a single precision random number */

#ifdef __STDC__
float get_rn_flt(int *igenptr)
#else
float get_rn_flt(igenptr)
int *igenptr;
#endif
{
  /* If you have a more efficient way of computing the random integer,
     then please replace the statement below with your scheme.        */

    return (float) get_rn_dbl(igenptr);
}






/*************************************************************************/
/*************************************************************************/
/*                  SPAWN_RNG: spawns new generators                     */
/*************************************************************************/
/*************************************************************************/

#ifdef __STDC__
int spawn_rng(int *igenptr, int nspawned, int ***newgens, int checkid)
#else
int spawn_rng(igenptr,nspawned, newgens, checkid)
int *igenptr,nspawned, ***newgens, checkid;
#endif
{
  struct rngen **genptr, *tempptr = (struct rngen *) igenptr;
  int i, j;
  
  if (nspawned <= 0) /* is nspawned valid ? */
  {
    nspawned = 1;
    fprintf(stderr,"WARNING - spawn_rng: nspawned <= 0. Default value of 1 used for nspawned\n");
  }
  
  genptr = (struct rngen **) mymalloc(nspawned*sizeof(struct rngen *));
  if(genptr == NULL)	   /* allocate memory for pointers to structures */
  {
    *newgens = NULL;
    return 0;
  }
  
  for(i=0; i<nspawned; i++)	/* create nspawned new streams */
  {
    int seed, gennum;
    
    gennum = tempptr->stream_number + tempptr->spawn_offset*(i+1);
  
    if(gennum > MAX_STREAMS)   /* change seed to avoid repeating sequence */
      seed = (tempptr->init_seed)^gennum; 
     else
      seed = tempptr->init_seed;
   
    /* Initialize a stream. This stream has incorrect spawning information.
       But we will correct it below. */

  
    genptr[i] = (struct rngen *) 
      init_rng(gennum, gennum+1, seed, tempptr->parameter);
    
  
    if(genptr[i] == NULL)	/* Was generator initiallized? */
    {
      nspawned = i;
      break;
    }
    genptr[i]->spawn_offset = (nspawned+1)*tempptr->spawn_offset;
  }
  
  tempptr->spawn_offset *= (nspawned+1);
  

  *newgens = (int **) genptr;
  
  
  if(checkid != 0)
    for(i=0; i<nspawned; i++)
      if(addID(( int *) genptr[i]) == NULL)
	return i;
  
  return nspawned;
}


/* Free memory allocated for data structure associated with stream */

#ifdef __STDC__
int free_rng(int *genptr)
#else
int free_rng(genptr)
int *genptr;
#endif
{
  struct rngen *q;
  int i;
  
  q = (struct rngen *) genptr;
  assert(q != NULL);
  
  for(i=0; i<q->narrays; i++)
    free(q->arrays[i]);

  if(q->narrays > 0)
  {
    free(q->array_sizes);
    free(q->arrays);
  }
  
  free(q);

  NGENS--;
  return NGENS;
}


#ifdef __STDC__
int pack_rng( int *genptr, char **buffer)
#else
int pack_rng(genptr,buffer)
int *genptr;
char **buffer;
#endif
{
  char *temp_buffer;
  int size, i;
  struct rngen *q;
  int pos;

  q = (struct rngen *) genptr;
  size = sizeof(struct rngen) + q->narrays*sizeof(int) + strlen(q->gentype)+1;
  for(i=0; i<q->narrays; i++)
    size += q->array_sizes[i]*sizeof(int);
  
  temp_buffer = (char *) mymalloc(size); /* allocate memory */
  if(temp_buffer == NULL)
  {
    *buffer = NULL;
    return 0;
  }
  
  memcpy(temp_buffer,q,sizeof(struct rngen));
  pos = sizeof(struct rngen);
  strcpy(temp_buffer+pos,q->gentype);
  pos += strlen(q->gentype)+1;
  
  if(q->narrays > 0)
  {
    memcpy(temp_buffer+pos,q->array_sizes,q->narrays*sizeof(int));
    pos += q->narrays*sizeof(int);
    for(i=0; i<q->narrays; i++)
    {
      memcpy(temp_buffer+pos,q->arrays[i],q->array_sizes[i]*sizeof(int));
      pos += q->array_sizes[i]*sizeof(int);
    }
  }
  
  assert(pos == size);
  
  *buffer = temp_buffer;
  return size;
}



#ifdef __STDC__
int *unpack_rng( char *packed)
#else
int *unpack_rng(packed)
char *packed;
#endif
{
  struct rngen *q;
  int i;
  int pos;

  q = (struct rngen *) mymalloc(sizeof(struct rngen));
  if(q == NULL)
    return NULL;

  memcpy(q,packed,sizeof(struct rngen));
  pos = sizeof(struct rngen);

  if(strcmp(packed+pos,GENTYPE) != 0)
  {
    fprintf(stderr,"ERROR: Unpacked ' %.24s ' instead of ' %s '\n",  
	    packed+pos, GENTYPE); 
    return NULL; 
  }
  else
    q->gentype = GENTYPE;
  pos += strlen(q->gentype)+1;
    
  if(q->narrays > 0)
  {
    q->array_sizes = (int *) mymalloc(q->narrays*sizeof(int));
    q->arrays = (int **) mymalloc(q->narrays*sizeof(int *));
    if(q->array_sizes == NULL || q->arrays == NULL)
      return NULL;
    memcpy(q->array_sizes,packed+pos,q->narrays*sizeof(int));
    pos += q->narrays*sizeof(int);
  
    for(i=0; i<q->narrays; i++)
    {
      q->arrays[i] = (int *) mymalloc(q->array_sizes[i]*sizeof(int));
      if(q->arrays[i] == NULL)
	return NULL;
      memcpy(q->arrays[i],packed+pos,q->array_sizes[i]*sizeof(int));
      pos += q->array_sizes[i]*sizeof(int);
    }   
  }
  else				/* narrays == 0 */
  {
    q->array_sizes = NULL;
    q->arrays = NULL;
  }
    
  NGENS++;
  
  return (int *) q;
}

      

#ifdef __STDC__
int get_seed_rng(int *gen)
#else
int get_seed_rng(gen)
int *gen;
#endif
{
  return ((struct rngen *) gen)->init_seed;
}



#ifdef __STDC__
int print_rng( int *igen)
#else
int print_rng(igen)
int *igen;
#endif
{
  struct rngen *gen;
  
  printf("\n%s\n", GENTYPE);
  
  gen = (struct rngen *) igen;
  printf("\n \tseed = %d, stream_number = %d\tparameter = %d\n\n", gen->init_seed, gen->stream_number, gen->parameter);

  return 1;
}


#include "../simple_.h"
#include "../fwrap_.h"
