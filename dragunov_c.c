#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK printf(".\n");
#include <sys/types.h>

//#include "SFMT-params.h"
#include "SFMT.h"
//#include "SFMT.c"
#define random genrand_real1

#define SVD_ENERGYI_PARTIAL        1
#define SVD_VERBOSE_1              2
#define SVD_VERBOSE_2              4
#define SVD_PAIRLIST_INCREMENTAL   8

struct SimData {
  double *q;
  double *boxsize;
  int flags;
  int N;
  int Nmax;
  int ndim;
  double beta;
  long *atomtypes;
  double *ei;
  long *pairlist;
  double trialMoveScale;
} ;

int init_mt(int seed)
{
  init_gen_rand(seed);
  return(0);
}

// #define FFid 3

#if FFid == 1
#include "ff/3scale-s2s.c"
#error not done yet

#elif FFid == 2 // this is 3 scale in uniform units
#error not done yet

#elif FFid == 3
#include "ff/3scale-s2s.c"

#elif FFid == 4
#include "ff/3scale-s3s.c"

#else
#error No Forcefield selected!
#endif


/* inline double eij(int i, int j, double d) { */
/*   int x; */
/*   if (j<i) { x=i; i=j; j=x; } // make i the lower index */

/* #define FFid 3 */

/* /\* Common parameters *\/ */


/* #if FFid == 1 */
/* #warning Using Force Field 2-scale Stanly-1S units. */
/* #define lambda 0.571428571428571 // 4/7 */
/*   if (i==0 && j==0) // Jagla-Jagla */
/*     { */
/*       //printf("0, 0\n"); */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*       else if (d < 1.) */
/* 	return (1. - d); */
/*       /\*else if (d < 1.) */
/* 	return (.5 - d); */
/*       else if (d < 1.5) */
/*       return (d - 1.5);*\/ */
/*       else */
/* 	return(0.0); */
/*     } */
/*   else if (i==0 && j==1) */
/*     { */
/*       //printf("0, 1\n"); */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*     } */
/*   else if (i==1 && j==1) */
/*     { */
/*       //printf("1, 1\n"); */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*     } */
/*   else */
/*     { */
/*       printf("unknown atom types (%d, %d), quitting.\n", i, j); */
/*       exit(145); */
/*     } */
/* #elif FFid == 2 */
/* #warning Using Force Field 3-scale Uniform units */
/* #define lambda_S2S_1 0.581395348837209 */
/* #define lambda_S2S_2 0.333333333333333 */
/* #define inv_lambda_S2S_1 1.72 */
/* #define inv_lambda_S2S_2 3. */
/*   if (i==0 && j==0) // Jagla-Jagla */
/*     { */
/*       if(d < 1.) */
/* 	return (1./0.); */
/*       else if (d < inv_lambda_S2S_1) */
/* 	return (   1/(inv_lambda_S2S_1 -1.  )  *d    + */
/* 		3.5/(4.5*(1-lambda_S2S_1))); */

/*       else if (d < inv_lambda_S2S_2) */
/* 	return ( -1./(4.5*(1-lambda_S2S_2)* */
/* 		     (inv_lambda_S2S_1-inv_lambda_S2S_2)) */
/* 		 * ( d-3. )); */
/*       else */
/* 	return(0.0); */
/*     } */
/*   else if (i==0 && j==1) */
/*     { */
/*       //printf("0, 1\n"); */
/*       if(d < 1) */
/* 	return (1./0.); */
/*     } */
/*   else if (i==1 && j==1) */
/*     { */
/*       //printf("1, 1\n"); */
/*       if(d < 1) */
/* 	return (1./0.); */
/*     } */
/*   else */
/*     { */
/*       printf("unknown atom types (%d, %d), quitting.\n", i, j); */
/*       exit(145); */
/*     } */
/* #elif FFid == 3 */
/* #warning Using 3-scale model with S-2S units. */
/* #define lambda 0.581395348837209  // 1/1.72 */
/* // gnuplot> plot [0:1.74] (-1/(4.5*(1-3*l)) * x + 3*l/(4.5*(1-3*l))), (1-x)-1/4.5, -1/4.5, 0  */
/*   if (i==0 && j==0) // Jagla-Jagla */
/*     { */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*       else if (d < 1.) */
/* 	return (1. - d) - 1./(4.5); */
/*       else if (d < 3.*lambda) */
/* 	return (-1./(4.5*(1-3*lambda))*d + 3*lambda/(4.5*(1-3*lambda)) ); */
/*       else */
/* 	return(0.0); */
/*     } */
/*   else if (i==0 && j==1) */
/*     { */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*     } */
/*   else if (i==1 && j==1) */
/*     { */
/*       if(d < lambda) */
/* 	return (1./0.); */
/*     } */
/*   else */
/*     { */
/*       printf("unknown atom types (%d, %d), quitting.\n", i, j); */
/*       exit(145); */
/*     } */
/* #else */
/* #error No Forcefield defined! */
/* #endif */
/*   return(0.0); // we shouldn't get here. */
/* } */



inline double distance(double *q, double *boxsize, int i, int j) {
  double d=0;
  double x;

  x = q[j*3  ] - q[i*3  ] ;
  x -= floor((x/boxsize[0] + .5) ) * boxsize[0];
  d += x*x;

  x = q[j*3+1] - q[i*3+1] ;
  x -= floor((x/boxsize[1] + .5) ) * boxsize[1];
  d += x*x;

  x = q[j*3+2] - q[i*3+2] ;
  x -= floor((x/boxsize[2] + .5) ) * boxsize[2];
  d += x*x;
  return (sqrt(d));
}

#define pairsRowsize 102
#define pairsMax 100
double energy_i(struct SimData *SD,
		   int i,
		   int flags
		   ) {
  // defaults to using all of them
  double E=0;
  double d;
  int j=0;
  int partial;
  double *q, *boxsize;
  long *atomtypes = SD->atomtypes;
  q = SD->q;
  boxsize = SD->boxsize;
  int usePairlist=0;
  partial = flags & SVD_ENERGYI_PARTIAL; // defaults to 
  if (partial) j=i+1;


  if (usePairlist) {
    int ni = SD->pairlist[i*pairsRowsize];
    for (ni -= 1;  ni>=0;  ni-- ) {
      j = SD->pairlist[i*pairsRowsize + 2 + ni];
      d = distance(q, boxsize, i, j);
      E += eij(atomtypes[i], atomtypes[j], d);
      if (E == 1./0.) return(E);
    }
  }
  else {
    for (; j<SD->N ; j++ ) {
      if (i == j) continue;
      d = distance(q, boxsize, i, j);
      E += eij(atomtypes[i], atomtypes[j], d);
      if (E == 1./0.) return(E);
    }
  }
  //printf("E = %f\n", E);
  return(E) ;
}


double forcedotr_i(struct SimData *SD,
		   int i,
		   int flags
		   ) {
  // defaults to using all of them
  double fdotr=0;
  int j=0;
  int partial;
  
  partial = flags & SVD_ENERGYI_PARTIAL; // defaults to all
  if (partial) j=i+1;

  double qi0 = SD->q[i*3  ];
  double qi1 = SD->q[i*3+1];
  double qi2 = SD->q[i*3+2];
  //printf("%f %f %f\n", qi0, qi1, qi2);
  for (; j<SD->N ; j++ ) {
    //printf("%d %d\n", i, j);
    double drx, dry, drz;
    double d;
    if (i == j) continue;
    d = 0;
    double rij[3];

    drx = SD->q[j*3  ] - qi0 ;
    drx -= floor((drx/SD->boxsize[0] + .5) ) * SD->boxsize[0];
    d += drx*drx;
    //if (d>1) return(0);  // greater than cutoff^2

    dry = SD->q[j*3+1] - qi1 ;
    dry -= floor((dry/SD->boxsize[1] + .5) ) * SD->boxsize[1];
    d += dry*dry;
    //if (d>1) return(0);  // greater than cutoff^2

    drz = SD->q[j*3+2] - qi2 ;
    drz -= floor((drz/SD->boxsize[2] + .5) ) * SD->boxsize[2];
    d += drz*drz;
    //if (d>1) return(0);  // greater than cutoff^2

    d = sqrt(d);
    double force[3];
    rij[0] = drx;   rij[1] = dry;   rij[2] = drz;
    fij(SD->atomtypes[i], SD->atomtypes[j], d, rij, force);
    fdotr += drx*force[0] + dry*force[1] + drx*force[2];

    if (flags&SVD_VERBOSE_1) {
      CHECK;
    }

  }
  //printf("E = %f\n", E);
  return(fdotr);
}


int trialMove(struct SimData *SD, int n) {
  /* Run N monte carlo trial moves, return number of moves accepted */
  double qi_old[3];
  double Eold, Enew;
  int ntry=0, naccept=0, i, counter1, accept;
  for(counter1=0; counter1<n; counter1++) {
    
    i = (int) floor(genrand_real2()*(SD->N));
    qi_old[0] = SD->q[i*3  ];
    qi_old[1] = SD->q[i*3+1];
    qi_old[2] = SD->q[i*3+2];
    //Eold = SD->ei[i];   // ei[i] not updated if something else moves closer
    Eold = energy_i(SD, i, 0);

    // v randn returns normal gaussion distributed points
    //#define trialmovescale .15
    double trialmovescale = SD->trialMoveScale;
    SD->q[i*3  ] += trialmovescale*(genrand_real1() - .5) ;
    SD->q[i*3+1] += trialmovescale*(genrand_real1() - .5) ;
    SD->q[i*3+2] += trialmovescale*(genrand_real1() - .5) ;

    Enew = energy_i(SD, i, 0);
    if (Enew <= Eold)     /* always accept, E decreases */
      accept = 1;
    else if (Enew == 1./0.) {
      accept = 0;
    }
    else {
      double x, ran;
      x = exp(SD->beta*(Eold-Enew));
      ran = genrand_real2();
      if (ran < x)
	accept = 1;
      else
	accept = 0;
    }
    ntry += 1;
    //SD->pairlist[i*pairsRowsize+1] += 1;
    if(accept) {
      SD->ei[i] = Enew;
      naccept += 1;
      //SD->pairlist[i*pairsRowsize+1] += 100;
    }
    else {
      SD->q[i*3  ] = qi_old[0];
      SD->q[i*3+1] = qi_old[1];
      SD->q[i*3+2] = qi_old[2];
    }
    /*if (0)
      printf("%6d %d, %5.3f->%5.3f accF: %.4f\n", 
      counter1, accept, Eold, Enew, (double)naccept/ntry);*/
  } // end loop over n, 
  return(naccept);
}

int pairlist_init(struct SimData *SD, double cutoff, int flags) {
  // SD->pairlist is a (natoms, pairsRowsize) array, capable of holding up to 15 pairs
  // - The first element in the row is for number of pairs
  // - second element is for number of moves that atom has done since
  //   the pairlist was created (not used yet)
  long i, j;
  long *pairlist=SD->pairlist;

  if (! (flags&SVD_PAIRLIST_INCREMENTAL)) {
    memset(pairlist, '\0', sizeof(long)*pairsRowsize*(SD->Nmax));
  }

  for (i=0; i<SD->N; i++) {
    //printf("i_pl_i: %ld\n", i);
    if (flags&SVD_PAIRLIST_INCREMENTAL) {
      if (pairlist[i*pairsRowsize+1] == 0) {
	// atom hasn't moved, don't regenerate it
	continue;
      }
      else { 
	// Zero pairlist for atom i
	memset(pairlist+(i*pairsRowsize), '\0', sizeof(long)*pairsRowsize);
      }
      j = 0;
    }
    else {
    }
    j=i+1;
    for (; j<SD->N; j++) {
      double d;
      if (i==j) continue;
      d = distance(SD->q, SD->boxsize, i, j);
      //printf("d: %f\n", d);
      if (d <= cutoff) {
	//printf("Adding atom %3d to pairlist of atom %3d (distance %f).\n", 
	//       i, j, d);
	int ni = pairlist[i*pairsRowsize];
	if ( ni >= pairsMax) {  
	  printf("Pairlist overfull for %ld.\n", i);
	  exit(5);
	}
	int nj = pairlist[j*pairsRowsize];
	if (nj == pairsMax) {  
	  printf("Pairlist overfull for %ld.\n", j);
	  exit(5);
	}
		     
	pairlist[i*pairsRowsize + 2 + ni] = j ;
	pairlist[j*pairsRowsize + 2 + nj] = i ;
	pairlist[i*pairsRowsize] ++;
	pairlist[j*pairsRowsize] ++;
      }
    }
  }
  return(0);
}

int pairlist_check(struct SimData *SD, double warn, int flags) {
  // SD->pairlist is a (natoms, pairsRowsize) array, capable of holding up to 15 pairs
  // - The first element in the row is for number of pairs
  // - second element is for number of moves that atom has done since
  //   the pairlist was created (not used yet)
  long *pairlist=SD->pairlist;
  long int i, count1;
  int nviolations=0;
  
  for (i=0;  i<SD->N;  i++) {
    for (count1=0;  count1<pairlist[pairsRowsize*i];  count1++) {
      double d;
      d = distance(SD->q, SD->boxsize, 
		   i, pairlist[pairsRowsize*i + 2 + count1]);
      if (d>warn) {
	long int j = pairlist[pairsRowsize*i + 2 + count1];
	printf("  Pairlist violation atoms %ld(%ld moves), %ld(%ld moves)\n", 
	       i, pairlist[i*pairsRowsize+1], 
	       j, pairlist[j*pairsRowsize+1]);
	nviolations++;
      }
    }
  }
  return(nviolations);
}
