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
#define SVD_PAIRLIST_INCREMENTAL   8 //means we only updt prlst if atm mvd
#define SVD_USE_PAIRLIST          16

#define pairsRowsize 152
#define pairsMax 150

struct SimData {
  double *q;
  double *qold;
  double *force;
  double *boxsize;
  int flags; // not used yet... rework flag system if I start using it.
  int N;
  int Nmax;
  int ndim;
  double dt;
  double beta;
  long *atomtypes;
  //double *ei;
  long *pairlist;
  double trialMoveScale;
  double trialMoveIsobaricScale;
  double pairlist_minDistance;
  double prob_PMove;
  double isobaricPressure;

  int ntry_shift;
  int naccept_shift;
  int ntry_shift_last;
  int naccept_shift_last;
  int ntry_isobaric;
  int naccept_isobaric;
  int ntry_isobaric_last;
  int naccept_isobaric_last;
} ;

int init_mt(int seed)
{
  init_gen_rand(seed);
  return(0);
}

// #use -D FFid=X on the command line.

// Note: fij returns the force pushing particles _away_ from each other.
//       it should basically return -(dV/dr)
#if FFid == 01
#include "ff/01_hardsphere.c"

#elif FFid == 02
#include "ff/02_lennardjones.c"

#elif FFid == 03
#include "ff/03_harmonic.c"

#elif FFid == 10
#include "ff/10_2scale-s2s.c"

#elif FFid == 11
#include "ff/11_3scale-s3s.c"

#elif FFid == 12
#include "ff/12_3scale-s2s.c"

#elif FFid == 13
#include "ff/13_2scale-jagla98.c"

/* #elif FFid == 2 // this is 3 scale in uniform units */
/* #error not done yet */

#else
#error No Forcefield selected!
#endif

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
  int usePairlist = flags & SVD_USE_PAIRLIST;
  partial = flags & SVD_ENERGYI_PARTIAL; // defaults to 
  if (partial) j=i+1;


  if (usePairlist) {
    //printf("using pairlist in energy_i\n");
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
      //printf("d=%f\n", d);
      E += eij(atomtypes[i], atomtypes[j], d);
      if (E == 1./0.) return(E);
    }
  }
  //printf("E = %f\n", E);
  return(E) ;
}
double energy(struct SimData *SD,
	      int flags) {
  double E=0;
  int i;
  flags |= SVD_ENERGYI_PARTIAL ;

  for(i=0;  i < SD->N;  i++) {
    E += energy_i(SD, i, flags);
  }
  return (E);
}


inline void force_ij(int i, int j,    // atomtypes
	     double d,        // distance between atoms
	     double *r,       // x, y, z coordinates
	     double *force) { // forces placed in here
  double F;
  F = fij(i, j, d);

  force[0] = F * r[0]/d;
  force[1] = F * r[1]/d;
  force[2] = F * r[2]/d;
  
}
/* double forcedotr_i_sub(double qi0, double qi1, double qi2, */
/* 			      int j, */
/* 			      double *q, */
/* 			      double *boxsize */
/* 			      ) { */
/*   double drx, dry, drz;   // vectors FROM i TO j, makes force FROM i ON j */
/*   double d; */
/*   d = 0; */
  
/*   drx = q[j*3  ] - qi0 ; */
/*   drx -= floor((drx/boxsize[0] + .5) ) * boxsize[0]; */
/*   d += drx*drx; */
/*   //if (d>1) return(0);  // greater than cutoff^2 */

/*   dry = q[j*3+1] - qi1 ; */
/*   dry -= floor((dry/boxsize[1] + .5) ) * boxsize[1]; */
/*   d += dry*dry; */
/*   //if (d>1) return(0);  // greater than cutoff^2 */

/*   drz = q[j*3+2] - qi2 ; */
/*   drz -= floor((drz/boxsize[2] + .5) ) * boxsize[2]; */
/*   d += drz*drz; */
/*   //if (d>1) return(0);  // greater than cutoff^2 */

/*   d = sqrt(d); */
/*   return(d); */
/* } */
double forcedotr_i(struct SimData *SD,
		   int i,
		   int flags
		   ) {
  // defaults to using all of them
  double fdotr=0;
  int j=0;
  int partial;
  double *q = SD->q;
  double *boxsize = SD->boxsize;
  
  int usePairlist = flags & SVD_USE_PAIRLIST;
  partial = flags & SVD_ENERGYI_PARTIAL; // defaults to all
  if (partial) j=i+1;

  //double qi0 = SD->q[i*3  ];
  //double qi1 = SD->q[i*3+1];
  //double qi2 = SD->q[i*3+2].;
  //printf("%f %f %f\n", qi0, qi1, qi2);
  if (usePairlist) {
    // pairlist
    int ni = SD->pairlist[i*pairsRowsize];
    for (ni -= 1;  ni>=0;  ni-- ) {
      j = SD->pairlist[i*pairsRowsize + 2 + ni];
      if (SVD_ENERGYI_PARTIAL && j <= i)
	continue;
      double d = distance(q, boxsize, i, j);
      double F = fij(SD->atomtypes[i], SD->atomtypes[j], d);
      fdotr += F * d ;
    }
  } else {
    // NO pairlist
    for (; j<SD->N ; j++ ) {
      if (i == j) continue;
      //double d = forcedotr_i_sub(qi0, qi1, qi2, j, SD->q, SD->boxsize);
      double d = distance(q, boxsize, i, j);
      double F = fij(SD->atomtypes[i], SD->atomtypes[j], d);
      fdotr += F * d ;
    }
  }
  return(fdotr);
}
double forcedotr_total(struct SimData *SD,
		       int flags
		       ) {
  /*
   * for i in range(0, self.N-1):
   *     fdotr += forcedotr_i(self.SD_p, i, flags)
   */
  int i;
  double fdotr=0.;

  flags |= SVD_ENERGYI_PARTIAL;

  // This is supposed to be N-1, since the last atom is N.
  for(i=0 ; i < (SD->N)-1 ; i++) {
    fdotr += forcedotr_i(SD, i, flags);
  }
  return(fdotr);
}




int nCloserThan(struct SimData *SD,
		double maxdist,
		int flags
		) {
  /* This function returns the number of pairs of atoms closer than
   * the distance `maxdist`.  It is designed to be used to calculate
   * impulsive corrections to the pressure.
   *
   * Things to do:
   * - make it away of different atom types
   * - make it take the paramaters r, dr, and calculate the histogram
   *   bins [r, r+dr), [r, r+2dr), [r, r+3dr).  Then do a simple linear
   *   regression to find the slope.  (maybe make it two points so that
   *   it's easier to calculate the slope)
   */
  // defaults to using all of them
  int n=0;
  int i=0, j=0;
  double *q = SD->q;
  double *boxsize = SD->boxsize;
  int usePairlist = flags & SVD_USE_PAIRLIST;
  
  for(i=0 ; i < (SD->N)-1 ; i++) {
    j=i+1;
  
    if (usePairlist) {
      // pairlist
      int ni = SD->pairlist[i*pairsRowsize];
      for (ni -= 1;  ni>=0;  ni-- ) {
        j = SD->pairlist[i*pairsRowsize + 2 + ni];
        if (j <= i)
  	continue;
        double d = distance(q, boxsize, i, j);
        if (d < maxdist) {
	  n++;
        }
      }
    } else {
      // NO pairlist
      for (; j<SD->N ; j++ ) {
        if (i == j) continue;
        //double d = forcedotr_i_sub(qi0, qi1, qi2, j, SD->q, SD->boxsize);
        double d = distance(q, boxsize, i, j);
        if (d < maxdist) {
	  n++;
        }
      }
    }
  }
  return(n);
}




double calcForce(struct SimData *SD, int flags) {
  // defaults to using all of them
  int i, j;

  double *F = SD->force;
  memset(F, 0, SD->N*3*sizeof(double));

  //printf("%f %f %f\n", qi0, qi1, qi2);
  for (i=0; i<(SD->N-1) ; i++ ) {

    double qi0 = SD->q[i*3  ];
    double qi1 = SD->q[i*3+1];
    double qi2 = SD->q[i*3+2];
    for (j=i+1; j<SD->N ; j++ ) {
      //printf("%d %d\n", i, j);
      double drx, dry, drz;   // vectors FROM i TO j, makes force FROM i ON j
      double d;
      if (i == j) printf("we should'nt be here raockhorke\n");
      d = 0;
  
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
      if (d < .25) {
	printf("%f \n", d);
      }
  
      double Fk[3], rk[3];
      double Fij;
      rk[0] = drx;
      rk[1] = dry;
      rk[2] = drz;
      /*force_ij(SD->atomtypes[i], SD->atomtypes[j], d, 
	       rk,
	       Fk);*/
      Fij = fij(SD->atomtypes[i], SD->atomtypes[j], d);
      Fk[0] = Fij * rk[0]/d;
      Fk[1] = Fij * rk[1]/d;
      Fk[2] = Fij * rk[2]/d;
      /*if (Fij > 10. ) {
	printf("%d, %d, %f, %f, %f %f %f \n", 
	         i,  j,  d, Fij, Fk[0], Fk[1], Fk[2]);
		 }*/

      F[i*3  ] -= Fk[0];
      F[i*3+1] -= Fk[1];
      F[i*3+2] -= Fk[2];
      F[j*3  ] += Fk[0];
      F[j*3+1] += Fk[1];
      F[j*3+2] += Fk[2];
      //F = fij(SD->atomtypes[i], SD->atomtypes[j], d);
      
      
      //fdotr += drx*(F*drx/d)  +  dry*(F*dry/d)  +  drz*(F*drz/d);
    }
  }
  
  //printf("E = %f\n", E);
  return(0);
}

double integrate(struct SimData *SD, int flags) {
  // integrate to the next time step, return the KE.
  double *q, *qold, *force;

  q = SD->q ;
  qold = SD->qold;
  force = SD->force;
  double mass = 1.;
  double qnew;
  double dt = SD->dt;
  double dt2 = dt*dt;  // delta time squared
  double v, KE = 0;

  int i;
  for(i=0;  i < (SD->N)*3;  i++) {
/*     v2 = 0; */
/*     for (k=0; k<3; k++) { */
      qnew = 2*q[i] - qold[i] + force[i]*dt2/mass;
      //qnew = 2* (*q)     - (*qold) + (*force)*dt2/mass;

      v = (qnew - qold[i])/(2*dt);
      //v = (qnew - (*qold))/(2*dt);

      //v = (qnew - q[i*3+k])/(dt);
      //v = (qnew - (*q))/(dt);
      //v2 += v * v ;

      KE += .5 * mass * v * v;
      
      qold[i] = q[i];
      q[i] = qnew;
      //(*qold) = (*q);
      //(*q) = qnew;
      //q++; qold++; force++;
/*     /} */
    
/*     KE += .5 * mass * v2; */
  }
  return(KE);
}

double mdStep(struct SimData *SD, int n, int flags) {
  int count ;
  double KE = 0 ;
  for (count=0; count<n; count++) {
    calcForce(SD, flags);
/*     printf("x1=%f x2=%f f1=%f f2=%f u1=%f u2=%f\n",  */
/* 	   SD->q[3*0+0],     SD->q[3*1+0], */
/* 	   SD->force[3*0+0], SD->force[3*1+0],  */
/* 	   energy_i(SD, 0, flags), energy_i(SD, 0, flags)); */
    KE = integrate(SD, flags);
/*     printf("\n"); */
  }  
  return(KE);
}

double kineticEnergy(struct SimData *SD) {
  double KE = 0;
  
  int i;
  for (i=0;  i < (SD->N)*3;  i++) {
    KE += 1;
  }
  return(KE);
}

/* double forcedotr_term(double *qi, int i, int j, int *atomtypes) { */
/* } */

double pressure_c(struct SimData *SD, int flags) {
  flags |= SVD_ENERGYI_PARTIAL;
  //  double fdotr = 0;

  
  return(-1);
}


int trialMove_isobaric(struct SimData *SD, int flags);
int trialMove(struct SimData *SD, int n, int flags) {
  /* Run N monte carlo trial moves, return number of moves accepted */
  double qi_old[3];
  double Eold, Enew;
  int ntry=0, naccept=0, i, counter1, accept;
  double x, ran;
  double trialmovescale = SD->trialMoveScale;
  for(counter1=0; counter1<n; counter1++) {

    /* What type of move do we need to do? */
    ran = genrand_real2(); /* on [0,1) */
    if (ran < SD->prob_PMove) {
      trialMove_isobaric(SD, flags);
      continue;
    }
    
    i = (int) floor(genrand_real2()*(SD->N));
    qi_old[0] = SD->q[i*3  ];
    qi_old[1] = SD->q[i*3+1];
    qi_old[2] = SD->q[i*3+2];
    //Eold = SD->ei[i];   // ei[i] not updated if something else moves closer
    Eold = energy_i(SD, i, flags);

    // v randn returns normal gaussion distributed points
    //#define trialmovescale .15
    SD->q[i*3  ] += trialmovescale*(genrand_real1() - .5) ;
    SD->q[i*3+1] += trialmovescale*(genrand_real1() - .5) ;
    SD->q[i*3+2] += trialmovescale*(genrand_real1() - .5) ;

    Enew = energy_i(SD, i, flags);
    if (Enew <= Eold)     /* always accept, E decreases */
      accept = 1;
    else if (Enew == 1./0.) {
      accept = 0;
    }
    else {
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
      //SD->ei[i] = Enew;
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
  SD->ntry_shift += ntry;
  SD->naccept_shift += naccept;
  return(naccept);
}
int trialMove_isobaric(struct SimData *SD, int flags) {
  int i;
  double *q = SD->q;
  double *boxsize = SD->boxsize;
  double ran;

  //def trialMove_isobaric_py(self, pressure, lnVScale):
  double lnVScale = SD->trialMoveIsobaricScale;
  double pressure = SD->isobaricPressure;

  double Vold = boxsize[0] * boxsize[1] * boxsize[2];
  double Eold = energy(SD, flags);

  //#Vnew = exp(ln(V + (random.random()-.5) * scale))
  //#lengthscale = (Vnew / Vold)**(1./3.) # this may not be right
  double linearScale = exp(((genrand_real1()/*[0,1]*/-.5) * lnVScale) / 3.);
  double volumeScale = linearScale * linearScale * linearScale;


  //self.q *= linearScale
  for (i=0 ; i < SD->N*3; i++) {
    q[i] *= linearScale;
  }
  //self.boxsize *= linearScale
  boxsize[0]*=linearScale; boxsize[1]*=linearScale; boxsize[2]*=linearScale;

  double Enew = energy(SD, flags);
  double Vnew = Vold * volumeScale;


  double x = - SD->beta * (Enew - Eold + pressure*(Vnew-Vold) - \
			   ((SD->N+1)/SD->beta)*log(volumeScale));
  //printf("%f\n", x);

  int accept;
  if (x > 0)
    accept = 1;
  else {
    x = exp(x);
    //printf("%f\n", x);
    ran = genrand_real2();
    if (ran < x)  accept = 1;
    else          accept = 0;
  }
  SD->ntry_isobaric += 1;
  if (accept) {
    //printf("+++ %.3f pressure move  +++\n", linearScale); 
    SD->naccept_isobaric += 1;
  }
  else {
    //printf("--- %.3f pressure move  ---\n", linearScale);
    //self.q /= linearScale;
    for (i=0 ; i < SD->N*3; i++) {
      q[i] /= linearScale;
    }
    //self.boxsize /= linearScale;
    boxsize[0]/=linearScale; boxsize[1]/=linearScale; boxsize[2]/=linearScale;
  }
  //if (flags & SVD_VERBOSE_1 || 1) {
  //    printf("Isobaric move accept=%d, vs=%f, x=%f, Eold/new=%.3f %.3f\n", accept, volumeScale, x, Eold, Enew);
  //}
  return(0);
}


int pairlistInit(struct SimData *SD, double cutoff, int flags) {
  // SD->pairlist is a (natoms, pairsRowsize) array, capable of holding up to 15 pairs
  // - The first element in the row is for number of pairs
  // - second element is for number of moves that atom has done since
  //   the pairlist was created (not used yet)
  // - remainder of the rows are for the pairlist pairs

  int i, j; // should be long ??
  long *pairlist=SD->pairlist;

  int incremental = flags&SVD_PAIRLIST_INCREMENTAL;

  if (! incremental ) {
    memset(pairlist, '\0', sizeof(long)*pairsRowsize*(SD->Nmax));
  }

  for (i=0; i<SD->N; i++) {
    //printf("i_pl_i: %d\n", i);
    if ( incremental ) {
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
	  printf("Pairlist overfull for %d.\n", i);
	  exit(5);
	}
	int nj = pairlist[j*pairsRowsize];
	if (nj == pairsMax) {  
	  printf("Pairlist overfull for %d.\n", j);
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

int pairlistCheck(struct SimData *SD, double warn, int flags) {
  // SD->pairlist is a (natoms, pairsRowsize) array, capable of holding up to 15 pairs
  // - The first element in the row is for number of pairs
  // - second element is for number of moves that atom has done since
  //   the pairlist was created (not used yet)
  long *pairlist=SD->pairlist;
  int i, j, count1;
  int nviolations=0;
  double minDistance = 1e6 ;
  
  for (i=0; i<SD->N; i++) {
    j=i+1;
    for (; j<SD->N; j++) {
      double d;
      if (i==j) continue;
      d = distance(SD->q, SD->boxsize, i, j);
      //printf("d: %f\n", d);
      if (d <= warn) {
	int foundInPairlist;

	// search for i in j's list:
	int ni = pairlist[i*pairsRowsize];

	foundInPairlist = 0;
	for (count1=0; count1<ni ; count1++) {
	  if (pairlist[i*pairsRowsize + 2 + count1] == j)
	    {
	      foundInPairlist = 1;
	      break;
	    }
	}
	if (! foundInPairlist) {
	  printf("Atom %d should be in atom %d's pairlist, but isn't (d=%f).\n",
		 j, i, d);
	  nviolations++ ;
	  // We want to record the closest distance of all atoms that were
	  // not in the original pairlist.
	  if (d < minDistance) { // store the closest approach for an 
	    minDistance = d;     //   adaptive pairlist cutoff algorithm
	  }
	}
      }
    }
  }
  SD->pairlist_minDistance = minDistance;
  return(nviolations);
}
