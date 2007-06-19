#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK printf(".\n");

#define energy_i energy_i_4c

//#include "SFMT-params.h"
#include "SFMT.h"
//#include "SFMT.c"
#define random genrand_real1

struct SimData {
  double *q;
  double *boxsize;
  int flags;
  int N;
  int Nmax;
  int ndim;
  double beta;
  int *atomtypes;
  double *ei;
} ;

int init_mt(int seed)
{
  init_gen_rand(seed);
  return(0);
}

inline double eij(int i, int j, double d)
{
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index
  if(d < 1.)
    return 1./0.;
  else if (d < .7)
    return (1.5 - d);
  else
    return(0.0);
}

#define SVD_ENERGYI_PARTIAL 1
#define SVD_VERBOSE_1       2
#define SVD_VERBOSE_2       4

double energy_i_4c(struct SimData *SD,
		   int i,
		   int flags
		   )
{ // defaults to using all of them
    double E=0;
    double d;
    int j=0;
    int partial;
    double x;

    partial = flags & SVD_ENERGYI_PARTIAL; // defaults to 
    if (partial) j=i+1;

    double qi0 = SD->q[i*3  ];
    double qi1 = SD->q[i*3+1];
    double qi2 = SD->q[i*3+2];
    //printf("%f %f %f\n", qi0, qi1, qi2);
    for (; j<SD->N ; j++ )
      {
	//printf("%d %d\n", i, j);

	if (i == j) continue;
	d = 0;

	x = SD->q[j*3  ] - qi0 ;
	x -= floor((x/SD->boxsize[0] + .5) ) * SD->boxsize[0];
	d += x*x;

	x = SD->q[j*3+1] - qi1 ;
	x -= floor((x/SD->boxsize[1] + .5) ) * SD->boxsize[1];
	d += x*x;

	x = SD->q[j*3+2] - qi2 ;
	x -= floor((x/SD->boxsize[2] + .5) ) * SD->boxsize[2];
	d += x*x;

	d = sqrt(d);

	/*if (flags&0x2 && ( d<1 || d>17.5))
	  {
	    CHECK;
	    printf("i = %d\n", i);
	    printf("j = %d\n", j);
	    printf("i = %f %f %f\n", qi0, qi1, qi2);
	    printf("j = %f %f %f\n", SD->q[j*3+0], SD->q[j*3+1], SD->q[j*3+2]);
	    printf("%f %f\n", d, eij(0, 0, d));
	    //sleep(5);
	    exit(4);
	    }*/

	E += eij(0, 0, d);
	if (E == 1./0.) return(E);
      }
    //printf("E = %f\n", E);
    return(E) ;
  }

int trialMove(struct SimData *SD, int n)
     /* Run N monte carlo trial moves, return number of moves accepted */
{
  double qi_old[3];
  double Eold, Enew;
  int ntry=0, naccept=0, i, counter1;
  for(counter1=0; counter1<n; counter1++)
  {
  
  i = (int) floor(genrand_real2()*(SD->N));
  //i = (int) floor((((double)random())/(RAND_MAX))*(SD->N) );
  //printf("  i:  %d \n", i);
  qi_old[0] = SD->q[i*3  ];
  qi_old[1] = SD->q[i*3+1];
  qi_old[2] = SD->q[i*3+2];
  Eold = SD->ei[i];

  // v randn returns normal gaussion distributed points
  SD->q[i*3  ] += genrand_real1() - .5 ;
  SD->q[i*3+1] += genrand_real1() - .5 ;
  SD->q[i*3+2] += genrand_real1() - .5 ;
  /*SD->q[i*3  ] += ((double)random()/RAND_MAX)-.5 ;
  SD->q[i*3+1] += ((double)random()/RAND_MAX)-.5 ;
  SD->q[i*3+2] += ((double)random()/RAND_MAX)-.5 ;*/

  Enew = energy_i(SD, i, 0);
  int accept;
  if (Enew <= Eold)     /* always accept, E decreases */
    accept = 1;
  else
    {
      double x, ran;
      x = exp(SD->beta*(Eold-Enew));
      ran = genrand_real2();
      //ran / (double)random()/RAND_MAX; 
      //printf("  ran, x: %f %f\n", ran, x);
      if (ran < x)
	accept = 1;
      else
	accept = 0;
    }
  ntry += 1;
  if(accept)
    {
      SD->ei[i] = Enew;
      naccept += 1;
    }
  else
    {
      SD->q[i*3  ] = qi_old[0];
      SD->q[i*3+1] = qi_old[1];
      SD->q[i*3+2] = qi_old[2];
    }
  /*if (0)
    printf("%6d %d, %5.3f->%5.3f accF: %.4f\n", 
    counter1, accept, Eold, Enew, (double)naccept/ntry);*/
  }
  return(naccept);
}

