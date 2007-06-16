#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK printf(".\n");

struct SimData {
  double *q;
  double *boxsize;
  int flags;
  int N;
  int Nmax;
  int ndim;
  double beta;
} ;

inline double eij(int i, int j, double d)
{
  if(d < 1.)
    return 1./0.;
  else if (d < .7)
    return (1.5 - d);
  else
    return(0.0);
}

#define dist(a) \
{ \
	x = qdata[a*3  ] - qi0 ; \
	x -= floor((x/boxsize[0] + .5) ) * boxsize[0]; \
	d += x*x; \
\
	x = qdata[a*3+1] - qi1 ; \
	x -= floor((x/boxsize[1] + .5) ) * boxsize[1]; \
	d += x*x; \
\
	x = qdata[a*3+2] - qi2 ; \
	x -= floor((x/boxsize[2] + .5) ) * boxsize[2]; \
	d += x*x; \
\
	d = sqrt(d);\
}


double energy_i_4c(struct SimData *SD,
		   int i,
		   int flags
		   )
  {
    double E=0;
    double d;
    int j=0;
    int all;
    double x;

    all = flags & 0x1;
    if (!all) j=i+1;

#define qdata SD->q
#define N SD->N
#define boxsize SD->boxsize

    double qi0 = qdata[i*3  ];
    double qi1 = qdata[i*3+1];
    double qi2 = qdata[i*3+2];
    //printf("%f %f %f\n", qi0, qi1, qi2);
    for (; j<N ; j++ )
      {
	//printf("%d %d\n", i, j);

	if (i == j) continue;
	d = 0;

	x = qdata[j*3  ] - qi0 ;
	x -= floor((x/boxsize[0] + .5) ) * boxsize[0];
	d += x*x;

	x = qdata[j*3+1] - qi1 ;
	x -= floor((x/boxsize[1] + .5) ) * boxsize[1];
	d += x*x;

	x = qdata[j*3+2] - qi2 ;
	x -= floor((x/boxsize[2] + .5) ) * boxsize[2];
	d += x*x;

	d = sqrt(d);

	/*if (flags&0x2 && ( d<1 || d>17.5))
	  {
	    CHECK;
	    printf("i = %d\n", i);
	    printf("j = %d\n", j);
	    printf("i = %f %f %f\n", qi0, qi1, qi2);
	    printf("j = %f %f %f\n", qdata[j*3+0], qdata[j*3+1], qdata[j*3+2]);
	    printf("%f %f\n", d, eij(0, 0, d));
	    //sleep(5);
	    exit(4);
	    }*/

	E += eij(0, 0, d);
      }
    //printf("E = %f\n", E);
    return(E) ;
  }

/* double energy_i_4c2(double *qdata, int N, */
/* 		   int i, */
/* 		   double *boxsize, */
/* 		   int flags */
/* 		   ) */
/*   { */
/*     double E=0; */
/*     double d; */
/*     int j=0; */
/*     int all; */
/*     double x; */

/*     all = flags & 0x1; */
/*     if (!all) j=i+1; */

/*     double qi0 = qdata[i*3  ]; */
/*     double qi1 = qdata[i*3+1]; */
/*     double qi2 = qdata[i*3+2]; */
/*     //printf("%f %f %f\n", qi0, qi1, qi2); */
/*     for (; j<N ; j++ ) */
/*       { */
/* 	//printf("%d %d\n", i, j); */

/* 	if (i == j) continue; */
/* 	d = 0; */
/* 	dist(j); */

/* 	/\*if (flags&0x2 && ( d<1 || d>17.5)) */
/* 	  { */
/* 	    CHECK; */
/* 	    printf("i = %d\n", i); */
/* 	    printf("j = %d\n", j); */
/* 	    printf("i = %f %f %f\n", qi0, qi1, qi2); */
/* 	    printf("j = %f %f %f\n", qdata[j*3+0], qdata[j*3+1], qdata[j*3+2]); */
/* 	    printf("%f %f\n", d, eij(0, 0, d)); */
/* 	    //sleep(5); */
/* 	    exit(4); */
/* 	    }*\/ */

/* 	E += eij(0, 0, d); */
/*       } */
/*     //printf("E = %f\n", E); */
/*     return(E) ; */
/*  } */
