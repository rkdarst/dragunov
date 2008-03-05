/* Richard Darst, 2007 */


//#warning Note: truncated Lennard Jones potential
#ifndef LJ_SIGMA
  #define LJ_SIGMA 1.0
  //#warning Using lj_sigma = 1.0
#endif

#ifndef LJ_EPSILON
  #define LJ_EPSILON 1.0
  //#warning Using lj_epsilon = 1.0
#endif

#define LJ_CUTOFF 2.5

inline double eij(int i, int j, double d) {
  //int x;
  //if (j<i) { x=i; i=j; j=x; } // make i the lower index

  float E;
  float d_orig = d;

  // d is the center to center distance.  

  // Note: dragunov, by definition, returns the force pushing
  // particles _away_ from ach other.
  if(d_orig < LJ_CUTOFF*LJ_SIGMA) {
    d = LJ_SIGMA / d ;
    float d6 = d*d*d;
    d6 *= d6;
    E = 4 * (d6*d6 - d6);
    E *= LJ_EPSILON;
    //E += 0.016316891136;        // U(r=2.5 sigma)
    //E += 0.00217478039165499;   // U(r=3.5 sigma)
    return (E);
  }
  else
    return(0.0);
}



inline double fij(int i, int j, double d) {
  //int x;
  //if (j<i) { x=i; i=j; j=x; } // make i the lower index
  float F;
  float d_orig = d;


  if(d_orig < LJ_CUTOFF*LJ_SIGMA) {
    d = LJ_SIGMA / d ;    float d6 = d*d*d;
    d6 *= d6;

    F = 24 * (LJ_EPSILON / d_orig) * (2 * d6*d6 - d6);
    return (F);
  }
  else
    return(0.0);
}
