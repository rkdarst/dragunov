/* Richard Darst, 2007 */

/* This is a simple inverse square particle.
 *  V = 1/r
 *  F = 1/r^2
 */

#define HO_CUTOFF 10.

inline double eij(int i, int j, double d) {

  double E;

  // d is the center to center distance.  

  // Note: dragunov, by definition, returns the force pushing
  // particles _away_ from ach other.
  if(d < HO_CUTOFF) {
    E = .5 *  d*d;
    return (E);
  }
  else
    return(0.0);
}



inline double fij(int i, int j, double d) {
  double F;

  if(d < HO_CUTOFF) {
    F = - d ;
    return (F);
  }
  else
    return(0.0);
}
