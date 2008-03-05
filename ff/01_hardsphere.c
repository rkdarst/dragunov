/* Richard Darst, 2007 */

#warning Note: Hard Sphere radius=.5 (diameter=1) potential

inline double eij(int i, int j, double d) {
  //int x;
  //if (j<i) { x=i; i=j; j=x; } // make i the lower index

  // d is the center to center distance.
  if(d < 1.)
    return (1./0.);
  else
    return(0.0);
}



inline double fij(int i, int j, double d) {
  //int x;
  //if (j<i) { x=i; i=j; j=x; } // make i the lower index

  // no forces in hard spheres.
  return (0.0);
}
