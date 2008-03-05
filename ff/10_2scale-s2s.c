//#warning Note: Using Force Field 2-scale Stanly-2S units.
#define lambda 0.571428571428571 // 4/7
inline double eij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

  if (i==0 && j==0) // Jagla-Jagla
    {
      if(d < lambda)
	return (1./0.);
      else if (d < 1.)
	return (1. - d);
      else
	return(0.0);
    }
  else if (i==0 && j==1)
    {
      if(d < lambda)
	return (1./0.);
    }
  else if (i==1 && j==1)
    {
      if(d < lambda)
	return (1./0.);
    }
  else
    {
      printf("unknown atom types (%d, %d), quitting.\n", i, j);
      exit(145);
    }
  return(0.0); // we shouldn't get here.
}

inline double fij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

  if (i==0 && j==0) {
    if(d < lambda) {
      return (0.0);
    }
    else if (d < 1.) {
      return 1. ;
    }
    else {
      return(0.0);
    }
  }
  else { /* Other atom types can never exert force */
    return (0.0);
  }
}
