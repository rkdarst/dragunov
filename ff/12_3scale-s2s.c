/* Richard Darst, 2007 */

//#warning Note: Using 3-scale model with S-2S units
#define lambda 0.581395348837209  // 1/1.72

inline double eij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

  // gnuplot> plot [.58:1.74] ((4.5/(10.75*(l-1)))*x - (3.5+l)/(10.75*(l-1))), (1/(10.75*(3*l-1))*(x-3*l)  ), 0, -1/10.75  

  if (i==0 && j==0) // Jagla-Jagla
    {
      if(d < lambda)
	return (1./0.);
      else if (d < 1.)
	return ((1.-d) - (1./(10.75)));
      else if (d < 3.*lambda)
	return ((1./(10.75*(3*lambda-1)))*(d - 3*lambda));
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
  return(0.0);
}



inline double fij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

  if (i==0 && j==0) { // Jagla-Jagla
    if(d < lambda)
      return (0.0);
    else if (d < 1.)
      return 1. ;
    else if (d < 3.*lambda)
      return( -1/(10.75*(3*lambda-1)) );
    else
      return(0.0);
  } // end jagla particle interaction
  else { // Other atom types can never exert force
    return (0.0);
  }
}
