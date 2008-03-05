/* Richard Darst, 2007 */

//#warning Note: Using 3-scale model with ***S-3S*** units
#define lambda 0.581395348837209  // 1/1.72

inline double eij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

// gnuplot> plot [0:1.74] (-1/(4.5*(1-3*l)) * x + 3*l/(4.5*(1-3*l))), +
//              (1-x)-1/4.5, -1/4.5, 0 
  if (i==0 && j==0) // Jagla-Jagla
    {
      // U_0 = 1 ; U_a == -1
      // a = 1
      if(d < 1.)
	return (1./0.);
      else if (d < 1.72)
	return (-1 + (-1. - 3.5)*(d-1.72)/(1.72-1.));
      else if (d < 3.)
	return (-1*(3.0-d)/(3.0-1.72));
      else
	return(0.0);
    }
  else if (i==0 && j==1)
    {
      if(d < 1.0)
	return (1./0.);
    }
  else if (i==1 && j==1)
    {
      if(d < 1.0)
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
    if(d < 1.) {
      return (0.0);
    }
    else if (d < 1.72) {
      return - (-(-1. - 3.5)/(1.72-1.));
    }
    else if (d < 3.) {
      return - (-(1/(3.0-1.72))); // this is actually a triple negative
    }
    else {
      return(0.0);
    }
  }
  else { /* Other atom types can never exert force */
    return (0.0);
  }
}

