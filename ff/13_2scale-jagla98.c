//#warning Note: Using Force Field Jagla98 
//Note to richard: (see paper J3, my notebook p35)
//#define lambda 0.606060606060606q // 1/1.65

# define GAMMA 1.76

// This was derived from 10_2scale-s2s
inline double eij(int i, int j, double d) {
  int x;
  if (j<i) { x=i; i=j; j=x; } // make i the lower index

  if (i==0 && j==0) // Jagla-Jagla
    {
      if(d < 1.)
	return (1./0.);
      else if (d < GAMMA)
	return ((GAMMA - d)/(GAMMA - 1.));
      else
	return(0.0);
    }
  else if (i==0 && j==1)
    {
      if(d < 1.)
	return (1./0.);
    }
  else if (i==1 && j==1)
    {
      if(d < 1.)
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
    if(d < 1.) {
      return (0.0);
    }
    else if (d < GAMMA) {
      return 1. / (GAMMA - 1.) ;
    }
    else {
      return(0.0);
    }
  }
  else { /* Other atom types can never exert force */
    return (0.0);
  }
}
