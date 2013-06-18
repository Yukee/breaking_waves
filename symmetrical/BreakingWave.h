#ifndef BREAKINGWAVE_H
#define BREAKINGWAVE_H
#include <math.h>
#include <iostream>

const static double Hup = 0.9; const static double Hdown = 0.1; const static double xa=0.8;

// velocity field in the travelling frame

double u(double z)
{
  return 2*z -1;
}

// initial concentration of small particules

double phi0(double x, double z)
{
  double phi0 = 0;
  double bound = Hup - (1/0.8)*(x + xa); 
  
  if(x <= -xa && z <= Hup) phi0 = 1;
  else if(x > -xa && x <= xa && z <= bound) phi0 = 1;
  else if(x > xa && z <= Hdown) phi0 = 1;

  return phi0;
}

#endif
