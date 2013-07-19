#ifndef BREAKINGWAVE_H
#define BREAKINGWAVE_H
#include <math.h>
#include <iostream>

const static double Hup = 0.9; const static double Hdown = 0.1; const static double xa=0.8;

// velocity field in the travelling frame

// Looks symmetrical to me:
// y = 0, ... 0.495, 0.5, 0.505, ... 1
// v = -1,...-0.001,   0, 0.001, ... 1    
double u(double z)
{
  return 2*z -1;
}

// initial concentration of small particules

double phi0(double x, double z)
{
  double phi0 = 0;
  double bound = 0.5*(Hup + Hdown - ((Hup - Hdown)/xa)*x); 
  
  if(x <= -xa && z < Hup) phi0 = 1;
  if(x > -xa && x <= xa && z <= bound) phi0 = 1;
  if(x > xa && z <= Hdown) phi0 = 1;

  return phi0;
}

#endif
