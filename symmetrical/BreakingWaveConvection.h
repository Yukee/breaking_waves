#ifndef BREAKINGWAVECONV_H
#define BREAKINGWAVECONV_H

#include "Flux.h"

class BreakingWaveConvection: public Flux
{
 private:

  SField m_velocity;

 public:

 BreakingWaveConvection(): Flux(2,1) {}

  inline void set_parameter(const SField & vel)
  {
    m_velocity = vel;
  }

  inline virtual VectorField evaluate(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_evaluated_flux_d[0] = u[0]*m_velocity;
      break;

    case 1:
      m_evaluated_flux_d[0] = -1*u[0]*(1-u[0]);
      break;

    default:
      throw std::invalid_argument("In BreakingWaveConvection::evaluate dimension must either be 0 or 1");
    }
    return m_evaluated_flux_d;
  }

  inline virtual SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_max_eigenvalue = m_velocity;
      break;

    case 1:
      m_max_eigenvalue = -1*(1-2*u[0]);
      break;

    default:
      throw std::invalid_argument("In BreakingWaveConvection::get_max_eigenvalue dimension must either be 0 or 1");
    }
    return m_max_eigenvalue;
  }

};

#endif
