#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>

// Solvers
#include "FD1Solver.h"
#include "timeSolver.h"
#include "EulerSolver.h"
//#include "RK2Solver.h"
#include "RK3Solver.h"

// Containers
#include "Vector.h"
#include "ScalarField.h"
#include "PrescribedField.h"
#include "Equation.h"

// Functions
#include "Flux.h"
#include "ZeroFlux.h"
#include "BreakingWaveConvection.h"

#include "WriteVectorField.h"

// Breaking wave pb
#include "BreakingWave.h"

using namespace std;
int main(int argc, char *argv[])
{
    // Prescribe the fluxes and source term, store them in an equation
    
   Flux *ptrCF = new BreakingWaveConvection();
   Flux *ptrDF = new ZeroFlux(2,1);
   Flux *ptrS = new ZeroFlux(2,1);
   Equation break_eq (ptrCF, ptrDF, ptrS);
   
   // Ask the user for the discretisation and timestep infos
   
   Vector<double> dx (2); Vector<double> xi (2); Vector<double> llc (2); double endtime; double timestep; double timebtwfiles; string filename;
   //cout << "Enter cell width dx:"; cin >> dx[0];
   dx[0] = 0.05;
   
   //cout << "Enter cell heigh dz:"; cin >> dx[1];
   dx[1] = 0.01;
   
   //cout << "Enter domain width:"; cin >> xi[0];
   xi[0] = 3;
   xi[1] = 1; llc[0] = -xi[0]/2; llc[1] = 0;
   
   cout << "Enter end time:"; cin >> endtime;
   
   //cout << "Enter time step:"; cin >> timestep;
   timestep = 0.05;
   
   cout << "Enter time between two consecutive file saves:"; cin >> timebtwfiles;
   
   cout << "Enter name for saved file:"; cin >> filename;
   
   FD1Solver *solv_ptr = new FD1Solver (dx, xi, &break_eq, llc);
   VectorField pos = solv_ptr->get_position();
   Vector<int> xr = solv_ptr->get_nxSteps();

   // Initialise value of the concentration in small particules, and velocity field

   VectorField phi (1, SField (xr));
   SField u0 (xr);
   
   for(int it=0;it<phi[0].get_size();++it)
   {
    phi[0][it] = phi0(pos[0][it], pos[1][it]);
    u0[it] = u(pos[1][it]);
   }
   
   // Set the west and east boundary conditions
   
   ScalarField phiWest (xr.drop(0)); 
   for(int it=0;it<phiWest.get_size();++it) phiWest[it] = phi0( llc[0]-dx[0] , llc[1]+it*dx[1] );
   phi[0].set_bound(0,-1,phiWest);

   ScalarField phiEast (xr.drop(0)); 
   for(int it=0;it<phiEast.get_size();++it) phiEast[it] = phi0( llc[0]+(xr[0]+1)*dx[0] , llc[1]+it*dx[1] );
   phi[0].set_bound(0,1,phiEast);

   // Sets the parameters of the convection flux: velocity field and segregation rate
   
   ptrCF->set_parameter(u0); 
   
   // Check if the initial setup is correct
   fstream phiinit;
   phiinit.open("Results/Check/phiinit.tsv",ios::out);
   phi[0].write_in_file(phiinit, dx, llc);
   
   fstream phiwest;
   phiwest.open("Results/Check/phiwest.tsv",ios::out);
   phiWest.write_in_file(phiwest, dx.drop(0), llc.drop(0));

   fstream phieast;
   phieast.open("Results/Check/phieast.tsv",ios::out);
   phiEast.write_in_file(phieast, dx.drop(0), llc.drop(0));
   
   fstream uinit;
   uinit.open("Results/Check/uinit.tsv",ios::out);
   u0.write_in_file(uinit, dx, llc);

   // Timestepper
   
   RK3Solver tsolv (timestep, endtime, solv_ptr, phi);
   
   // Compute the time evolution
   
   //tsolv.get_solution(filename, timebtwfiles);

   return 0;
}
