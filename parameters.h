/* Simple solver for the conformal BDNK equations in slab symmetry.
 * 
 * Written by Alex Pandya, last updated 4/4/2024; numerical method
 * described in:
 * [1] https://arxiv.org/abs/2201.12317
 *
 * Choose options for the simulation in this file.
* */

/*============================define flags====================================*/
//NOTE: do not change these.  They just define an arbitrary number associated
//with the types of initial data and boundary conditions, so that we can
//write easily readable statements like ID_TYPE (GAUSSIAN) below rather than
//opaque ones like ID_TYPE (11)

//types of initial data
#define GAUSSIAN     (11)                                                       
#define STEP         (12)                                                       
#define SMOOTH_SHOCK (13)

//types of boundary conditions
#define GHOST    (20)                                                              
#define PERIODIC (21)

//include math.h because we use pow() below; define pi if it's undefined
#include <math.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/*============================set parameters==================================*/

//specify hydrodynamic frame.  This is done according to [1] eqs (19-20).
#define EPS0 (10.)
const double eta0    = 1.*pow(EPS0,0.25)/(3.*M_PI);
const double lambda0 = 25./7.*eta0;
const double chi0    = 25./4.*eta0;

//set simulation spatial resolution as RES_MULTIPLE*BASE_RESOLUTION.
//for a convergence test, one should only change RES_MULTIPLE; cell
//centers should be appropriately aligned
#define RES_MULTIPLE    (1)
#define BASE_RESOLUTION (128)

//this is equal to the number of timesteps when RES_MULTIPLE = 1;
//see the definition of MAX_TIMESTEP below
#define BASE_NUM_TIMESTEP (1024)

//set CFL number (simulation time resolution)
#define CFL (0.1)

//save data to file every TS_STEP timesteps
#define TS_STEP (10)

//set domain boundaries
#define X_MIN (-200)
#define X_MAX (200)

//choose the type of initial data from options GAUSSIAN, STEP (step
//function), and SMOOTH_SHOCK (approximate BDNK steady-state shock
//solution)
#define ID_TYPE (GAUSSIAN)

//set parameters for GAUSSIAN ID; unused if ID != GAUSSIAN
#define GAUSSIAN_AMPLITUDE (1)
#define GAUSSIAN_MEAN      (0)
#define GAUSSIAN_SPREAD    (25)
#define GAUSSIAN_CONST     (1e-1)

//set parameters for STEP ID; unused if ID != STEP
#define STEP_EPS_L (1)
#define STEP_EPS_R (0.1)

//set parameters for SMOOTH_SHOCK ID; unused if ID != SMOOTH_SHOCK
#define SMOOTH_SHOCK_EPS_L (1)
#define SMOOTH_SHOCK_V_L   (0.8)

//choose the type of boundary condition from options GHOST (ghost
//cells; roughly, outflow conditions) and PERIODIC (periodic
//boundaries)
#define BC (GHOST)

//choose the WENO/CWENO parameter epsilon > 0; determines how much 
//the derivative/reconstruction stencils avoid sharp features
#define epsW (1e-3)

//directory to save the data files in.  Data is saved with a new
//file for each variable being output, and data in the file is
//arranged such that the nth row contains the variable at timestep n.
//NOTE: this folder must already exist; if not, the code will segfault
#define DIREC ("datafiles/")

//tolerance below which we use perfect fluid primitive solve in a
//given cell rather than the BDNK primitive solve.  Ideally, the
//perfect fluid primitive solve would only ever be used if the
//BDNK solve is unstable at low resolution (in this case, perfect
//fluid solve should give answers that are close).  If we have
//sufficient resolution (or if the BDNK primitive solve is stable)
//we should always use BDNK primitive solve.  
//Set TOL < 0 to always use the BDNK primitive solve.
#define TOL (-1)

/*======================compute derived parameters============================*/

//solve system for N cells; without the +1, cell centers are 
//misaligned when doubling RES_MULTIPLE
#define N (RES_MULTIPLE*BASE_RESOLUTION+1)

//set MAX_TIMESTEP.  This is set up so that doubling RES_MULTIPLE
//doubles the resolution both in space and in time; once again,
//the +1 is included so that time levels are aligned when we
//double RES_MULTIPLE
#define MAX_TIMESTEP      (RES_MULTIPLE*BASE_NUM_TIMESTEP+1)
