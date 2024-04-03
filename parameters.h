/* Simple solver for the conformal BDNK equations in slab symmetry.
 * 
 * Written by Alex Pandya, last updated 4/3/2024; numerical method
 * described in https://arxiv.org/abs/2201.12317
 *
 * Choose options for the simulation in this file.
* */

/*============================define flags====================================*/

//types of initial data
#define GAUSSIAN     (11)                                                       
#define SHOCK        (12)                                                       
#define SMOOTH_SHOCK (17)

//types of boundary conditions
#define GHOST (20)                                                              
#define PERIODIC (21)

/*============================set parameters==================================*/

//set simulation spatial resolution as RES_MULTIPLE*BASE_RESOLUTION.
//for a convergence test, one should only change RES_MULTIPLE; cell
//centers should be appropriately aligned
#define RES_MULTIPLE    (1)
#define BASE_RESOLUTION (128)

//set CFL number (simulation time resolution)
#define CFL (0.1)

//set maximum number of timesteps
#define max_timestep (1*1024+1)

//save data to file every TS_STEP timesteps
#define TS_STEP (10)

//tolerance below which we use perfect fluid primitive solve in a
//given cell rather than the BDNK primitive solve.  Ideally, the
//perfect fluid primitive solve would only ever be used if the
//BDNK solve is unstable at low resolution (in this case, perfect
//fluid solve should give answers that are close).  If we have
//sufficient resolution (or if the BDNK primitive solve is stable)
//we should always use BDNK primitive solve.  
//Set TOL < 0 to always use the BDNK primitive solve.
#define TOL (-1)

//choose the type of initial data from options GAUSSIAN, SHOCK (step
//function), and SMOOTH_SHOCK (approximate BDNK steady-state shock
//solution)
#define ID_TYPE (GAUSSIAN)

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
#define DIREC "low/"

/*======================compute derived parameters============================*/

//solve system for N cells; without the +1, cell centers are 
//misaligned when doubling RES_MULTIPLE
#define N (RES_MULTIPLE*BASE_RESOLUTION+1)
