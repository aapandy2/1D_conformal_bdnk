# Simple solver for the conformal BDNK equations

**Alex Pandya**

_April 3rd, 2024_

This is a simplified `C` code that solves the conformal BDNK equations in flat spacetime in Cartesian coordinates under the assumption of slab symmetry (solutions invariant under translations in $y, z$).  For details on the algorithm, see:

[1] https://arxiv.org/abs/2201.12317

## Instructions

1. compile the code with `make`
2. make the directory to save datafiles in (by default, the directory is called `datafiles` and should be placed in the same directory as the code): `mkdir datafiles/`
3. choose simulation parameters in `parameters.c` (default set of parameters gives the evolution of a Gaussian clump of energy density)
4. run the code with `./1d_solver`
5. data is stored in plaintext format with a separate file for each variable (e.g., `eps.txt` contains the solution for $\epsilon(t, x)$, the energy density).  Within a datafile, rows are at a constant time and columns are at constant spatial coordinate $x$.
