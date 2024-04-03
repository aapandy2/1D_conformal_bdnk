1D_solver: solver_1D.c parameters.h
	gcc -std=c11 -o 1D_solver solver_1D.c parameters.h -lm -I . -O2
