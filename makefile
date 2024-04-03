1D_solver: solver_1D.c
	gcc -std=c11 -o 1D_solver solver_1D.c -lm -I . -Ofast
