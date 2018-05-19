#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cfd.h"
#include "Problem.h"

#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
#define DR(i,a,b) for (int i=(b); i-->(a);)
#define D(i,n) DR(i,0,n)

int main(int argc, char *argv[]){
    //Problem *problem = Default_Problem();
    
    // New_Problem(int Nx, int dt, int saveIter, int useMixer)
    Problem *problem;
    if (argc == 5)
		problem = New_Problem(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
	else 
		problem = New_Problem(100,100,20,0);
	problem->Description(problem);

	int Nx = problem->Nx;
	int Ny = problem->Ny;
	double **T = problem->T;
	double **u = problem->u;
	double **v = problem->v;
	double **HnpT = problem->HnpT;
	double h = problem->h;
	        
    // Initial
    problem->BoundaryConditions(problem);
    problem->Mixer(problem);
    problem->Diagnostics(problem);
    problem->WriteData(problem);
    problem->IterationInfo(problem);

    while(problem->averageTemp < 3e-3 && problem->iter < 100000 && problem->Reh <= 40 && problem->Rehw <= 25){
        problem->BoundaryConditions(problem);
        problem->Mixer(problem);
        problem->Momentum(problem);
        problem->Energy(problem);
        problem->Poisson(problem);
        problem->SpeedUpdate(problem);
        problem->PressureUpdate(problem);
        problem->iter++;
        problem->WriteData(problem);
        problem->IterationInfo(problem);
    }
    
    if (problem->Reh > 40)
		printf("Simulation stopped: Reh = %.2lf > 40\n", problem->Reh);
		
	if (problem->Rehw > 25)
		printf("Simulation stopped: Rehw = %.2lf > 25\n", problem->Rehw);
		
    
    return(0);
}
