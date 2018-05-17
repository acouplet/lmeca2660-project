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

    // Initial
    problem->BoundaryConditions(problem);
    problem->Mixer(problem);
    problem->Diagnostics(problem);
    problem->IterationInfo(problem);
    problem->WriteData(problem);

    while(problem->averageTemp < 3e-3){
        problem->BoundaryConditions(problem);
        problem->Mixer(problem);
        problem->Momentum(problem);
        problem->Energy(problem);
        problem->Poisson(problem);
        problem->SpeedUpdate(problem);
        problem->PressureUpdate(problem);
        problem->iter++;
        problem->IterationInfo(problem);
        problem->WriteData(problem);
    }
    return(0);
}
