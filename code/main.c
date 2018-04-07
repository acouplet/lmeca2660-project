#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cfd.h"

#define g 9.81

int main(){	
	clock_t begin = clock();

    int i,j;
    int Nx = 50;
    int Ny = 1.5*50;
    double h = 1.0/Nx;


    
    double **P = matrix(Ny,Nx);
    double **T = matrix(Ny,Nx);
    double **u = matrix(Ny,Nx+2);
    double **v = matrix(Ny+2,Nx);

    for(i=0;i<Ny;i++){
        for(j=0;j<Nx;j++){
            P[i][j] = g*((Ny-i-0.5)*h);
        }
        printf("%lf\n",P[i][0]);
    }



    free_matrix(P,Ny);
    free_matrix(T,Ny);
    free_matrix(u,Ny);
    free_matrix(v,Ny+2);
    
	
	clock_t end = clock();
	printf("%lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
}
