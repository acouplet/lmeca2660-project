#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cfd.h"

#define P(x) {if (debug) printf("%lf\n",x);}
#define H(x) {if (debug) printf("%s = %lf\n",#x,x);}

const int debug = 1;
const double g = 9.81;

int main(){	
	clock_t begin = clock();

    int i,j;
    double L = 1.0;
    double H = 1.5;
    int Nx = L*50;
    int Ny = H*Nx;
    double h = 1.0/Nx;
    double dt = 0.01;
    double dtau = dt/1000;
    double l0 = 0.001*H;


    
    double **P     = matrix(Ny,Nx);
    double **T     = matrix(Ny,Nx);
    double **u     = matrix(Ny,Nx+2);
    double **v     = matrix(Ny+2,Nx);
    double **ustar = matrix(Ny,Nx+2);
    double **vstar = matrix(Ny+2,Nx);

    for(i=0;i<Ny;i++){
        for(j=0;j<Nx;j++){
            P[i][j] = g*((Ny-i-0.5)*h);
            T[i][j] = 0.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            v[Ny][j] = 0.0; v[Ny+1][j] = 0.0;
        }
        u[i][Nx] = 0.0; u[i][Nx+1] = 0.0;
        H(P[i][0]);
    }

    // LOOP for ustar
    for(i=0;i<Ny;i++){
        for(j=0;j<Nx+2;j++){
        }
    }


    free_matrix(P,Ny);
    free_matrix(T,Ny);
    free_matrix(u,Ny);
    free_matrix(v,Ny+2);
    
	
	clock_t end = clock();
	printf("%lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
}
