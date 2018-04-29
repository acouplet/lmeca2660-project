#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cfd.h"

#define P(x) {if (debug) printf("%lf\n",x);}
#define H(x) {if (debug) printf("%s = %lf\n",#x,x);}
#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
#define DR(i,a,b) for (int i=(b); i-->(a);)
#define D(i,n) DR(i,0,n)

const int debug = 1;
const double g = 9.81;

int main(){	
	clock_t begin = clock();

    int Nx       = 40;
    int Ny       = 60;
    double h     = 1.0/Ny;
    double dt    = 1e-2;
    double dtau  = dt/1e3;
    double l0    = 1e-3;
    double Gr    = 2.0e10;
    double alpha = 1.975;

    double sumR = 0.0;
    double globe = 0.0;
    double **P       = matrix(Ny,Nx);
    double **T       = matrix(Ny,Nx);
    double **u       = matrix(Ny+2,Nx+1);
    double **v       = matrix(Ny+1,Nx+2);
    double **ustar   = matrix(Ny+2,Nx+1);
    double **vstar   = matrix(Ny+1,Nx+2);
    double **Hnpx    = matrix(Ny,Nx);
    double **Hnpy    = matrix(Ny,Nx);
    double **Phi     = matrix(Ny+2,Nx+2);
    double **Phistar = matrix(Ny+2,Nx+2);
    double **R       = matrix(Ny,Nx);

    F(i,Ny) F(j,Nx) P[i][j] = g*(Ny-i-0.5)*h; // Probably wrong, any idea how to make initial P dimensionless

    //================//
    // First equation //
    //================//
    F(i,Ny+2) { u[i][0] = 0; u[i][Nx] = 0; } // No-slip condition 
    F(j,Nx+2) { v[0][j] = 0; v[Ny][j] = 0; }
    F(j,Nx+1) { u[0][j] = u[1][j]; u[Ny+1][j] = u[Ny][j]; } // Ghost points, no-through flow
    F(i,Ny+1) { v[i][0]    = -0.2*(v[i][3]    - 5*v[i][2]    + 15*v[i][1]);
                v[i][Nx+1] = -0.2*(v[i][Nx-2] - 5*v[i][Nx-1] + 15*v[i][Nx]); }

    FR(i,1,Ny+1){
        FR(j,1,Nx){
            ustar[i][j]    = -1*NS_pressurex(i,j,P,h);
            ustar[i][j]   += (1/sqrt(Gr))*NS_diffusionx(i,j,u,h);
            double H       = NS_convectionx(i,j,u,v,h);
            ustar[i][j]   -= 0.5*(3*H - Hnpx[i-1][j-1]);
            Hnpx[i-1][j-1] = H;
            ustar[i][j]    = dt*ustar[i][j] + u[i][j];
        }
    }
    F(j,Nx+1) { ustar[0][j] = ustar[1][j]; ustar[Ny+1][j] = ustar[Ny][j]; } // Ghost points

    FR(i,1,Ny){
        FR(j,1,Nx+1){
            vstar[i][j]    = -1*NS_pressurey(i,j,P,h);
            vstar[i][j]   += (1/sqrt(Gr))*NS_diffusiony(i,j,v,h);
            vstar[i][j]   -= NS_buoyancy(i,j,T);
            double H       = NS_convectiony(i,j,u,v,h);
            vstar[i][j]   -= 0.5*(3*H - Hnpy[i-1][j-1]);
            Hnpy[i-1][j-1] = H;
            vstar[i][j]    = dt*vstar[i][j] + v[i][j];
        }
    }
    F(i,Ny+1) { vstar[i][0]    = -0.2*(vstar[i][3]    - 5*vstar[i][2]    + 15*vstar[i][1]);
                vstar[i][Nx+1] = -0.2*(vstar[i][Nx-2] - 5*vstar[i][Nx-1] + 15*vstar[i][Nx]); }

    globe = 1.0;
	clock_t SORtb = clock();
    while(globe > 1e-3){
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                Phistar[i][j] = (h*h)/4*(-1/dt*((ustar[i][j]-ustar[i][j-1])/h + (vstar[i-1][j]-vstar[i][j])/h) + (Phi[i+1][j] + Phi[i-1][j])/(h*h) + (Phi[i][j+1] + Phi[i][j-1])/(h*h));
                Phi[i][j] = alpha*Phistar[i][j] + (1-alpha)*Phi[i][j];
            }
        }
        sumR = 0;
        F(i,Ny){
            F(j,Nx){
                R[i][j] = (Phi[i+2][j+1]-2*Phi[i+1][j+1]+Phi[i][j+1])/(h*h) + (Phi[i+1][j+2]-2*Phi[i+1][j+1]+Phi[i+1][j])/(h*h) - 1/dt*((ustar[i+1][j+1]-ustar[i+1][j])/h + (vstar[i][j+1]-vstar[i+1][j+1])/h);
                sumR += R[i][j]*R[i][j];
            }
        }
        globe = dt*sqrt((3*h*h)/2*sumR);

        printf("Sum local residual = %lf\n",sumR);
        printf("Global error = %lf\n",globe);
    }
	clock_t SORte = clock();


    F(i,Ny+2){ 
        F(j,Nx+1)
            printf("%1.2lf\t ",ustar[i][j]);
        printf("\n");
    }
    F(i,Ny+1){ 
        F(j,Nx+2)
            printf("%1.2lf\t ",vstar[i][j]);
        printf("\n");
    }

    F(i,Ny+2){
        F(j,Nx+2)
            printf("%1.2lf\t ",Phistar[i][j]);
        printf("\n");
    }

        


    free_matrix(P,Ny);
    free_matrix(T,Ny);
    free_matrix(u,Ny+2);
    free_matrix(v,Ny+1);
    free_matrix(ustar,Ny+2);
    free_matrix(vstar,Ny+1);
    free_matrix(Hnpx,Ny);
    free_matrix(Hnpy,Ny);
    
	
	clock_t end = clock();
	printf("Total time %lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
	printf("SOR time %lf ms (%.2lf of total time)\n",(double)(SORte - SORtb) / CLOCKS_PER_SEC * 1000, (double)(SORte-SORtb)/(double)(end-begin));
}
