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

const int debug = 0;
const double g = 9.81;

int main(){	
	clock_t begin = clock();

    int Nx       = 400;
    int Ny       = 1.5*Nx;
    double h     = 1.0/(Ny-1); // 
    double Pr    = 2.0;
    double dt    = 0.1;
    double dtau  = dt/1e3;
    double l0    = 1e-3;
    double Tinf  = -5e-3;
    double Gr    = 2.0e10;
    double alpha = 1.97;
	double H;

    double sumR      = 0.0;
    double globe     = 0.0;
    double Tavg      = 0.0;
    double Trms      = 0.0;
	clock_t SORtb    = clock();
	clock_t SORte    = clock();
    double SORtime   = 0.0;
    double **P       = matrix(Ny,Nx);
    double **T       = matrix(Ny+2,Nx+2);
    double **u       = matrix(Ny+2,Nx+1);
    double **v       = matrix(Ny+1,Nx+2);
    double **ustar   = matrix(Ny+2,Nx+1);
    double **vstar   = matrix(Ny+1,Nx+2);
    double **Hnpx    = matrix(Ny,Nx); // Change to Nx-1
    double **Hnpy    = matrix(Ny,Nx); // Change to Ny-1
    double **Phi     = matrix(Ny+2,Nx+2); // Change to (Ny,Nx)
    double **Phistar = matrix(Ny+2,Nx+2); // Change to (Ny,Nx)
    double **R       = matrix(Ny,Nx);
    double **HnpT    = matrix(Ny,Nx);

	FILE *fP = fopen("data/P_400_600.csv","w");
	FILE *fT = fopen("data/T_400_600.csv","w");
	FILE *fu = fopen("data/u_400_600.csv","w");
	FILE *fv = fopen("data/v_400_600.csv","w");
    
    F(k,1000){
    //================//
    // First equation //
    //================//
	// No-slip conditions not necessary, never modified
    //F(i,Ny+2) { u[i][0] = 0; u[i][Nx] = 0; } // No-slip condition 
    //F(j,Nx+2) { v[0][j] = 0; v[Ny][j] = 0; }
    //F(j,Nx+1) { u[0][j] = u[1][j]; // Ghost points, no-through flow
	//			u[Ny+1][j] = -0.2*(u[Ny-2][j] - 5*u[Ny-1][j] + 15*u[Ny][j]); } // CORRECTED Ghost point (prev: du/dy = 0, now: u_gamma = 0)
    //F(i,Ny+1) { v[i][0]    = -0.2*(v[i][3]    - 5*v[i][2]    + 15*v[i][1]);
    //            v[i][Nx+1] = -0.2*(v[i][Nx-2] - 5*v[i][Nx-1] + 15*v[i][Nx]); }

    FR(i,1,Ny+1){
        FR(j,1,Nx){
            ustar[i][j]    = -1*NS_pressurex(i,j,P,h);
            ustar[i][j]   += (1/sqrt(Gr))*NS_diffusionx(i,j,u,h);
            H       	   = NS_convectionx(i,j,u,v,h);
            ustar[i][j]   -= 0.5*(3*H - Hnpx[i-1][j-1]);
            Hnpx[i-1][j-1] = H;
            ustar[i][j]    = dt*ustar[i][j] + u[i][j];
        }
    }
    F(j,Nx+1) { ustar[0][j] = ustar[1][j]; ustar[Ny+1][j] = ustar[Ny][j]; } // WRONG, SEE CORRECTED ABOVE

    FR(i,1,Ny){
        FR(j,1,Nx+1){
            vstar[i][j]    = -1*NS_pressurey(i,j,P,h);
            vstar[i][j]   += (1/sqrt(Gr))*NS_diffusiony(i,j,v,h);
            vstar[i][j]   -= NS_buoyancy(i,j,T);
            H       	   = NS_convectiony(i,j,u,v,h);
            vstar[i][j]   -= 0.5*(3*H - Hnpy[i-1][j-1]);
            Hnpy[i-1][j-1] = H;
            vstar[i][j]    = dt*vstar[i][j] + v[i][j];
        }
    }
    F(i,Ny+1) { vstar[i][0]    = -0.2*(vstar[i][3]    - 5*vstar[i][2]    + 15*vstar[i][1]);
                vstar[i][Nx+1] = -0.2*(vstar[i][Nx-2] - 5*vstar[i][Nx-1] + 15*vstar[i][Nx]); }

    //============//
    // Equation 5 //
    //============//
    FR(j,1,Nx+1){
        T[Ny+1][j] = T[Ny][j] + h; // Ghost point, bottom heat transfer
        T[0][j] = (-h/l0)*(1.5*T[1][j] - 0.5*T[2][j] - Tinf) + T[1][j];
    }

    Tavg = 0.0;
    Trms = 0.0;
    FR(i,1,Ny+1){
        FR(j,1,Nx+1){
            H       	   = (T[i][j+1]-T[i][j-1])*(u[i][j]+u[i][j-1])/(2*h) + (T[i-1][j] - T[i+1][j])*(v[i-1][j]+v[i][j])/(2*h);
            double Tn1     = -0.5*(3*H - HnpT[i-1][j-1]);
            HnpT[i-1][j-1] = H;
            Tn1           += (1/(Pr*sqrt(Gr)*h*h))*((T[i][j+1]-2*T[i][j]+T[i][j-1])+(T[i-1][j]-2*T[i][j]+T[i+1][j]));
            T[i][j]        = dt*Tn1 + T[i][j];
            Tavg          += T[i][j];
        }
    }
    Tavg /= Nx*Ny;
    FR(i,1,Ny+1) FR(j,1,Nx+1) Trms += (T[i][j] - Tavg)*(T[i][j] - Tavg);
    Trms /= Nx*Ny;
    FR(i,1,Ny+1){
        T[i][0] = T[i][1]; // Ghost points, adiabatic
        T[i][Nx+1] = T[i][Nx];
    }



    //============//
    // SOR Method //
    //============//
    globe = 1.0;
    int iter = 0;
	SORtb = clock();
    while(globe > 1e-3){
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
				// Ghost points for PHI ? 
                Phistar[i][j] = (h*h)/4*(-1/dt*((ustar[i][j]-ustar[i][j-1])/h + (vstar[i-1][j]-vstar[i][j])/h) + (Phi[i+1][j] + Phi[i-1][j])/(h*h) + (Phi[i][j+1] + Phi[i][j-1])/(h*h));
                Phi[i][j] = alpha*Phistar[i][j] + (1-alpha)*Phi[i][j];
            }
        }
        sumR = 0;
        F(i,Ny){
            F(j,Nx){
				// VERIFY INDICES
                R[i][j] = (Phi[i+2][j+1]-2*Phi[i+1][j+1]+Phi[i][j+1])/(h*h) + (Phi[i+1][j+2]-2*Phi[i+1][j+1]+Phi[i+1][j])/(h*h) - 1/dt*((ustar[i+1][j+1]-ustar[i+1][j])/h + (vstar[i][j+1]-vstar[i+1][j+1])/h);
                sumR += R[i][j]*R[i][j];
            }
        }
        globe = dt*sqrt((3*h*h)/2*sumR);
        iter += 1;
        //printf("Sum local residual = %lf\n",sumR);
        //printf("Global error = %lf\n",globe);
    }
	SORte = clock();
    SORtime += (double)(SORte - SORtb) / CLOCKS_PER_SEC * 1000;

    //============//
    // Equation 3 //
    //============//
    FR(i,1,Ny+1)
        FR(j,1,Nx)
            u[i][j] = -(dt/h)*(Phi[i][j+1]-Phi[i][j]) + ustar[i][j];

    FR(i,1,Ny)
        FR(j,1,Nx+1)
            v[i][j] = -(dt/h)*(Phi[i][j]-Phi[i+1][j]) + vstar[i][j];

	// No-slip conditions not necessary, never modified
    //F(i,Ny+2) { u[i][0] = 0; u[i][Nx] = 0; } // No-slip condition 
    //F(j,Nx+2) { v[0][j] = 0; v[Ny][j] = 0; }
    F(j,Nx+1) { u[0][j] = u[1][j]; // Ghost points, no-through flow
				u[Ny+1][j] = -0.2*(u[Ny-2][j] - 5*u[Ny-1][j] + 15*u[Ny][j]); } // CORRECTED Ghost point (prev: du/dy = 0, now: u_gamma = 0)
    F(i,Ny+1) { v[i][0]    = -0.2*(v[i][3]    - 5*v[i][2]    + 15*v[i][1]);
                v[i][Nx+1] = -0.2*(v[i][Nx-2] - 5*v[i][Nx-1] + 15*v[i][Nx]); }

    //============//
    // Equation 4 //
    //============//
    F(i,Ny)
        F(j,Nx)
            P[i][j] += Phi[i+1][j+1];
    

    printf("SOR Iterations = %d\n",iter);
    printf("Tavg = %.4lf\n",Tavg);
    printf("Trms = %.4lf\n",Trms);

	/*
    FR(i,0,Ny+2){
        FR(j,0,Nx+2)
            printf("%.4lf\t",T[i][j]);
        printf("\n");
    }
    */
    
    
    MatrixToFile(fP, P, Ny, Nx);
    MatrixToFile(fT, T, Ny+2, Nx+2);
    MatrixToFile(fu, u, Ny+2, Nx+1);
    MatrixToFile(fv, v, Ny+1, Nx+2);

    }


	fclose(fP);
	fclose(fT);
	fclose(fu);
	fclose(fv);
    free_matrix(P,Ny);
    free_matrix(T,Ny+2);
    free_matrix(u,Ny+2);
    free_matrix(v,Ny+1);
    free_matrix(ustar,Ny+2);
    free_matrix(vstar,Ny+1);
    free_matrix(Hnpx,Ny);
    free_matrix(Hnpy,Ny);
    free_matrix(Phi,Ny+2);
    free_matrix(Phistar,Ny+2);
    free_matrix(R,Ny);
    free_matrix(HnpT,Ny);
    
	
	clock_t end = clock();
	printf("Total time %lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
	printf("SOR time %lf ms (%.2lf of total time)\n",SORtime, SORtime/(double)(end-begin));
}
