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

double divergence(double **u, double **v, double h, int Nx, int Ny){
    // du/dx + dv/dy = 0
    double div, maxdiv = -42;
    double maxu=-42, maxv=-42;
    int i,j;
    FR(i,1,Ny+1){
        FR(j,1,Nx+1){
            div = ((u[i][j]-u[i][j-1])/h + (v[i-1][j]-v[i][j])/h);
            if(div > maxdiv)
                maxdiv = div;
            if(u[i][j] > maxu) maxu = u[i][j];
            if(v[i][j] > maxv) maxv = v[i][j];
        }
    }
    printf("maxu=%.6f, maxv=%.6f\n", maxu, maxv);
    return maxdiv;
}

int main(int argc, char *argv[]){
	clock_t begin = clock();

    int Nx       = 100;
    int Ny       = 1.5*Nx;
    double h     = 1.0/Ny;
    double Pr    = 2.0;
    double dt    = 0.01;
    double dtau  = dt/1e3;
    double l0    = 1e-3;
    double Tinf  = -5e-3;
    double Gr    = 2.0e10;
    double alpha = 1.97;
    int miter    = 2000;
    int usemixer = 1;
    int saveIter = 50;
	double H;
    // mixer radius
    double a = 3./10. * Nx * h;
    double omega = 0.1; // U/H = 1 right? TO CHECK 
    double xg = 1.0/3.0; 
    double yg = 1.0/3.0;
	
	if (argc == 6){ 
		Nx    = atoi(argv[1]);
		dt    = 1/(double)atoi(argv[2]);
		miter = atoi(argv[3]);
        saveIter = atoi(argv[4]);
        usemixer = atoi(argv[5]);
		Ny    = 1.5*Nx;
		h     = 1.0/Ny;
	}

    double sumR      = 0.0;
    double globe     = 0.0;
    double Tavg      = 0.0;
    double Trms      = 0.0;
	clock_t SORtb    = clock();
	clock_t SORte    = clock();
    double SORtime   = 0.0;
    double **P       = matrix(Ny,Nx);
    double **T       = matrix(Ny+2,Nx+2);
    double **T_tmp   = matrix(Ny+2,Nx+2);
    double **u       = matrix(Ny+2,Nx+1);
    double **v       = matrix(Ny+1,Nx+2);
    double **normv   = matrix(Ny+1,Nx+1);
    double **vorty   = matrix(Ny+1,Nx+1);
    // terme source du mixer
    double **umixer  = matrix(Ny+2,Nx+1);
    double **vmixer  = matrix(Ny+1,Nx+2);
    double **ustar   = matrix(Ny+2,Nx+1);
    double **vstar   = matrix(Ny+1,Nx+2);
    double **xiu     = matrix(Ny+2,Nx+1);
    double **xiv     = matrix(Ny+1,Nx+2);
    double **xiT     = matrix(Ny+2,Nx+2);
    double **Hnpx    = matrix(Ny,Nx);
    double **Hnpy    = matrix(Ny,Nx);
    double **Phi     = matrix(Ny+2,Nx+2); 
    double **Phistar = matrix(Ny+2,Nx+2);
    double **R       = matrix(Ny,Nx);
    double **HnpT    = matrix(Ny,Nx);
    double mixer_angle = 0;
    
    char filename[54];
    FILE *fTb,*fwb,*fvb,*fdiag;
    sprintf(filename,"data/diagnostics_Nx%d_dt%d_mixing%d.bin",Nx,(int)(1/dt),usemixer);
    fdiag = fopen(filename,"wb");


    // diagnostics
    double average_temp;  
    double rms_temp; 
    double avg_temp_mixer;
    double avg_heat_flux;
    


    int k = 0;
    while((miter != 0 && k < miter) || (miter == 0 && Tavg < 3e-3)){
        //=======//
        // Mixer //
        //=======//
        double x, y, theta, d;
        double average_mixer_temp = 0;
        int nmixer_temp = 0;
        
        // xiu
        FR(i,1,Ny+1){
            FR(j,1,Nx){
                x = j*h; y = (Ny-i+0.5)*h;
                theta = atan2(y-yg,x-xg);
                d = sqrt(pow(x-xg,2) + pow(y-yg,2));
                xiu[i][j] = ( d <= a * cos(3*theta + mixer_angle) ) && usemixer;
                umixer[i][j] = -sin(theta) * omega * d;
            }
        }

        // xiv
        FR(i,1,Ny){
            FR(j,1,Nx+1){
                x = (j-0.5)*h; y = (Ny-i)*h;
                theta = atan2(y-yg,x-xg);
                d = sqrt(pow(x-xg,2) + pow(y-yg,2));
                xiv[i][j] = (d <= a*cos(3*theta + mixer_angle)) && usemixer;
                vmixer[i][j] =  cos(theta) * omega * d;
            }
        }

        // xiT and average_mixer_temp
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                x = (j-0.5)*h; y = (Ny-i+0.5)*h;
                theta = atan2(y-yg,x-xg);
                d = sqrt(pow(x-xg,2) + pow(y-yg,2));
                xiT[i][j] = (d <= a*cos(3*theta + mixer_angle)) && usemixer;
                average_mixer_temp += T[i][j] * (d <= a);
                nmixer_temp += (d <= a);
            }
        }

        mixer_angle += dt * omega;
        avg_temp_mixer = average_mixer_temp / nmixer_temp;

        
        //================//
        // First equation //
        //================//
        FR(i,1,Ny+1){
            FR(j,1,Nx){
                ustar[i][j]    = -1*NS_pressurex(i,j,P,h);
                ustar[i][j]   += (1/sqrt(Gr))*NS_diffusionx(i,j,u,h);
                H       	   = NS_convectionx(i,j,u,v,h);
                ustar[i][j]   -= 0.5*(3*H - Hnpx[i-1][j-1]);
                Hnpx[i-1][j-1] = H;
                // vstar = 1 / (1 + xhi Dt/Dtau)   *  (vn + dt * blabla + xhi Dt/Dtau v_s)
                // xsi calculé sur l'interieur du domaine
                ustar[i][j]    = 1/(1 + xiu[i][j] * dt/dtau) * (u[i][j] + dt * ustar[i][j]  + xiu[i][j] * dt/dtau * umixer[i][j] );
            }
        }
        
        F(j,Nx+1) { ustar[0][j] = ustar[1][j]; // Ghost points, no-through flow
    				ustar[Ny+1][j] = -0.2*(ustar[Ny-2][j] - 5*ustar[Ny-1][j] + 15*ustar[Ny][j]); }

        FR(i,1,Ny){
            FR(j,1,Nx+1){
                vstar[i][j]    = -1*NS_pressurey(i,j,P,h);
                vstar[i][j]   += (1/sqrt(Gr))*NS_diffusiony(i,j,v,h);
                vstar[i][j]   += NS_buoyancy(i,j,T);
                H       	   = NS_convectiony(i,j,u,v,h);
                vstar[i][j]   -= 0.5*(3*H - Hnpy[i-1][j-1]);
                Hnpy[i-1][j-1] = H;
                // vstar = 1 / (1 + xhi Dt/Dtau)   *  (vn + dt * blabla + xhi Dt/Dtau v_s)
                // xsi calculé sur l'interieur du domaine
                vstar[i][j]    = 1/(1 + xiv[i][j] * dt/dtau) * (v[i][j] + dt * vstar[i][j]  + xiv[i][j] * dt/dtau * vmixer[i][j] );
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

        // computation of the average flux at the open boundary
        double avgflux = 0;
        FR(j,1,Nx+1){
            avgflux += l0 * h * (T[0][j] - T[1][j]) / h;
        }
        avg_heat_flux = avgflux / Nx;

        Tavg = 0.0;
        Trms = 0.0;
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                // TODO: fix this iteration problem
                // T[i+1][j] depends on the PREVIOUS value of T[i][j], not the NEXT value
                //  as is, this loop "co-updates" the values one after the other, depending on the iteration order.
                // this should be 4*h, not 2*h, as it compounds a mean (x+x)/2 and a derivative of step 2*h --> 4*h
                H       	   = (T[i][j+1]-T[i][j-1])*(u[i][j]+u[i][j-1])/(4*h) + (T[i-1][j] - T[i+1][j])*(v[i-1][j]+v[i][j])/(4*h);
                double Tn1     = -0.5*(3*H - HnpT[i-1][j-1]);
                HnpT[i-1][j-1] = H;
                Tn1           += (1/(Pr*sqrt(Gr)*h*h))*((T[i][j+1]-2*T[i][j]+T[i][j-1])+(T[i-1][j]-2*T[i][j]+T[i+1][j]));
                // T^(n+1) = 1/(1+xhi dt/Dtau) * (T^n + dt blabla + xhi *dt/dtau * Temp_src)
                // source term for the mixer is the average temperature of the mixer
                T_tmp[i][j] = 1/(1 + xiT[i][j] * dt/dtau) * (T[i][j] + dt*Tn1 + xiT[i][j] * dt/dtau * avg_temp_mixer );
                Tavg          += T_tmp[i][j];
            }
        }

        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                T[i][j] = T_tmp[i][j];
            }
        }
        Tavg /= Nx*Ny;
        FR(i,1,Ny+1) FR(j,1,Nx+1) Trms += (T[i][j] - Tavg)*(T[i][j] - Tavg);
        Trms /= Nx*Ny;
        FR(i,1,Ny+1){
            T[i][0] = T[i][1]; // Ghost points, adiabatic
            T[i][Nx+1] = T[i][Nx];
        }


        average_temp = Tavg;
        rms_temp = Trms;

        //============//
        // SOR Method //
        //============//
        globe = 1.0;
        int iter = 0;
    	SORtb = clock();
        while(globe > 1e-9){
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
            iter += 1;
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

        F(j,Nx+1) { u[0][j] = u[1][j]; // Ghost points, no-through flow
    				u[Ny+1][j] = -0.2*(u[Ny-2][j] - 5*u[Ny-1][j] + 15*u[Ny][j]); }
        F(i,Ny+1) { v[i][0]    = -0.2*(v[i][3]    - 5*v[i][2]    + 15*v[i][1]);
                    v[i][Nx+1] = -0.2*(v[i][Nx-2] - 5*v[i][Nx-1] + 15*v[i][Nx]); }

        //============//
        // Equation 4 //
        //============//
        F(i,Ny)
            F(j,Nx)
                P[i][j] += Phi[i+1][j+1];
        

        // normv
        F(i,Ny+1)
            F(j,Nx+1)
                normv[i][j] = sqrt(pow((u[i][j] + u[i+1][j])/2,2) + pow((v[i][j]+v[i][j+1])/2,2));

        // vorty
        F(i,Ny+1)
            F(j,Nx+1)
                vorty[i][j] = (v[i][j+1]-v[i][j])/h - (u[i][j]-u[i+1][j])/h;

        printf("temp mixer %.6f\n", avg_temp_mixer);
        printf("SOR Iterations = %d\n",iter);
        printf("Tavg = %.4lf\n",Tavg);
        printf("Trms = %.4lf\n",Trms);
        printf("%d\n",k);

        // analysis of the divergence maybe?
        printf("Divergence of [u*,v*]: %.16f\n", divergence(ustar, vstar, h, Nx, Ny));
        printf("Divergence of [u, v ]: %.16f\n", divergence(u, v, h, Nx, Ny));

        fwrite(&average_temp,sizeof(average_temp),1,fdiag);
        fwrite(&avg_heat_flux,sizeof(avg_heat_flux),1,fdiag);
        fwrite(&rms_temp,sizeof(rms_temp),1,fdiag);
        fwrite(&avg_temp_mixer,sizeof(avg_temp_mixer),1,fdiag);

        if ((k%saveIter) == 0) {
            if (k != 0){
                fclose(fTb);
                fclose(fvb);
                fclose(fwb);
            }
            sprintf(filename,"data/T_Nx%d_dt%d_iter%d_mixing%d.bin",Nx,(int)(1/dt),k,usemixer);
            fTb = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename,"data/v_Nx%d_dt%d_iter%d_mixing%d.bin",Nx,(int)(1/dt),k,usemixer);
            fvb = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename,"data/w_Nx%d_dt%d_iter%d_mixing%d.bin",Nx,(int)(1/dt),k,usemixer);
            fwb = fopen(filename,"wb");
            F(i,Ny+1){
                fwrite(T[i],sizeof(T[i][0]),Nx+2,fTb);
                fwrite(normv[i],sizeof(normv[i][0]),Nx+1,fvb);
                fwrite(vorty[i],sizeof(vorty[i][0]),Nx+1,fwb);
            }
            fwrite(T[Ny+1],sizeof(T[Ny+1][0]),Nx+2,fTb);
        }
        k++;
    }
    
	fclose(fTb);
	fclose(fvb);
    fclose(fwb);
    fclose(fdiag);


    /*
    FILE * fd1, *fd2, *fd3, *fd4;
    sprintf(filename, "data/average_temp");
    fd1 = fopen(filename, "wb");
    fwrite(average_temp, sizeof(average_temp[0]), miter, fd1);
    fclose(fd1);
    sprintf(filename, "data/avg_heat_flux");
    fd2 = fopen(filename, "wb");
    fwrite(avg_heat_flux, sizeof(avg_heat_flux[0]), miter, fd2);
    fclose(fd2);
    sprintf(filename, "data/rms_temp");
    fd3 = fopen(filename, "wb");
    fwrite(rms_temp, sizeof(rms_temp[0]), miter, fd3);
    fclose(fd3);
    sprintf(filename, "data/avg_temp_mixer");
    fd4 = fopen(filename, "wb");
    fwrite(avg_temp_mixer, sizeof(avg_temp_mixer[0]), miter, fd4);
    fclose(fd4);


    free_vector(average_temp);
    free_vector(rms_temp);
    free_vector(avg_heat_flux);
    free_vector(avg_temp_mixer);
    */
    free_matrix(P,Ny);
    free_matrix(T,Ny+2);
    free_matrix(T_tmp,Ny+2);
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
    free_matrix(umixer,Ny);
    free_matrix(vmixer,Ny);
    
	
	clock_t end = clock();
	printf("Total time %lf ms\n",(double)(end - begin) / CLOCKS_PER_SEC * 1000);
	printf("SOR time %lf ms (%.3lf of total time)\n",SORtime, SORtime/((double)(end-begin) / CLOCKS_PER_SEC * 1000));
}
