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
    int usemixer = 0;
	double H;
    // mixer radius
    double a = 3./10. * Nx * h;
    double omega = 0.1; // U/H = 1 right? TO CHECK 
    double L = h*Nx;
	
	if (argc == 4){ 
		Nx    = atoi(argv[1]);
		Ny    = 1.5*Nx;
		h     = 1.0/Ny;
		dt    = 1/(double)atoi(argv[2]);
		miter = atoi(argv[3]);
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
    // terme source du mixer
    double **umixer  = matrix(Ny,Nx);
    double **vmixer  = matrix(Ny,Nx);
    double **ustar   = matrix(Ny+2,Nx+1);
    double **vstar   = matrix(Ny+1,Nx+2);
    double **Hnpx    = matrix(Ny,Nx);
    double **Hnpy    = matrix(Ny,Nx);
    double **Phi     = matrix(Ny+2,Nx+2); 
    double **Phistar = matrix(Ny+2,Nx+2);
    double **R       = matrix(Ny,Nx);
    double **HnpT    = matrix(Ny,Nx);
    double **xhi     = matrix(Ny,Nx);
    double mixer_angle = 0;
    
    char filename[54];
	FILE *fPb,*fTb,*fub,*fvb,*fmix;
    FILE *debug_vstar, *debug_ustar;

    // diagnostics
    double * average_temp   = vector(miter);
    double * rms_temp       = vector(miter);
    double * avg_temp_mixer = vector(miter);
    double * avg_heat_flux  = vector(miter);
    
    /*********************************************************/
    /*DEBUG***************************************************/
    
    // Initialiser temperature : carré chaud puis reste à zéro            
    FR(i,2*Ny/5,3*Ny/5){
        FR(j,2*Nx/5,3*Nx/5)
        T[i][j] = 0.01; 
    }
        

    // Imposer u=1, v=0            
    FR(i,1,Ny+2){
        FR(j,1,Nx+1){
            u[i][j] = 1; 
    } }

    FR(i,1,Ny+1){
        FR(j,1,Nx+2){
            v[i][j] = 0; 
    } }

    

    
    /*********************************************************/
    /*********************************************************/

    F(k,miter){

        // compute the mixer mask
        double x, y, theta, d;
        double average_mixer_temp = 0;
        int nmixer_temp = 0;
        FR(j, 0, Nx){
            FR(i, 0, Ny){
                x = (j+1/2) * h; // 1/2 ? do we have to center the cell point
                y = (i+1/2) * h;
                theta = atan2( y - L, x - L/2 );
                d = sqrt( (x-L/2)*(x-L/2) + (y-L)*(y-L) );
                // compute mask
                xhi[i][j] = ( d <= a * cos(3*theta + mixer_angle) ) && usemixer;
                // ADD xhi_u and xhi_v for the speed (which is located on the borders of cells)
                //   (x != (j+1/2) * h, etc.)
                // speed source of the mixer
                // TO CHECK TO CHECK
                umixer[i][j] = -sin(theta) * omega * d;
                vmixer[i][j] =  cos(theta) * omega * d;
                // check the inside of the cylinder
                average_mixer_temp += T[i+1][j+1] * (d <= a);
                nmixer_temp += (d <= a);
            }
        } 

        // update the mixer angle
        mixer_angle += dt * omega;   // H = 1 right?
        
        printf("%d circle\n", nmixer_temp);
        avg_temp_mixer[k] = average_mixer_temp / nmixer_temp;
        printf("temp mixer %.6f\n", avg_temp_mixer[k]);  

        
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

        /*
        FR(i,1,Ny+1){
            FR(j,1,Nx){
                ustar[i][j]    = -1*NS_pressurex(i,j,P,h);
                ustar[i][j]   += (1/sqrt(Gr))*NS_diffusionx(i,j,u,h);
                H       	   = NS_convectionx(i,j,u,v,h);
                ustar[i][j]   -= 0.5*(3*H - Hnpx[i-1][j-1]);
                Hnpx[i-1][j-1] = H;
                // vstar = 1 / (1 + xhi Dt/Dtau)   *  (vn + dt * blabla + xhi Dt/Dtau v_s)
                // xsi calculé sur l'interieur du domaine
                ustar[i][j]    = 1/(1 + xhi[i-1][j-1] * dt/dtau) * (u[i][j] + dt * ustar[i][j]  + xhi[i-1][j-1] * dt/dtau * umixer[i-1][j-1] );
                // ustar[i][j]    = dt*ustar[i][j] + u[i][j];
            }
        }
        
        F(j,Nx+1) { ustar[0][j] = ustar[1][j]; // Ghost points, no-through flow
    				ustar[Ny+1][j] = -0.2*(ustar[Ny-2][j] - 5*ustar[Ny-1][j] + 15*ustar[Ny][j]); }
        //F(j,Nx+1) { ustar[0][j] = ustar[1][j]; ustar[Ny+1][j] = ustar[Ny][j]; } // WRONG, SEE CORRECTED ABOVE

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
                vstar[i][j]    = 1/(1 + xhi[i-1][j-1] * dt/dtau) * (v[i][j] + dt * vstar[i][j]  + xhi[i-1][j-1] * dt/dtau * vmixer[i-1][j-1] );
                // vstar[i][j]    = dt*vstar[i][j] + v[i][j];
            }
        }
        F(i,Ny+1) { vstar[i][0]    = -0.2*(vstar[i][3]    - 5*vstar[i][2]    + 15*vstar[i][1]);
                    vstar[i][Nx+1] = -0.2*(vstar[i][Nx-2] - 5*vstar[i][Nx-1] + 15*vstar[i][Nx]); }    */
                    
                    

                    

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
            // flux = - k dT/dy
            avgflux += l0 * h * (T[0][j] - T[1][j]) / h;
        }
        avg_heat_flux[k] = avgflux / Nx;

        Tavg = 0.0;
        Trms = 0.0;
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                // TODO: fix this iteration problem
                // T[i+1][j] depends on the PREVIOUS value of T[i][j], not the NEXT value
                //  as is, this loop "co-updates" the values one after the other, depending on the iteration order.

                // this should be 4*h, not 2*h, as it compounds a mean (x+x)/2 and a derivative of step 2*h --> 4*h
                H = 0.5*((u[i][j]*(T[i][j+1]-T[i][j])/h) + (u[i][j-1]*(T[i][j]-T[i][j-1])/h) + (v[i-1][j]*(T[i-1][j]-T[i][j])/h) + (v[i][j]*(T[i][j]-T[i+1][j])/h));
                //H = ((u[i][j-1] + u[i][j])/2 * (T[i][j+1] - T[i][j-1])/(2*h)) + ((v[i-1][j] + v[i][j])/2 * (T[i-1][j] - T[i+1][j])/(2*h)) + (1/2)*((T[i][j]*(u[i][j] - u[i][j-1])/h) + T[i][j]*(v[i-1][j] - v[i][j])/h);
                //H       	   = (T[i][j+1]-T[i][j-1])*(u[i][j]+u[i][j-1])/(4*h) + (T[i-1][j] - T[i+1][j])*(v[i-1][j]+v[i][j])/(4*h);
                
                // Idea : HnpT pas de valeurs au début de la boucle --> may cause problems 
                double Tn1     =   -0.5*(3*H - HnpT[i-1][j-1]); //= 0; 
                HnpT[i-1][j-1] = H;
                //Tn1           += (1/(Pr*sqrt(Gr)*h*h))*((T[i][j+1]-2*T[i][j]+T[i][j-1])+(T[i-1][j]-2*T[i][j]+T[i+1][j]));
                // T^(n+1) = 1/(1+xhi dt/Dtau) * (T^n + dt blabla + xhi *dt/dtau * Temp_src)
                // source term for the mixer is the average temperature of the mixer
                T_tmp[i][j] = 1/(1 + xhi[i-1][j-1] * dt/dtau) * (T[i][j] + dt*Tn1 + xhi[i-1][j-1] * dt/dtau * avg_temp_mixer[k] );
                //T_tmp[i][j]    = dt*Tn1 + T[i][j];
                Tavg          += T[i][j];
            }
        }

        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                T[i][j] = T_tmp[i][j];
            }
        }
        Tavg /= 1.5;
        FR(i,1,Ny+1) FR(j,1,Nx+1) Trms += (T[i][j] - Tavg)*(T[i][j] - Tavg);
        Trms /= Nx*Ny;
        FR(i,1,Ny+1){
            T[i][0] = T[i][1]; // Ghost points, adiabatic
            T[i][Nx+1] = T[i][Nx];
        }


        average_temp[k] = Tavg;
        rms_temp[k] = Trms;

        //============//
        // SOR Method //
        //============//
        globe = 1.0;
        int iter = 0;
    	SORtb = clock();
       /* while(globe > 1e-9){
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
                    R[i][j] = (Phi[i+2][j+1]-2*Phi[i+1][j+1]+Phi[i][j+1])/(h*h) + (Phi[i+1][j+2]-2*Phi[i+1][j+1]+Phi[i+1][j])/(h*h) - 1/dt*((ustar[i+1][j+1]-ustar[i+1][j])/h + (vstar[i][j+1]-vstar[i+1][j+1])/h);
                    sumR += R[i][j]*R[i][j];
                }
            }
            globe = dt*sqrt((3*h*h)/2*sumR);
            iter += 1;
            //printf("Sum local residual = %lf\n",sumR);
            //printf("Global error = %lf\n",globe);
        }  */
    	SORte = clock();
        SORtime += (double)(SORte - SORtb) / CLOCKS_PER_SEC * 1000;

        // residual of the method


        //============//
        // Equation 3 //
        //============//
        /* FR(i,1,Ny+1)
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
                    v[i][Nx+1] = -0.2*(v[i][Nx-2] - 5*v[i][Nx-1] + 15*v[i][Nx]); }   */
 
        //============//
        // Equation 4 //
        //============//
        F(i,Ny)
            F(j,Nx)
                P[i][j] += Phi[i+1][j+1];
        

        printf("SOR Iterations = %d\n",iter);
        printf("Tavg = %.4lf\n",Tavg);
        printf("Trms = %.4lf\n",Trms);
        printf("%d\n",k);

        // analysis of the divergence maybe?
        printf("Divergence of [u*,v*]: %.16f\n", divergence(ustar, vstar, h, Nx, Ny));
        printf("Divergence of [u, v ]: %.16f\n", divergence(u, v, h, Nx, Ny));

        if ((k%1) == 0) {
            if (k != 0){
                fclose(fPb);
                fclose(fTb);
                fclose(fub);
                fclose(fvb);
                fclose(fmix);
            }
            sprintf(filename,"data/P_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
    	    fPb = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename,"data/T_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            fTb = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename,"data/u_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            fub = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename,"data/v_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            fvb = fopen(filename,"wb");  // w for write, b for binary
            sprintf(filename, "data/mix_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            fmix = fopen(filename,"wb");
            sprintf(filename,"data/ustar_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            debug_ustar = fopen(filename, "wb");
            sprintf(filename,"data/vstar_Nx%d_dt%d_iter%d.bin",Nx,(int)(1/dt),k);
            debug_vstar = fopen(filename, "wb");
            F(i,Ny+2){
                fwrite(T[i],sizeof(T[i][0]),Nx+2,fTb);
                fwrite(u[i],sizeof(u[i][0]),Nx+1,fub);
                fwrite(ustar[i], sizeof(ustar[i][0]),Nx+1,debug_ustar);
                if(i < Ny+1){
                    fwrite(v[i],sizeof(T[i][0]),Nx+2,fvb);
                    fwrite(vstar[i],sizeof(T[i][0]),Nx+2,debug_vstar);
                }
            }
            F(i,Ny){
                fwrite(xhi[i], sizeof(xhi[0][0]), Nx, fmix);
            }
        }

    	//F(i,Ny)
    	//	fwrite(P[i],sizeof(P[i][0]),Nx,fPb);

        /*
    	F(i,Ny+2)
    		fwrite(u[i],sizeof(u[i][0]),Nx+1,fub);
    	F(i,Ny+1)
    		fwrite(v[i],sizeof(v[i][0]),Nx+2,fvb);
        */
        
    }
    
	fclose(fPb);
	fclose(fTb);
	fclose(fub);
	fclose(fvb);
    fclose(fmix);


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
