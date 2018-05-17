#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "cfd.h"

#include "Problem.h"

#define FR(i,a,b) for (int i=(a); i<(b); i++)
#define F(i,n) FR(i,0,n)
#define DR(i,a,b) for (int i=(b); i-->(a);)
#define D(i,n) DR(i,0,n)

static void Problem_Init(Problem*, int Nx, int ddt, int saveIter, int useMixer);

Problem* Default_Problem(){
    return New_Problem(200,100,20,0);
}

Problem* New_Problem(int Nx, int ddt, int saveIter, int useMixer){
    Problem *This = malloc(sizeof(Problem));
    if(!This) return NULL;
    Problem_Init(This, Nx, ddt, saveIter, useMixer);
    This->Free = Problem_Free;
    return This;
}

static void Problem_Init(Problem *This, int Nx, int ddt, int saveIter, int useMixer){
    This->Description = Problem_Description;
    This->BoundaryConditions = Problem_BoundaryConditions;
    This->Diagnostics = Problem_Diagnostics;
    This->WriteData = Problem_WriteData;
    This->IterationInfo = Problem_IterationInfo;
    This->Momentum = Problem_Momentum;
    This->Energy = Problem_Energy;
    This->Poisson = Problem_Poisson;
    This->SpeedUpdate = Problem_SpeedUpdate;
    This->PressureUpdate = Problem_PressureUpdate;
    This->Mixer = Problem_Mixer;

    This->Nx          = Nx;
    This->Ny          = 1.5*Nx;
    This->useMixer    = useMixer;
    This->saveIter    = saveIter;
    This->iter        = 0;
    This->h           = 1.0/This->Ny;
    This->Pr          = 2.0;
    This->dt          = 1.0/ddt;
    This->dtau        = This->dt/1e3;
    This->l0          = 1e-3;
    This->Tinf        = -5e-3;
    This->Gr          = 2.0e10;
    This->alpha       = 1.97;
    This->mixerRadius = 0.3*This->Nx*This->h;
    This->omega       = 0.1;
    This->xg          = 1./3.;
    This->yg          = 1./3.;
    This->mixerAngle  = 0;

    This->sumR         = 0.0;
    This->globe        = 0.0;
    This->averageTemp  = 0.0;
    This->rmsTemp      = 0.0;
    This->avgTempMixer = 0.0;
    This->avgHeatFlux  = 0.0;
    This->SORtime = 0.0;

    int Ny = This->Ny;
    This->P = matrix(Ny,Nx);
    This->T = matrix(Ny+2,Nx+2);
    This->Ts = matrix(Ny+2,Nx+2);
    This->u = matrix(Ny+2,Nx+1);
    This->v = matrix(Ny+1,Nx+2);
    This->normv = matrix(Ny+1,Nx+1);
    This->vorty = matrix(Ny+1,Nx+1);
    This->umixer = matrix(Ny+2,Nx+1);
    This->vmixer = matrix(Ny+1,Nx+2);
    This->ustar = matrix(Ny+2,Nx+1);
    This->vstar = matrix(Ny+1,Nx+2);
    This->xiu = matrix(Ny+2,Nx+1);
    This->xiv = matrix(Ny+1,Nx+2);
    This->xiT = matrix(Ny+2,Nx+2);
    This->Hnpu = matrix(Ny+2,Nx+1);
    This->Hnpv = matrix(Ny+1,Nx+2);
    This->HnpT = matrix(Ny+2,Nx+2);
    This->Phi = matrix(Ny+2,Nx+2);
    This->Phistar = matrix(Ny+2,Nx+2);
    This->R = matrix(Ny,Nx);

    char filename[100];
    char dirname[100];

    struct stat st = {0};
    sprintf(dirname,"data/Nx%d_dt%d_mixing%d",Nx,(int)(1/This->dt),This->useMixer);
    if (stat(dirname, &st) == -1) {
        mkdir(dirname, 0700);
    }

    sprintf(filename,"data/Nx%d_dt%d_mixing%d/diagnostics.bin",Nx,(int)(1/This->dt),This->useMixer);
    This->fdiag = fopen(filename,"wb");
}

void Problem_Momentum(Problem *This){
    double **P = This->P;
    double **u = This->u;
    double **v = This->v;
    int Nx = This->Nx;
    int Ny = This->Ny;
    double h = This->h;
    double H, adv1, adv2, div1, div2;
    double pressure, diffusion, convection, buoyancy, penalization;

    // ustar
    FR(i,1,Ny+1){
        FR(j,1,Nx){
            pressure = (P[i-1][j] - P[i-1][j-1])/h;
            diffusion = (u[i][j+1]-2*u[i][j]+u[i][j-1])/(h*h) + (u[i-1][j]-2*u[i][j]+u[i+1][j])/(h*h);
            adv1 = (1/(4*h))*((u[i][j-1]+u[i][j])*(u[i][j]-u[i][j-1]) + (u[i][j]+u[i][j+1])*(u[i][j+1]-u[i][j]));
            adv2 = (1/(4*h))*((v[i-1][j+1]+v[i-1][j])*(u[i-1][j]-u[i][j]) + (v[i][j+1]+v[i][j])*(u[i][j]-u[i+1][j]));
            div1 = (1/(4*h))*((u[i][j+1]+u[i][j])*(u[i][j+1]+u[i][j]) - (u[i][j]+u[i][j-1])*(u[i][j]+u[i][j-1]));
            div2 = (1/(4*h))*((v[i-1][j+1]+v[i-1][j])*(u[i-1][j]+u[i][j]) - (v[i][j+1]+v[i][j])*(u[i][j]+u[i+1][j]));
            H = 0.5*(adv1+adv2) + 0.5*(div1+div2);
            convection = 0.5*(3*H - This->Hnpu[i][j]);
            penalization = This->xiu[i][j]*This->umixer[i][j]/This->dtau + u[i][j]/This->dt;
            This->ustar[i][j] = 1./(1./This->dt + This->xiu[i][j]/This->dtau)*((1/sqrt(This->Gr))*diffusion - pressure - convection + penalization);
            This->Hnpu[i][j] = H;
        }
    }
    F(j,Nx+1) { 
        This->ustar[0][j] = This->ustar[1][j]; // top, no-through flow
        This->ustar[Ny+1][j] = -0.2*(This->ustar[Ny-2][j] - 5*This->ustar[Ny-1][j] + 15*This->ustar[Ny][j]); // bottom, no-slip
    }

    // vstar
    FR(i,1,Ny){
        FR(j,1,Nx+1){
            pressure = (P[i-1][j-1]-P[i][j-1])/h;
            diffusion = (v[i][j+1]-2*v[i][j]+v[i][j-1])/(h*h) + (v[i-1][j]-2*v[i][j]+v[i+1][j])/(h*h);
            adv1 = (1/(4*h))*((u[i][j]+u[i+1][j])*(v[i][j+1]-v[i][j]) + (u[i][j-1]+u[i+1][j-1])*(v[i][j]-v[i][j-1]));
            adv2 = (1/(4*h))*((v[i-1][j]+v[i][j])*(v[i-1][j]-v[i][j]) + (v[i][j]+v[i+1][j])*(v[i][j]-v[i+1][j]));
            div1 = (1/(4*h))*((u[i][j]+u[i+1][j])*(v[i][j+1]+v[i][j]) - (u[i][j-1]+u[i+1][j-1])*(v[i][j]+v[i][j-1]));
            div2 = (1/(4*h))*((v[i-1][j]+v[i][j])*(v[i-1][j]+v[i][j]) - (v[i][j]+v[i+1][j])*(v[i][j]+v[i+1][j]));
            H = 0.5*(adv1+adv2) + 0.5*(div1+div2);
            convection = 0.5*(3*H - This->Hnpv[i][j]);
            buoyancy = (This->T[i][j] + This->T[i+1][j])/2;
            penalization = This->xiv[i][j]*This->vmixer[i][j]/This->dtau + v[i][j]/This->dt;
            This->vstar[i][j] = 1./(1./This->dt + This->xiv[i][j]/This->dtau)*((1/sqrt(This->Gr))*diffusion - pressure - convection + buoyancy + penalization);
            This->Hnpv[i][j] = H;
        }
    }
    F(i,Ny+1) { 
        This->vstar[i][0]    = -0.2*(This->vstar[i][3]    - 5*This->vstar[i][2]    + 15*This->vstar[i][1]);
        This->vstar[i][Nx+1] = -0.2*(This->vstar[i][Nx-2] - 5*This->vstar[i][Nx-1] + 15*This->vstar[i][Nx]); 
    }
}

void Problem_Energy(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;
    double h = This->h;
    double **T = This->T;
    double H, convection, diffusion, penalization;
    This->avgHeatFlux = 0;
    FR(j,1,Nx+1)
        This->avgHeatFlux += This->l0*h*(T[0][j] - T[1][j])/h;
    This->avgHeatFlux /= Nx;

    This->averageTemp = 0;
    This->rmsTemp = 0;
    FR(i,1,Ny+1){
        FR(j,1,Nx+1){
            H = 0.5*((This->u[i][j]*(T[i][j+1]-T[i][j])/h) + (This->u[i][j-1]*(T[i][j]-T[i][j-1])/h) + (This->v[i-1][j]*(T[i-1][j]-T[i][j])/h) + (This->v[i][j]*(T[i][j]-T[i+1][j])/h));
            convection = 0.5*(3*H - This->HnpT[i][j]);
            diffusion = (1/(This->Pr*sqrt(This->Gr)*h*h))*((T[i][j+1]-2*T[i][j]+T[i][j-1])+(T[i-1][j]-2*T[i][j]+T[i+1][j]));
            penalization = This->xiT[i][j]*(This->avgTempMixer/This->dtau) + T[i][j]/This->dt;
            This->Ts[i][j] = 1./(1./This->dt + This->xiT[i][j]/This->dtau)*(diffusion - convection + penalization);
            This->averageTemp += This->Ts[i][j];
            This->HnpT[i][j] = H;
        }
    }

    FR(i,1,Ny+1){
        FR(j,1,Nx+1){
            T[i][j] = This->Ts[i][j];
        }
    }
    This->averageTemp /= (Nx*Ny);
    FR(i,1,Ny+1) 
        FR(j,1,Nx+1) 
            This->rmsTemp += pow(T[i][j] - This->averageTemp,2);
    This->rmsTemp /= (Nx*Ny);
}

void Problem_Poisson(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;
    double h = This->h;
    double **Phistar = This->Phistar;
    double **Phi = This->Phi;
    double **ustar = This->ustar;
    double **vstar = This->vstar;
    This->globe = 1.0;
    int iter = 0;
    This->SORtb = clock();
    while(This->globe > 1e-9){
        FR(i,1,Ny+1){
            FR(j,1,Nx+1){
                Phistar[i][j] = (h*h)/4*(-1/This->dt*((ustar[i][j]-ustar[i][j-1])/h + (vstar[i-1][j]-vstar[i][j])/h) + (Phi[i+1][j] + Phi[i-1][j])/(h*h) + (Phi[i][j+1] + Phi[i][j-1])/(h*h));
                Phi[i][j] = This->alpha*Phistar[i][j] + (1-This->alpha)*Phi[i][j];
            }
        }
        This->sumR = 0;
        F(i,Ny){
            F(j,Nx){
                This->R[i][j] = (Phi[i+2][j+1]-2*Phi[i+1][j+1]+Phi[i][j+1])/(h*h) + (Phi[i+1][j+2]-2*Phi[i+1][j+1]+Phi[i+1][j])/(h*h) - 1/This->dt*((ustar[i+1][j+1]-ustar[i+1][j])/h + (vstar[i][j+1]-vstar[i+1][j+1])/h);
                This->sumR += pow(This->R[i][j],2);
            }
        }
        This->globe = This->dt*sqrt((3*h*h)/2*This->sumR);
        iter++;
    }
    This->SORte = clock();
    This->SORtime += (double)(This->SORte - This->SORtb)/CLOCKS_PER_SEC*1000;
}

void Problem_SpeedUpdate(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;
    double **Phi = This->Phi;
    double h = This->h;

    FR(i,1,Ny+1)
        FR(j,1,Nx)
            This->u[i][j] = -(This->dt/h)*(Phi[i][j+1]-Phi[i][j]) + This->ustar[i][j];

    FR(i,1,Ny)
        FR(j,1,Nx+1)
            This->v[i][j] = -(This->dt/h)*(Phi[i][j]-Phi[i+1][j]) + This->vstar[i][j];
}

void Problem_PressureUpdate(Problem *This){
    F(i,This->Ny)
        F(j,This->Nx)
            This->P[i][j] += This->Phi[i+1][j+1];
}

void Problem_Diagnostics(Problem *This){
    This->avgHeatFlux = 0;
    int Nx = This->Nx;
    int Ny = This->Ny;
    FR(j,1,Nx+1){
        This->avgHeatFlux += This->l0 * This->h * (This->T[0][j] - This->T[1][j]) / This->h;
    }
    This->avgHeatFlux = This->avgHeatFlux / Nx;

    This->averageTemp = 0;
    FR(i,1,Ny+1)
        FR(j,1,Nx+1)
            This->averageTemp += This->T[i][j];
    This->averageTemp /= (Nx*Ny); 

    This->rmsTemp = 0;
    FR(i,1,Ny+1)
        FR(j,1,Nx+1)
            This->rmsTemp += pow(This->T[i][j]-This->averageTemp,2);
    This->rmsTemp /= (Nx*Ny); 
}

void Problem_WriteData(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;
    int iter = This->iter;
    double **u = This->u;
    double **v = This->v;
    double h = This->h;
    char filename[100];

    fwrite(&(This->averageTemp),sizeof(This->averageTemp),1,This->fdiag);
    fwrite(&(This->avgHeatFlux),sizeof(This->avgHeatFlux),1,This->fdiag);
    fwrite(&(This->rmsTemp),sizeof(This->averageTemp),1,This->fdiag);
    fwrite(&(This->avgTempMixer),sizeof(This->avgTempMixer),1,This->fdiag);


    if (iter%This->saveIter == 0){
        // normv & vorty
        F(i,Ny+1){
            F(j,Nx+1){
                This->normv[i][j] = sqrt(pow((u[i][j] + u[i+1][j])/2,2) + pow((v[i][j]+v[i][j+1])/2,2));
                This->vorty[i][j] = (v[i][j+1]-v[i][j])/h - (u[i][j]-u[i+1][j])/h;
            }
        }

        sprintf(filename,"data/Nx%d_dt%d_mixing%d/T_iter%d.bin",Nx,(int)(1/This->dt),This->useMixer,iter);
        FILE *fT = fopen(filename,"wb");
        sprintf(filename,"data/Nx%d_dt%d_mixing%d/v_iter%d.bin",Nx,(int)(1/This->dt),This->useMixer,iter);
        FILE *fv = fopen(filename,"wb");
        sprintf(filename,"data/Nx%d_dt%d_mixing%d/w_iter%d.bin",Nx,(int)(1/This->dt),This->useMixer,iter);
        FILE *fw = fopen(filename,"wb");

        F(i,Ny+1){
            fwrite(This->T[i],sizeof(This->T[i][0]),Nx+2,fT);
            fwrite(This->normv[i],sizeof(This->normv[i][0]),Nx+1,fv);
            fwrite(This->vorty[i],sizeof(This->vorty[i][0]),Nx+1,fw);
        }
        fwrite(This->T[Ny+1],sizeof(This->T[Ny+1][0]),Nx+2,fT);

        fclose(fT);
        fclose(fv);
        fclose(fw);
    }
}

void Problem_BoundaryConditions(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;

    // T BC
    FR(j,1,Nx+1){
        This->T[Ny+1][j] = This->T[Ny][j] + This->h; // bottom heat transfer
        This->T[0][j] = (-This->h/This->l0)*(1.5*This->T[1][j] - 0.5*This->T[2][j] - This->Tinf) + This->T[1][j]; // top heat transfer
    }
    FR(i,1,Ny+1){
        This->T[i][0] = This->T[i][1]; // left adiabatic wall
        This->T[i][Nx+1] = This->T[i][Nx]; // right adiabatic wall
    }

    // u BC
    F(j,Nx+1) { 
        This->u[0][j] = This->u[1][j]; // top no-through flow
        This->u[Ny+1][j] = -0.2*(This->u[Ny-2][j] - 5*This->u[Ny-1][j] + 15*This->u[Ny][j]);  // bottom, no-slip
    }

    // v BC
    F(i,Ny+1) { 
        This->v[i][0]    = -0.2*(This->v[i][3]    - 5*This->v[i][2]    + 15*This->v[i][1]); // left no-slip
        This->v[i][Nx+1] = -0.2*(This->v[i][Nx-2] - 5*This->v[i][Nx-1] + 15*This->v[i][Nx]); // right no-slip
    }
}

void Problem_Mixer(Problem *This){
    int Nx = This->Nx;
    int Ny = This->Ny;
    double x, y, theta, d;
    double h = This->h;

    // Xi mask on u nodes
    FR(i,1,Ny+1){
        FR(j,1,Nx){
            x = j*h;
            y = (Ny-i+0.5)*h;
            theta = atan2(y-This->yg,x-This->xg);
            d = sqrt(pow(x-This->xg,2) + pow(y-This->yg,2));
            This->xiu[i][j] = (d<=This->mixerRadius*(cos(3*theta + This->mixerAngle))) && This->useMixer;
            This->umixer[i][j] = -sin(theta)*This->omega*d;
        }
    }

    // Xi mask on v nodes
    int nMixerNodes = 0;
    FR(i,1,Ny){
        FR(j,1,Nx+1){
            x = (j-0.5)*h;
            y = (Ny-i)*h;
            theta = atan2(y-This->yg,x-This->xg);
            d = sqrt(pow(x-This->xg,2) + pow(y-This->yg,2));
            This->xiv[i][j] = (d<=This->mixerRadius*(cos(3*theta + This->mixerAngle))) && This->useMixer;
            This->vmixer[i][j] = cos(theta)*This->omega*d;

        }
    }

    // Xi mask on T nodes + avgTempMixer
    FR(i,1,Ny+1){
        FR(j,1,Nx+1){
            x = (j-0.5)*h;
            y = (Ny-i+0.5)*h;
            theta = atan2(y-This->yg,x-This->xg);
            d = sqrt(pow(x-This->xg,2) + pow(y-This->yg,2));
            This->xiT[i][j] = (d<=This->mixerRadius*cos(3*theta + This->mixerAngle)) && This->useMixer;
            if (d <= This->mixerRadius){
                This->avgTempMixer = This->T[i][j];
                nMixerNodes++;
            }
        }
    }
    This->mixerAngle += This->dt*This->omega;
    This->avgTempMixer /= nMixerNodes;
}

void Problem_IterationInfo(Problem *This){
    printf("Iteration: %d\tAverage temperature: %lf\tRMS temperature: %lf\tAverage mixer temperature %lf\tAverage heat flux %lf\n",This->iter,This->averageTemp,This->rmsTemp,This->avgTempMixer,This->avgHeatFlux);
}

void Problem_Free(Problem *This){
    // WE NEED TO COMPLETE THIS FUNCTION
    fclose(This->fdiag);
    free(This);
}

void Problem_Description(Problem *This){
    printf("--- Problem Description ---\n");
    printf("Domain size:\t(%dx%d) \t| Spatial step:\t%.5f\t| Time step:\t%.3f\n",This->Ny,This->Nx,This->h,This->dt);
    printf("save every:\t%d iterations\t| use of mixer:\t%d\t|\n",This->saveIter,This->useMixer);
    printf("\n");
}
