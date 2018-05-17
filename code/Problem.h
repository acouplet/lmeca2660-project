#ifndef PROBLEM_H
#define PROBLEM_H

#ifdef __cplusplus
    extern "C" {
#endif

typedef struct Problem {
    void(*Description)(struct Problem*);
    void(*BoundaryConditions)(struct Problem*);
    void(*Diagnostics)(struct Problem*);
    void(*WriteData)(struct Problem*);
    void(*IterationInfo)(struct Problem*);
    void(*Mixer)(struct Problem*);
    void(*Momentum)(struct Problem*);
    void(*Energy)(struct Problem*);
    void(*Poisson)(struct Problem*);
    void(*SpeedUpdate)(struct Problem*);
    void(*PressureUpdate)(struct Problem*);
    void(*Free)(struct Problem*);


    int Nx;
    int Ny;
    int useMixer;
    int saveIter;
    int iter;

    double h;
    double Pr;
    double dt;
    double dtau;
    double l0;
    double Tinf;
    double Gr;
    double alpha;
    double mixerRadius;
    double omega;
    double xg;
    double yg;
    double mixerAngle;
    double averageTemp;
    double rmsTemp;
    double avgTempMixer;
    double avgHeatFlux;

    double sumR;
    double globe;
    double Tavg;
    double Trm;
    clock_t SORtb;
    clock_t SORte;
    double SORtime;
    double **P;
    double **T;
    double **Ts;
    double **u;
    double **v;
    double **normv;
    double **vorty;
    double **umixer;
    double **vmixer;
    double **ustar;
    double **vstar;
    double **xiu;
    double **xiv;
    double **xiT;
    double **Hnpu;
    double **Hnpv;
    double **HnpT;
    double **Phi;
    double **Phistar;
    double **R;

    FILE *fdiag;

} Problem;

Problem* New_Problem(int Nx, int ddt, int saveIter, int useMixer);
Problem* Default_Problem();
void Problem_Description(Problem*);
void Problem_BoundaryConditions(Problem*);
void Problem_Mixer(Problem*);
void Problem_Diagnostics(Problem*);
void Problem_WriteData(Problem*);
void Problem_IterationInfo(Problem*);
void Problem_Momentum(Problem*);
void Problem_Energy(Problem*);
void Problem_Poisson(Problem*);
void Problem_SpeedUpdate(Problem*);
void Problem_PressureUpdate(Problem*);
void Problem_Free(Problem*);
    
#ifdef __cplusplus
}
#endif

#endif
