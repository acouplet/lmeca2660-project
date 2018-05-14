#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void nrerror(char error_text[]){
	/* ADAPTED FROM NUMERICAL RECIPES IN C
	 * CFD error handler 
	 */
	fprintf(stderr,"** CFD run-time error...\n");
	fprintf(stderr,"** %s\n",error_text);
	fprintf(stderr,"** now exiting to system...\n");
	exit(1);
}

double *vector(long N){
	/* 
	 * allocate a double vector with subscript range v[0..N-1] 
	 */
	double *v;

	v = (double *)calloc(N, sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v;
}

int *ivector(long N){
	/* 
	 * allocate a double vector with subscript range v[0..N-1] 
	 */
	int *v;

	v = (int *)calloc(N, sizeof(int));
	if (!v) nrerror("allocation failure in vector()");
	return v;
}

void free_vector(double *v){
	free(v);
}

void free_ivector(int *v){
	free(v);
}

void copy_vector(double *u, double *v, int N){
	int i;
	for (i=0; i<N; i++){
		u[i] = v[i];
	}
}
	

double **matrix(long NR, long NC){
	/* 
	 * allocate a double matrix with subscript range m[0..NR-1][0..NC-1]
	 * Could be optimized https://stackoverflow.com/questions/2128728/allocate-matrix-in-c
	 */
	long i;
	double **m;

	m = (double **)calloc(NR, sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	
	for(i=0; i<NR; i++){
		m[i] = (double *)calloc(NC,sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
	}
	
	return m;
}

void free_matrix(double **m, long NR){
	long i;
	for(i=0; i<NR; i++) free(m[i]);
	free(m);
}

void copy_matrix(double **m, double **o, long NR, long NC){
    long i,j;
    for(i=0;i<NR;i++)
        for(j=0;j<NC;j++)
            o[i][j] = m[i][j];
}



double *linspace(double start, double end, long N){
	int i;
	double *v, delta;
	v = vector(N);
	delta = (end-start)/(N-1);
	for (i=0; i<N; i++){
		v[i] = start + i*delta;
	}
	return v;
}

void VectorToFile(FILE *f, double *u, int N){
	int i;
	for(i=0;i<N;i++) fprintf(f,"%lf\n",u[i]);
}

void MatrixToFile(FILE *f, double **m, int NR, int NC){
	int i,j;
	for(i=0;i<NR;i++) {
		for(j=0;j<NC-1;j++)
			fprintf(f,"%lf,",m[i][j]);
		fprintf(f,"%lf\n",m[i][j]);
	}
}

double NS_pressurex(int i,int j,double **P,double h){
    return (P[i-1][j]-P[i-1][j-1])/h;
}

double NS_pressurey(int i, int j, double **P, double h){
    return (P[i-1][j-1]-P[i][j-1])/h;
}


double NS_diffusionx(int i, int j, double **u, double h){
    return (u[i][j+1]-2*u[i][j]+u[i][j-1])/(h*h) + (u[i-1][j]-2*u[i][j]+u[i+1][j])/(h*h);
}

double NS_diffusiony(int i, int j, double **v, double h){
    return (v[i][j+1]-2*v[i][j]+v[i][j-1])/(h*h) + (v[i-1][j]-2*v[i][j]+v[i+1][j])/(h*h);
}

double NS_buoyancy(int i, int j, double **T){
    return (T[i][j] + T[i+1][j])/2;
}

double NS_convectionx(int i, int j, double **u, double **v, double h){
	// shouldn't it be (u[i][j]-u[i][j-1]) instead ? (and same later)
    double adv1 = (1/(4*h))*((u[i][j-1]+u[i][j])*(u[i][j]-u[i][j-1]) + (u[i][j]+u[i][j+1])*(u[i][j+1]-u[i][j]));
    double adv2 = (1/(4*h))*((v[i-1][j+1]+v[i-1][j])*(u[i-1][j]-u[i][j]) + (v[i][j+1]+v[i][j])*(u[i][j]-u[i+1][j]));
    double div1 = (1/(4*h))*((u[i][j+1]+u[i][j])*(u[i][j+1]+u[i][j]) - (u[i][j]+u[i][j-1])*(u[i][j]+u[i][j-1]));
    double div2 = (1/(4*h))*((v[i-1][j+1]+v[i-1][j])*(u[i-1][j]+u[i][j]) - (v[i][j+1]+v[i][j])*(u[i][j]+u[i+1][j]));
    return 0.5*(adv1+adv2) + 0.5*(div1+div2);
}

double NS_convectiony(int i, int j, double **u, double **v, double h){
	// error: (v[i][j]-v[i-1][j]))    should be    (v[i][j] - v[i][j-1])
    double adv1 = (1/(4*h))*((u[i][j]+u[i+1][j])*(v[i][j+1]-v[i][j]) + (u[i][j-1]+u[i+1][j-1])*(v[i][j]-v[i][j-1]));
    double adv2 = (1/(4*h))*((v[i-1][j]+v[i][j])*(v[i-1][j]-v[i][j]) + (v[i][j]+v[i+1][j])*(v[i][j]-v[i+1][j]));
    double div1 = (1/(4*h))*((u[i][j]+u[i+1][j])*(v[i][j+1]+v[i][j]) - (u[i][j-1]+u[i+1][j-1])*(v[i][j]+v[i][j-1]));
    double div2 = (1/(4*h))*((v[i-1][j]+v[i][j])*(v[i-1][j]+v[i][j]) - (v[i][j]+v[i+1][j])*(v[i][j]+v[i+1][j]));
    return 0.5*(adv1+adv2) + 0.5*(div1+div2);
}

int xi(double x, double y, double a, double mixerAngle) {
    double xg = 1.0/3.0;
    double yg = 1.0/3.0;
    double theta = atan2(y-yg,x-xg);
    double d = sqrt(pow(x-xg,2) + pow(y-yg,2));
    return (d <= a*cos(3*theta + mixerAngle));
}

	
// TODO: recheck divergence form
