#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
	
	
