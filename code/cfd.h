void 	nrerror(char error_text[]);

double 	*vector(long N);
int 	*ivector(long N);
void 	free_vector(double *v);
void 	free_ivector(int *v);
void 	copy_vector(double *u, double *v, int N);

double 	**matrix(long NR, long NC);
void 	free_matrix(double **m, long NR);

double 	*linspace(double start, double end, long N);

void 	VectorToFile(FILE *f, double *u, int N);
void 	MatrixToFile(FILE *f, double **m, int NR, int NC);
