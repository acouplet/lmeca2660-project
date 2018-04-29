void 	nrerror(char error_text[]);

double 	*vector(long N);
int 	*ivector(long N);
void 	free_vector(double *v);
void 	free_ivector(int *v);
void 	copy_vector(double *u, double *v, int N);

double 	**matrix(long NR, long NC);
void 	free_matrix(double **m, long NR);
void    copy_matrix(double **m, double **o, long NR, long NC);

double 	*linspace(double start, double end, long N);

void 	VectorToFile(FILE *f, double *u, int N);
void 	MatrixToFile(FILE *f, double **m, int NR, int NC);

double  NS_pressurex(int i,int j,double **P,double h);
double  NS_pressurey(int i,int j,double **P,double h);
double  NS_diffusionx(int i, int j, double **u, double h);
double  NS_diffusiony(int i, int j, double **v, double h);
double  NS_buoyancy(int i, int j, double **T);
double  NS_convectionx(int i, int j, double **u, double **v, double h);
double  NS_convectiony(int i, int j, double **u, double **v, double h);
