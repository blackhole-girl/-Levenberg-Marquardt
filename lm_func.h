void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **covar, float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [], int), float *alamda);
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int));
void covsrt(float **covar, int ma, int ia[], int mfit);



