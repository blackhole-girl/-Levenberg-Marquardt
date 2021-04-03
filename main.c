#include<stdio.h>
#include<err.h>
#include<stdlib.h>
#include<math.h>
#include "fit_funcs.h"
#include "gaussj.h"
#include "nrutil.h"
#include "lm_func.h"



void fgauss(float x, float a[], float *y, float dyda[], int na);
int main(void)
{
	/*VARIABLES*/
	float *x = NULL, *y = NULL, *sig = NULL, *a= NULL;
	/*float dyda[4] ;*/
	float **covar = NULL, **alpha = NULL;
	FILE* data=NULL; /*data file*/
        char file[100] = "3c454_line.txt"; /*this is a pointer*/
	int ia[4] = {1,1,1,1};
	char lines[500];
	int i = 0, ndata = 0, na = 3, ma = 3, k = 0, n = 0, iter = 0;
	float alamda = 0. , chisq= 0. , ochisq = 0;
	



 	/*INSTRUCTIONS*/
	
	/*opening file to read data*/
        data = fopen("3c454_line.txt", "r");
	/*check*/
        if ( data == NULL)
        {
                printf("error in opening file %s\n", file);
                return 1;
        }

	/*number of lines*/
        /*to obtain the no of lines of the data file, assuming it doesnt surpass 500*/
        while(fgets(lines,500, data) != 0)
        {       
                ndata += 1;    
        }
        printf("%d\n",ndata);

	/*opening file again*/
        rewind(data);

	/*memory allocation*/
	/*memory allocation*/
        x = (float *)malloc((ndata+1) * sizeof (float));
        if ( x == NULL )
                err(-1, "Malloc failed in %s, line%d", __FILE__, __LINE__);
        y = (float *)malloc((ndata+1) * sizeof (float));
        if ( y == NULL )
                err(-1, "Malloc failed in %s, line%d", __FILE__, __LINE__);
        sig = (float *)malloc((ndata +1)* sizeof (float));
        if ( sig == NULL )
                err(-1, "Malloc failed in %s, line%d", __FILE__, __LINE__);
        a = (float *)malloc((ma+1)* sizeof (float));
        if ( a == NULL )
                err(-1, "Malloc failed in %s, line%d", __FILE__, __LINE__);
	/* for covar , alpha and chisq to use min function*/
	covar = (float**)malloc((ma+1)*sizeof(float*));
        if(covar==NULL){err(-1,"Malloc  failed  in %s,line %d",__FILE__ ,__LINE__ );}
        for(i=0;i<= ma+1;i++)
        {
                covar[i] = (float *)malloc ((ma+1)*sizeof(float));
                if(covar[i]==NULL){err(-1,"Malloc  failed  in %s,line %d",__FILE__ ,__LINE__ );}
        }

	alpha = (float**)malloc((ndata+1)*sizeof(float*));
	if(alpha==NULL){err(-1,"Malloc  failed  in %s,line %d",__FILE__ ,__LINE__ );}
	for(i=0;i<= ma+1;i++)
	{
		alpha[i] = (float *)malloc ((ma+1)*sizeof(float));
		if(alpha[i]==NULL){err(-1,"Malloc  failed  in %s,line %d",__FILE__ ,__LINE__ );}
	}	
	

        
	/*store data in arrays: x, y and sigma*/
	for (i = 1; i <= ndata ; i++)
        {
                fscanf(data, "%e %e %e ", &x[i], &y[i], &sig[i]);
	}

/*	display
	for (i = 1; i <= ndata ; i ++)
        {
                printf("x %e y %e sig  %e\n", x[i], y[i], sig[i]);
        }*/


	
	/* Reads in initial estimates for parameters 
	* chosen "educated" guess: [5, 5200, 20]*/
	printf("enter intial params of a \n");
	for(i=1;i<=na;i++) scanf("%f",&a[i]);
	
	/*display a  params
        for (i = 0; i < na ; i ++)
        {
                printf(" a matrix is %f\n", a[i]);
        }*/

	

	/*calling the guassian  outside main.... */

						       
	
	/*intialising lambda and calling mrqmin*/
	alamda = -1.0;
	ochisq = 0.0;
	mrqmin( x,  y, sig,  ndata,  a,  ia,  ma, covar, alpha,  &chisq, fgauss,  &alamda);
	
	/*ask user for no of interations*/
	printf("How many interations?\n");
        scanf("%d", &k);
	
	/*start loop*/
	for(i=1; i<=k; i++) 
{	iter++;
	ochisq = chisq;
	mrqmin( x,  y, sig,  ndata,  a,  ia,  ma, covar, alpha,  &chisq, fgauss,  &alamda);
	if (chisq > ochisq)										
		n = 0;
	else if (fabs(ochisq - chisq) < 0.001)
		n++;
	if(n > 100)
	{
		printf("fit failed\n");
		return -1;}		
}
	/*if good fit is found: set lam = 0 and call func again to get covariance*/
	alamda = 0.0;
	mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,&chisq,fgauss,&alamda);


	/*printing params a along with their associated errors using covar matrix*/
	printf("params of a are (Norm,centre,sigma) respectively:\n");
	for (i=1; i<= ma; i++) {printf("%f,",a[i]);}
	printf("\n");
	/*errors in a : sqrt of diag elements of the covar */
	printf("their respective errors are:\n");
	for (i=1; i<= ma; i++) {printf("%f,",sqrt( a[i] * covar[i][i]));}


	/*printing the chi squared result*/
	printf("\n");
	printf("chi squared test: %f\n", chisq);
	

										

	/*close file*/
        fclose(data);
        data = NULL;

	/*memory desallocation*/
	for( i=1 ; i<=ma ; i++ )
	{
		free( covar[i] );
		free ( alpha[i] );
		covar[i] = NULL;
		alpha[i] = NULL;
	}

	free(x);
	free(y);
	free(sig);
	free(a);
	

        x = NULL;
        y = NULL;
        sig = NULL;
	a = NULL;
	


	return 0;

}




