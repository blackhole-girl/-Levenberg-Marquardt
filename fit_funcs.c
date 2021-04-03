#include <math.h>
void fgauss(float x, float a[], float *y, float dyda[], int na)
/*y(x; a) is the sum of na/3 Gaussians (15.5.16). The amplitude, center, and width of the Gaussians are stored in consecutive locations of a: a[i] = Bk, a[i+1] = Ek, a[i+2] = Gk, k = 1,...,na/3. The dimensions of the arrays are a[1..na], dyda[1..na]. */
{
   int i;
   float fac,ex,arg;
   
   *y=0.0;
   for (i=1;i<=na-1;i+=3) {
       arg=(x-a[i+1])/a[i+2];
       ex=exp(-0.5*arg*arg);
       fac=a[i]*ex*arg;
       *y += a[i]*ex;
       dyda[i]=ex;
       dyda[i+1]=fac/a[i+2];
       dyda[i+2]=fac*arg/a[i+2];
   }
}
