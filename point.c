#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_qrng.h>
int main (int argc, char **argv)
{
int i,j;
int dim;
dim=atoi(argv[1]);
gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, dim);

	
FILE *fp;
 fp=fopen(argv[2],"w"); // Careful with dim
for (i = 0; i < 512; i++)
{
double v[dim];
gsl_qrng_get (q, v);
for (j=0;j<dim;j++){
//printf ("%.5f %.5f \n", v[0], v[1]);
		  fprintf(fp,"%lf ",v[j]);
	  
  }
  fprintf(fp,"\n");
  
}
fclose(fp);
gsl_qrng_free (q);
return 0;
}
