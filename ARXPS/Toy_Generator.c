// This code produces the toy depth profile data, including
// a the data with relative error fixed by sigma, the kernel,
// and a model (and an equivalent initial image) with a variable
// cutoff.
//
// to compile gcc Toy_generator.c -lm
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//
// some fundemental parameters
//
#define PI 3.1415927
//
// some parameters which define our problem
//
#define NDATA 6		//the number of toy data

double rang(void);

int main()
{
	int iseed,i,j,nf,ndata;
	double data[NDATA],r1,r2,cutoff,dz0,beta,sigma,lambda;
//	double theta[NDATA]={70, 63, 55, 40, 10, 0};
	double theta[NDATA]={77, 70, 63, 55, 40, 10};
	
// declare some pointers
	FILE *ModelFile, *DataFile, *InitFile, *ToyfFile;
	double *z,*dz,*kernel,*model,*initimage,*f;
	
	ModelFile= fopen("model","w");
	DataFile= fopen("kerneldata","w");
	InitFile= fopen("initimage","w");
	ToyfFile= fopen("toyf","w");

	printf("Enter a random number seed (int): ");
	scanf("%d",&iseed);
	srand(iseed); // seed the random number generator
	printf("Enter the cutoff for the model: ");
	scanf("%lf",&cutoff);
	printf("Enter the dz, the depth step: ");
	scanf("%lf",&dz0);
	printf("Enter the nf, the number of depths: ");
	scanf("%d",&nf);
	printf("Enter beta, the inverse width of f: ");
	scanf("%lf",&beta);
	printf("Enter sigma, the relative error of data: ");
	scanf("%lf",&sigma);
	printf("Enter lambda, the penetration depth: ");
	scanf("%lf",&lambda);
	
// allocate some memory

	f = 		(double *) 	malloc(sizeof(double)*nf);
	z = 		(double *) 	malloc(sizeof(double)*nf);
	dz = 	(double *) 	malloc(sizeof(double)*nf);
	model = 	(double *) 	malloc(sizeof(double)*nf);
	initimage = 	(double *) 	malloc(sizeof(double)*nf);
	f = 		(double *) 	malloc(sizeof(double)*nf);
	kernel = 	(double *) 	malloc(sizeof(double)*nf*NDATA);
	
// Test the Gaussian RNG
	r1=0.0;
	for(j=0;j<100000;j++)
	{
	r1=r1+pow(rang(),2);
	}
	r1=r1/(double)100000;
	printf("<R^2>=%lf",r1);

//convert the angles into radians
	for(j=0;j<NDATA;j++)
	{
//		printf("theta[j]=%lf",theta[j]);
		theta[j]=theta[j]*PI/(double)180;
//		printf("  %lf\n",theta[j]);
	}

// Calculate the integration weight, dz[i]

	for(i=1;i<nf-1;i++)
	{
	dz[i]=dz0;
	}
	dz[0]=dz0/2.0;
	dz[nf-1]=dz[0];
// Calculate the toy depth distribution	
	for (i=0;i<nf;i++)
	{
		z[i]=((float)i)*dz0;
// sigmoid		
//
//		  if(z[i]>1.0)
//		  {
//		    f[i]=dz[i]*exp(-beta*(z[i]-1.0))/(exp(-beta*(z[i]-1.0))+1.0);
//		  }
//		  else
//		  {
//		    f[i]=dz[i]*1.0/(exp(beta*(z[i]-1.0))+1.0);
//		  }
// gaussian
//		f[i]=dz[i]*exp(-pow(beta*(z[i]-1.0),2));
// exponential
		f[i]=dz[i]*exp(-beta*z[i]);
//
	fprintf(ToyfFile,"%lf,  %lf\n",z[i],f[i]/dz[i]);
	}
	
	
// calculate the kernel
	for(i=0;i<nf;i++)
	{
		for(j=0;j<NDATA;j++)
		{
			kernel[i+j*nf]=exp(-z[i]/(lambda*cos(theta[j])));
		}
	}
	
// calculate the pure data from f and the kernel
	for(j=0;j<NDATA;j++)
	{
		data[j]=0.0;
		for(i=0;i<nf;i++)
		{
			data[j]+=f[i]*kernel[i+j*nf];
		}
	}

//
// Now start writing thigs to the files
//
// first the header
	i=NDATA; 
	j=nf;
	fprintf(DataFile,"%d , %d , 1\n",i,j);
//
// now the data and its error
//
	for(j=0;j<NDATA;j++)
	{
		r1=(1.0+sigma*rang())*data[j] ; 
		r2=sigma*data[j];
		fprintf(DataFile,"%lf , %lf \n",r1,r2);
	}
//
//now write the kernel
//
	for(j=0;j<NDATA;j++)
	{
		for(i=0;i<nf;i++)
		{
			fprintf(DataFile,"%lf\n",kernel[i+j*nf]);
		}
	}

//
//now form and write the model
//
	for(i=0;i<nf;i++)
	{
//triangle
		model[i]=2.0*(cutoff-z[i])/(cutoff*cutoff);
// flat
//                model[i]=1.0/cutoff;
//box
//              if((z[i]-cutoff)<-0.5*dz0)
//		{
//		model[i]=1.0/cutoff;
//		}
//		else
//		if((fabs(z[i]-cutoff))<0.001)
//		{
//		model[i]=0.5/cutoff;
//		}
//		else
//		{
//		model[i]=0.001;
//		}
//

		if(model[i]<0.001) 
		{
			model[i]=0.001;
		}
		fprintf(ModelFile,"%lf, %lf, %lf\n",z[i],model[i],dz[i]);
		fprintf(InitFile,"%lf, %lf\n ",z[i],model[i]);
	}

return 0;
} 

double rang(void)
{
// we assume that the random number seed is set (srand)
//this code generates a gaussianly distributed random number with mean zero
//and a standard deviation of one.  http://www.taygeta.com/random/gaussian.html
         double x1, x2, w, y1, y2;
 
         do {
                 x1 = 2.0 * ((double)rand()/(double)RAND_MAX )- 1.0;
                 x2 = 2.0 * ((double)rand()/(double)RAND_MAX ) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;
	 return y1;
}
