#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include "arraymalloc.h"
#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double error, bnorm;
  double tolerance=0.0; //tolerance for convergence. <=0 means do not check

  //main arrays
  double **psi, **zet;
  //temporary versions of main arrays
  double **psitmp, **zettmp;

  //command line arguments
  int scalefactor, numiter;

  double re; // Reynold's number - must be less than 3.7

  //simulation sizes
  int bbase=10;
  int hbase=15;
  int wbase=5;
  int mbase=32;
  int nbase=32;

  int irrotational = 1, checkerr = 0;

  int m,n,b,h,w;
  int iter;
  int i,j;

  // OpenMP variables

  int nthread;

  double tstart, tstop, ttot, titer;

  //do we stop because of tolerance?
  if (tolerance > 0) {checkerr=1;}

  //check command line parameters and parse them

  if (argc <3|| argc >4)
    {
      printf("Usage: cfd <scale> <numiter> [reynolds]\n");
      return 0;
    }

  scalefactor=atoi(argv[1]);
  numiter=atoi(argv[2]);

  if (argc == 4)
    {
      re=atof(argv[3]);
      irrotational=0;
    }
  else
    {
      re=-1.0;
    }

  if(!checkerr)
    {
      printf("Scale Factor = %i, iterations = %i\n",scalefactor, numiter);
    }
  else
    {
      printf("Scale Factor = %i, iterations = %i, tolerance= %g\n",scalefactor,numiter,tolerance);
    }

  if (irrotational)
    {
      printf("Irrotational flow\n");
    }
  else
    {
      printf("Reynolds number = %f\n",re);
    }

  //Calculate b, h & w and m & n
  b = bbase*scalefactor;
  h = hbase*scalefactor;
  w = wbase*scalefactor;
  m = mbase*scalefactor;
  n = nbase*scalefactor;

  re = re / (double)scalefactor;

  nthread = omp_get_max_threads();
  
  printf("Running CFD on %d x %d grid using %d thread(s)\n",m,n,nthread);

  //allocate arrays

  psi    = (double **) arraymalloc2d(m+2,n+2,sizeof(double));
  psitmp = (double **) arraymalloc2d(m+2,n+2,sizeof(double));

  //zero the psi array
  for (i=0;i<m+2;i++)
    {
      for(j=0;j<n+2;j++)
	{
	  psi[i][j]=0.0;
	}
    }

  if (!irrotational)
    {
      //allocate arrays

      zet =   (double **) arraymalloc2d(m+2,n+2,sizeof(double));
      zettmp =(double **) arraymalloc2d(m+2,n+2,sizeof(double));

      //zero the zeta array

      for (i=0;i<m+2;i++)
	{
	  for(j=0;j<n+2;j++)
	    {
	      zet[i][j]=0.0;
	    }
	}
    }
  
  //set the psi boundary conditions

  boundarypsi(psi,m,n,b,h,w);

  //compute normalisation factor for error

  bnorm=0.0;

  for (i=0;i<m+2;i++)
    {
      for (j=0;j<n+2;j++)
	{
	  bnorm += psi[i][j]*psi[i][j];
	}
    }

  if (!irrotational)
    {
      //update zeta BCs that depend on psi
      boundaryzet(zet,psi,m,n);

      //update normalisation

      for (i=0;i<m+2;i++)
	{
	  for (j=0;j<n+2;j++)
	    {
	      bnorm += zet[i][j]*zet[i][j];
	    }
	}
    }

  bnorm=sqrt(bnorm);

  //begin iterative Jacobi loop

  printf("\nStarting main loop...\n\n");
  
  tstart=gettime();

  for(iter=1;iter<=numiter;iter++)
    {
      //calculate psi for next iteration

      if (irrotational)
	{
	  jacobistep(psitmp,psi,m,n);
	}
      else
	{
	  jacobistepvort(zettmp,psitmp,zet,psi,m,n,re);
	}

      //calculate current error if required

      if (checkerr || iter == numiter)
	{
	  error = deltasq(psitmp,psi,m,n);

	  if(!irrotational)
	    {
	      error += deltasq(zettmp,zet,m,n);
	    }

	  error=sqrt(error);
	  error=error/bnorm;
	}

      //copy back

#pragma omp parallel for default(none) private(i,j) shared(psi,psitmp,m,n)
      for(i=1;i<=m;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      psi[i][j]=psitmp[i][j];
	    }
	}

      if (!irrotational)
	{
#pragma omp parallel for default(none) private(i,j) shared(zet,zettmp,m,n)
	  for(i=1;i<=m;i++)
	    {
	      for(j=1;j<=n;j++)
		{
		  zet[i][j]=zettmp[i][j];
		}
	    }
	}

      if (!irrotational)
	{
	  //update zeta BCs that depend on psi
	  boundaryzet(zet,psi,m,n);
	}

      //quit early if we have reached required tolerance

      if (checkerr)
	{
	  if (error < tolerance)
	    {
	      printf("Converged on iteration %d\n",iter);
	      break;
	    }
	}

      //print loop information

      if(iter%printfreq == 0)
	{
	  if (!checkerr)
	    {
	      printf("Completed iteration %d\n",iter);
	    }
	  else
	    {
	      printf("Completed iteration %d, error = %g\n",iter,error);
	    }
	}
    }

  if (iter > numiter) iter=numiter;

  tstop=gettime();

  ttot=tstop-tstart;
  titer=ttot/(double)iter;


  //print out some stats

  printf("\n... finished\n");
  printf("After %d iterations, the error is %g\n",iter,error);
  printf("Time for %d iterations was %g seconds\n",iter,ttot);
  printf("Each iteration took %g seconds\n",titer);

  //output results

  writedatafiles(psi,m,n, scalefactor);

  writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  free(psi);
  free(psitmp);

  if (!irrotational)
    {
      free(zet);
      free(zettmp);
    }

  printf("... finished\n");

  return 0;
}
