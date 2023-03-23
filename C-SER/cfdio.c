#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cfdio.h"
#include "arraymalloc.h"

void writedatafiles(double **psi, int m, int n, int scale)
{
  typedef double Vecvel[2];
  typedef int    Vecrgb[3];

  Vecvel **vel;
  Vecrgb **rgb;

  FILE *cfile, *vfile;

  double modvsq, hue;
  int i,j, ix, iy;
  int nvel, nrgb;

  printf("\n\nWriting data files ...\n");

  vel = (Vecvel **) arraymalloc2d(m,n,sizeof(Vecvel));
  rgb = (Vecrgb **) arraymalloc2d(m,n,sizeof(Vecrgb));

  //calculate velocities and hues

  double v1, v2;

  for (i=0;i<m;i++)
    {
      for (j=0;j<n;j++)
	{
	  vel[i][j][0] =  (psi[i+1][j+2]-psi[i+1][j])/2.0;
	  vel[i][j][1] = -(psi[i+2][j+1]-psi[i][j+1])/2.0;

	  v1 = vel[i][j][0];
	  v2=  vel[i][j][1];

	  modvsq = v1*v1 + v2*v2;

	  hue = pow(modvsq,0.4);

	  hue2rgb(hue,&(rgb[i][j][0]),&(rgb[i][j][1]),&(rgb[i][j][2]));
	}
    }

  //write data

  cfile=fopen("colourmap.dat","w");
  vfile=fopen("velocity.dat","w");

  for (i=0;i<m;i++)
    {
      ix = i+1;

      for (j=0;j<n;j++)
	{
	  iy = j+1;

	  fprintf(cfile,"%i %i %i %i %i\n", ix, iy,
		  rgb[i][j][0], rgb[i][j][1],rgb[i][j][2]);

	  if ((ix-1)%scale == (scale-1)/2 &&
	      (iy-1)%scale == (scale-1)/2    )
	    {
	      fprintf(vfile,"%i %i %f %f\n",
		      ix,iy,vel[i][j][0],vel[i][j][1]);
	    }
	}
    }

  fclose(vfile);
  fclose(cfile);

  free(rgb);
  free(vel);

  printf("... done!\n");
}

void writeplotfile(int m, int n, int scale)
{
  FILE *gnuplot;

  gnuplot = fopen("cfd.plt","w");

  fprintf(gnuplot,"set size square\n");
  fprintf(gnuplot,"set key off\n");
  fprintf(gnuplot,"unset xtics\n");
  fprintf(gnuplot,"unset ytics\n");

  fprintf(gnuplot,"set xrange [%i:%i]\n",1-scale,m+scale);
  fprintf(gnuplot,"set yrange [%i:%i]\n",1-scale,n+scale);

  fprintf(gnuplot,"plot \"colourmap.dat\" w rgbimage, \"velocity.dat\" u 1:2:(%d*0.75*$3/sqrt($3**2+$4**2)):(%d*0.75*$4/sqrt($3**2+$4**2)) with vectors  lc rgb \"#7F7F7F\"",scale,scale);

  fclose(gnuplot);

  printf("\nWritten gnuplot script 'cfd.plt'\n");
}


void hue2rgb(double hue, int *r, int *g, int *b)
{
  int rgbmax = 255;

  *r = (int)(rgbmax*colfunc(hue-1.0));
  *g = (int)(rgbmax*colfunc(hue-0.5));
  *b = (int)(rgbmax*colfunc(hue    ));
}


double colfunc(double x)
{
  double absx;

  double x1=0.2;
  double x2=0.5;

  absx=fabs(x);

  if (absx > x2)
    {
      return 0.0;
    }
  else if (absx < x1)
    {
      return 1.0;
    }
  else
    {
      return 1.0-pow((absx-x1)/(x2-x1),2);
    }
}


#include <sys/time.h>

/* wall-clock time */

double gettime(void)
{
  struct timeval tp;
  gettimeofday (&tp, NULL);
  return tp.tv_sec + tp.tv_usec/(double)1.0e6;
}
