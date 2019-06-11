#include <math.h>
#define MAX(x,y) ( ( (x) > (y) ) ? (x) : (y) )
/*********************************************
 *
 *	FUNCTION : pyrGetSz
 *
 * Calcules the image sizes for an nlevel pyramid
 * and return the total number of images of the structure
 *
 *********************************************/
int pyr_getsize(int nx, int ny, int nz, int nlevels, long nxa[], long nya[], long nza[]);
int pyr_getsize(int nx, int ny, int nz, int nlevels, long nxa[], long nya[], long nza[])
{ int reduc,j;

	  reduc=1;
	  if (nlevels>0) reduc=pow(2.,(double)(nlevels-1));
	  nxa[0]=MAX(1,(nx/reduc)*reduc);
	  nya[0]=MAX(1,(ny/reduc)*reduc);
	  nza[0]=MAX(1,(nz/reduc)*reduc);
	  for (j=1;j<nlevels;j++) {
		 nxa[j]=MAX(1,nxa[j-1]/2);
		 nya[j]=MAX(1,nya[j-1]/2);
		 nza[j]=MAX(1,nza[j-1]/2);
	  }
	  return (nlevels);
}
	  
