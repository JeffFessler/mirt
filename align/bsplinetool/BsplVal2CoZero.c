/*******************************************************************************
      B-spline coefficient computatation for given signal, especially for
      3D image

      Date: June 13, 2003
      Programmer: Jeongtae Kim
                  (EECS:Systems, U of Michigan)
      Note: 1. Parts of the program for image interpolation are originated
         from M. Unser's group.

      Bug report: jeongtae@umich.edu, +1-734-647-8390

      Extended to zero end condition case, 
      Se Young Chun, May 21, 2007, the University of Michigan

*******************************************************************************/

#include	<float.h>
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include        "mex.h"
#include        "matrix.h"

static void     ConvertToInterpolationCoefficientsZero (
                        double  c[],            /* samples --> coefficients */
                        long    DataLength,     /* number of samples/coeffs */
                        double  z[],            /* poles */
                        long    NbPoles,        /* number of poles */
                        double  Tolerance       /* admissible relative error */
                );

static void     GetColumn (
                        float   *Image,         /* input image array */
                        long    x,              /* x coord of selected line */
                        long    z,              /* z coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height          /* width of the image */
                );

static void     GetRow (
                        float   *Image,         /* input image array */
                        long    y,              /* y coord of selected line */
                        long    z,              /* z coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height
                );

static void     GetSlice (
                        float   *Image,         /* input image array */
                        long    x,              /* x coord of selected line */
                        long    y,              /* y coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height,
                        long    Slice
                );

#if 0
static double   InitialCausalCoefficient (
                        double  c[],            /* coefficients */
                        long    DataLength,     /* number of coefficients */
                        double  z,              /* actual pole */
                        double  Tolerance       /* admissible relative error */
                                );

static double   InitialAntiCausalCoefficient (
                        double  c[],            /* coefficients */
                        long    DataLength,     /* number of samples/coeffs */
                        double  z                       /* actual pole */
                                );
#endif

static void     PutColumn (
                        float   *Image,         /* input image array */
                        long    x,              /* x coord of selected line */
                        long    z,              /* z coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height          /* width of the image */
                );

static void     PutRow (
                        float   *Image,         /* input image array */
                        long    y,              /* y coord of selected line */
                        long    z,              /* z coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height
                );

static void     PutSlice (
                        float   *Image,         /* input image array */
                        long    x,              /* x coord of selected line */
                        long    y,              /* y coord of selected line */
                        double  Line[],         /* output linear array */
                        long    Width,          /* length of the line */
                        long    Height,
                        long    Slice
                );




static void	ConvertToInterpolationCoefficientsZero (
			double	c[],
			long	DataLength,
			double	z[],	
			long	NbPoles,
			double	Tolerance
				)
{ /* begin ConvertToInterpolationCoefficientsZero */

	double	*ct, Lambda = 1.0;
	long	n, k, Horizon;
	
	/* special case required by mirror boundaries */
	if (DataLength == 1L) {
		return;
	}
	/* compute the overall gain */
	for (k = 0L; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	}
	/* apply the gain */
	for (n = 0L; n < DataLength; n++) {
		c[n] *= Lambda;
	}
	/* loop over all poles */
	for (k = 0L; k < NbPoles; k++) {
		/* causal initialization; none for zero end condition, c[0] */
		/* causal recursion */
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		/* anticausal initialization; no for zero end condition */

        	if (Tolerance > 0.0) {
                	Horizon = (long)ceil( (log(Tolerance)-log(
				c[DataLength-1L])) / log(fabs(z[k])));
        	}
		else
		{
			Horizon = (long)ceil((-15-log(c[DataLength-1L])) 
				/ log(fabs(z[k])));
		}
		if (Horizon > 1) {
			ct = (double*)malloc((size_t)Horizon 
					* (long)sizeof(double));
			ct[0] = z[k]*c[DataLength - 1L];
                	for (n = 1L; n < Horizon; n++) {
				ct[n] = z[k]*ct[n - 1L];
                	}
			ct[Horizon - 1L] = -z[k]*ct[Horizon - 1L];
			for (n = Horizon - 2L; 0 <= n; n--) {
				ct[n] = z[k] * (ct[n + 1L] - ct[n]);
			}

			c[DataLength - 1L] = z[k] * (ct[0] - c[DataLength - 1L]);
			free(ct);
		}
		else {
			c[DataLength - 1L] = -z[k] * c[DataLength - 1L];			
		}
		/* anticausal recursion */
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
} /* end ConvertToInterpolationCoefficientsZero */

static void	GetColumn (
			float	*Image,
			long	x,
			long    z,    
			double	Line[],
			long	Width,
			long	Height	
		)
{ /* begin GetColumn */
	long	y;
	Image = Image + (ptrdiff_t)(x + z * Width * Height);
	for (y = 0L; y < Height; y++) {
		Line[y] = (double)*Image;
		Image += (ptrdiff_t)Width;
	}
} /* end GetColumn */

static void	GetRow (
			float	*Image,
			long	y,
			long    z,      
			double	Line[],	
			long	Width,	
			long    Height
		)
{ /* begin GetRow */
	long	x;
	Image = Image + (ptrdiff_t)(y * Width + z * Width * Height);
	for (x = 0L; x < Width; x++) {
		Line[x] = (double)*Image++;
	}
} /* end GetRow */

static void	GetSlice (
			float	*Image,
			long    x,     
			long	y,
			double	Line[],	
			long	Width,
			long    Height,
			long    Slice
		)
{ /* begin GetSlice */
	long	z;
	Image = Image + (ptrdiff_t)(x + y * Width);
	for (z = 0L; z < Slice; z++) {
		Line[z] = (double)*Image;
		Image += (ptrdiff_t)(Width * Height);
	}
} /* end GetSlice */

static void	PutColumn (
			float	*Image,	
			long	x,	
			long    z,      
			double	Line[],	
			long	Width,	
			long	Height	
		)
{ /* begin PutColumn */
	long	y;
	Image = Image + (ptrdiff_t)(x + z * Width * Height);
	for (y = 0L; y < Height; y++) {
		*Image = (float)Line[y];
		Image += (ptrdiff_t)Width;
	}
} /* end PutColumn */

static void	PutRow (
			float	*Image,	
			long	y,
			long    z,  
			double	Line[],	
			long	Width,	 
                        long    Height
		)
{ /* begin PutRow */
	long	x;
	Image = Image + (ptrdiff_t)(y * Width + z * Width * Height);
	for (x = 0L; x < Width; x++) {
		*Image++ = (float)Line[x];
	}
} /* end PutRow */

static void	PutSlice (
			float	*Image,	
			long    x,      
			long	y,	
			double	Line[],	
			long	Width,	
			long    Height,
			long    Slice
		)
{ /* begin PutSlice */
	long	z;
	Image = Image + (ptrdiff_t)(x + y * Width);
	for (z = 0L; z < Slice; z++) {
		*Image = (float)Line[z] ;
		Image += (ptrdiff_t)(Width * Height);
	}
} /* end PutSlice */


extern int	SamplesToCoefficientsZero (
			float	*Image,
			long	Width,
			long	Height,	
			long    Slice,     
			long	SplineDegree
		)
{ /* begin SamplesToCoefficientsZero */
	double	*Line;
	double	Pole[2];
	long	NbPoles;
	long	x, y, z;
	/* recover the poles from a lookup table */
	switch (SplineDegree) {
		case 2L:
			NbPoles = 1L;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1L;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2L;
			Pole[0] = sqrt(664.0 - sqrt(438976.0))+sqrt(304.0)-19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0))-sqrt(304.0)-19.0;
			break;
		case 5L:
			NbPoles = 2L;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) 
					+ sqrt(105.0 / 4.0) - 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) 
					- sqrt(105.0 / 4.0) - 13.0 / 2.0;
			break;
		default:
			printf("Invalid spline degree\n");
			return(1);
	}
	/* convert the image samples into interpolation coefficients */
	
	/* in-place separable process, along x */

	Line = (double *)malloc((size_t)(Width * (long)sizeof(double)));
	
	if (Line == (double *)NULL) {
		printf("Row allocation failed\n");
		return(1);
	}
	for (z = 0L; z < Slice; z++){
	  for (y = 0L; y < Height; y++) {
	    GetRow(Image, y, z, Line, Width, Height);
	    ConvertToInterpolationCoefficientsZero(Line, Width, Pole, NbPoles, 
		DBL_EPSILON);
	    PutRow(Image, y, z, Line, Width, Height);
	  }
	}
	free(Line);

	/* in-place separable process, along y */
	Line = (double *)malloc((size_t)(Height * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Column allocation failed\n");
		return(1);
	}
        for (z = 0L; z < Slice; z++){
	  for (x = 0L; x < Width; x++) {
		GetColumn(Image, x, z, Line, Width, Height);
		ConvertToInterpolationCoefficientsZero(Line, Height, Pole, 
			NbPoles, DBL_EPSILON);
		PutColumn(Image, x, z, Line, Width, Height);
	  }
	}
	free(Line);

	/* in-place separable process, along z */
	if (Slice >= 2){
	  Line = (double *)malloc((size_t)(Slice * (long)sizeof(double)));
	  if (Line == (double *)NULL) {
		printf("Slice allocation failed\n");
		return(1);
	  }
	  for (y = 0L; y < Height; y++){
	    for (x = 0L; x < Width; x++) {
	      GetSlice(Image, x, y, Line, Width, Height, Slice);
	      ConvertToInterpolationCoefficientsZero(Line, Slice, Pole, NbPoles,
 			DBL_EPSILON); 
	      PutSlice(Image, x, y, Line, Width, Height, Slice);
	    }
	  }
	  free(Line); 
	  } 
	return(0);
} /* end SamplesToCoefficientsZero */


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])   
{
      	const    float           *im_in;
      	float                    *im_out;
      	long     width, height, slice, splinedegree;

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nlhs != 1)
	{
		mexPrintf("need nlhs=1\n");
		return;
	}

        if (nrhs == 1)
        {
                splinedegree = 3;
        }
        else if (nrhs == 2)
        {
                splinedegree=mxGetScalar(prhs[1]);
        }
	else
	{
              mexPrintf("COEFF = BsplVal2CoZero(Val, deg);\n");
              mexPrintf("[in]\n");
              mexPrintf("\tVal    : original values - single type\n");
              mexPrintf("\tdeg    : (opt) degree of spline basis {default:3}\n");
              mexPrintf("[out]\n");
              mexPrintf("\tCOEFF  : bspline coefficients array\n");
              return;
      	}

	if (mxIsSingle(prhs[0]) != 1)
        {
                mexErrMsgTxt("First argument must be of single type\n");
                return;
        }

       	if (mxGetNumberOfDimensions(prhs[0]) == 2)
        {
      		width  = mxGetDimensions(prhs[0])[0];
      		height = mxGetDimensions(prhs[0])[1];
      		slice  = 1;
        }
        else if (mxGetNumberOfDimensions(prhs[0]) == 3)
        {
      		width  = mxGetDimensions(prhs[0])[0];
      		height = mxGetDimensions(prhs[0])[1];
      		slice  = mxGetDimensions(prhs[0])[2];
        }
        else
        {
                mexErrMsgTxt("First argument must be 2D or 3D\n");
                return;
        }

      	im_in = mxGetData(prhs[0]);
     
      	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), 
			mxGetDimensions(prhs[0]), mxSINGLE_CLASS, mxREAL);      
      	im_out = mxGetData(plhs[0]);

     	memcpy(im_out, im_in, width*height*slice*sizeof(float));
 
      	if (SamplesToCoefficientsZero(im_out,width, height, slice, splinedegree))
	{
		mxFree(plhs[0]);
                mexErrMsgTxt("Change of basis failed");
                return;
	}
}
