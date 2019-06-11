/*******************************************************************************
      B-spline coefficient computatation for given signal, especially for 
      3D image
  
      Date: June 13, 2003
      Programmer: Jeongtae Kim
		  (EECS:Systems, U of Michigan)
      Note: 1. Parts of the program for image interpolation are originated 
	 from M. Unser's group.
     
      Bug report: jeongtae@umich.edu, +1-734-647-8390


      MEX interface added by Se Young Chun, Oct 22, 2006, 
      the University of Michigan
  
*******************************************************************************/

#include	<float.h>
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	"mex.h"
#include	"matrix.h"

#define		 BASIS_SHIFT \
				((double)(0.5 - sqrt(1.0 / 12.0)))

static void	ConvertToInterpolationCoefficientsMirr (
			double	c[],		/* samples --> coefficients */
			long	DataLength,	/* number of samples/coeffs */
			double	z[],		/* poles */
			long	NbPoles,	/* number of poles */
			double	Tolerance	/* admissible relative error */
		);

static void	GetColumn (
			float	*Image,		/* input image array */
			long	x,		/* x coord of selected line */
			long    z,	      /* z coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long	Height		/* width of the image */
		);

static void	GetRow (
		 	float	*Image,		/* input image array */
			long	y,		/* y coord of selected line */
			long    z,	      /* z coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long    Height
		);

static void	GetSlice (
			float	*Image,		/* input image array */
			long    x,	      /* x coord of selected line */
			long	y,		/* y coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long    Height,
			long    Slice
		);

static double	InitialCausalCoefficient (
			double	c[],		/* coefficients */
			long	DataLength,	/* number of coefficients */
			double	z,		/* actual pole */
			double	Tolerance	/* admissible relative error */
				);

static double	InitialAntiCausalCoefficient (
			double	c[],		/* coefficients */
			long	DataLength,	/* number of samples/coeffs */
			double	z		/* actual pole */
				);

static void	PutColumn (
			float	*Image,		/* input image array */
			long	x,		/* x coord of selected line */
			long    z,	      /* z coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long	Height		/* width of the image */
		);

static void	PutRow (
			float	*Image,		/* input image array */
			long	y,		/* y coord of selected line */
			long    z,	      /* z coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long    Height
		);

static void	PutSlice (
			float	*Image,		/* input image array */
			long    x,	      /* x coord of selected line */
			long	y,		/* y coord of selected line */
			double	Line[],		/* output linear array */
			long	Width,		/* length of the line */
			long    Height,
			long    Slice
		);

static int      ShiftedSamplesToCoefficients (
			float   *Image,	 /* in-place processing */
			long    Width,	  /* width of the image */
			long    Height,	 /* height of the image */
			long    Slice	   /* height of the image */
		);

/* 
   This part came from:
   T. Blu, P. Thevenaz, M. Unser, "Linear Interpolation Revitalized,"
   IEEE Transactions on Image Processing, vol. 13, no. 5, pp. 710-719, May 2004
*/
static void     ShiftedBasisCoefficients1D (
			double  c[],	    /* in-place proc; 1D array */
			long    DataLength      /* length of the 1D array */
		);

static double   ShiftedBasisFirstCoefficient1D (
			double  c[]	     /* input 1D array */
		);

static void	ConvertToInterpolationCoefficientsMirr (
			double	c[],	
			long	DataLength,
			double	z[],	
			long	NbPoles,
			double	Tolerance
		)
{ /* begin ConvertToInterpolationCoefficientsMirr */

	double	Lambda = 1.0;
	long	n, k;

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
		/* causal initialization */
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
		/* causal recursion */
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		/* anticausal initialization */
		c[DataLength - 1L] = InitialAntiCausalCoefficient(
			c, DataLength, z[k]);
		/* anticausal recursion */
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}

} /* end ConvertToInterpolationCoefficientsMirr */

static double	InitialCausalCoefficient (
			double	c[],	
			long	DataLength,
			double	z,		
			double	Tolerance
		)
{ /* begin InitialCausalCoefficient */

	double	Sum, zn, z2n, iz;
	long	n, Horizon;

	/* this initialization corresponds to mirror boundaries */
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		/* accelerated loop */
		zn = z;
		Sum = c[0];
		for (n = 1L; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(Sum);
	}
	else {
		/* full loop */
		zn = z;
		iz = 1.0 / z;
		z2n = pow(z, (double)(DataLength - 1L));
		Sum = c[0] + z2n * c[DataLength - 1L];
		z2n *= z2n * iz;
		for (n = 1L; n <= DataLength - 2L; n++) {
			Sum += (zn + z2n) * c[n];
			zn *= z;
			z2n *= iz;
		}
		return(Sum / (1.0 - zn * zn));
	}

}
 /* end InitialCausalCoefficient */

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

static double	InitialAntiCausalCoefficient (
			double	c[],
			long	DataLength,	
			double	z	
		)
{ /* begin InitialAntiCausalCoefficient */

	/* this initialization corresponds to mirror boundaries */
	return((z/(z*z-1.0))*(z*c[DataLength-2L]+c[DataLength-1L]));

} /* end InitialAntiCausalCoefficient */

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

static int	SamplesToCoefficientsMirr (
			float	*Image,	
			long	Width,
			long	Height,
			long    Slice, 
			long	SplineDegree
		)
{ /* begin SamplesToCoefficientsMirr */

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
			Pole[0] = sqrt(664.0-sqrt(438976.0))+sqrt(304.0)-19.0;
			Pole[1] = sqrt(664.0+sqrt(438976.0))-sqrt(304.0)-19.0;
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
	    ConvertToInterpolationCoefficientsMirr(Line, Width, Pole, 
		NbPoles, DBL_EPSILON);
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
		ConvertToInterpolationCoefficientsMirr(Line, Height, Pole, 
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
	      ConvertToInterpolationCoefficientsMirr(Line, Slice, Pole, 
			NbPoles, DBL_EPSILON); 
	      PutSlice(Image, x, y, Line, Width, Height, Slice);
	    }
	  }
	  free(Line); 
	} 

	return(0);

} /* end SamplesToCoefficientsMirr */

static void  	ShiftedBasisCoefficients1D (
			double  c[],
			long    DataLength 
		)

{ /* begin ShiftedBasisCoefficients1D */

	double  Factor = 1.0 / (1.0 - BASIS_SHIFT);
	double  Pole = BASIS_SHIFT / (BASIS_SHIFT - 1.0);
	long    n;

	c[0] = ShiftedBasisFirstCoefficient1D(c);
	for (n = 1L; (n < DataLength); n++) {
		c[n] = Factor * c[n] + Pole * c[n - 1L];
	}
} /* end ShiftedBasisCoefficients1D */

static double   ShiftedBasisFirstCoefficient1D (
			double  c[] 
		)
{ /* begin ShiftedBasisFirstCoefficient1D */
	return(c[0]);
} /* end ShiftedBasisFirstCoefficient1D */

static int      ShiftedSamplesToCoefficients (
			float   *Image,	 /* in-place processing */
			long    Width,	  /* width of the image */
			long    Height,	 /* height of the image */
			long    Slice	  /* height of the image */
				)

{ /* begin ShiftedSamplesToCoefficients */

	double  *Line;
	long    x, y, z;

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
	       	ShiftedBasisCoefficients1D(Line, Width);
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
		ShiftedBasisCoefficients1D(Line, Height);
		PutColumn(Image, x, z, Line, Width, Height);
	  }
	}
	free(Line);

	/* in-place separable process, along z */
	if (Slice >= 2) {
	  Line = (double *)malloc((size_t)(Slice * (long)sizeof(double)));
	  if (Line == (double *)NULL) {
		printf("Slice allocation failed\n");
		return(1);
	  }
	  for (y = 0L; y < Height; y++){
	    for (x = 0L; x < Width; x++) {
	      GetSlice(Image, x, y, Line, Width, Height, Slice);
	      ShiftedBasisCoefficients1D(Line, Slice);
	      PutSlice(Image, x, y, Line, Width, Height, Slice);
	    }
	  }
	  free(Line); 
	} 

	return(0);
} /* end ShiftedSamplesToCoefficients */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])   
{
      	const float *im_in;
      	float *im_out;
      	long width, height, slice, splinedeg;

	if (nlhs == 0 && nrhs == 1 && mxIsChar(prhs[0])) // check
		return;

	if (nlhs != 1)
	{
		mexPrintf("need nlhs=1\n");
		return;
	}

	if (nrhs == 1)
	{
		splinedeg = 3;
	}
	else if (nrhs == 2)
	{
		splinedeg=mxGetScalar(prhs[1]);
	}
	else
	{
	     mexPrintf("COEFF = BsplVal2CoMirr(Val, deg);\n");
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

	if (splinedeg == 1)
	{
      	  if (ShiftedSamplesToCoefficients(im_out,width,height,slice))
	  {
		mxFree(plhs[0]);
		mexErrMsgTxt("Change of basis failed");
		return;
	  }
	}
	else 
	{ 
      	  if (SamplesToCoefficientsMirr(im_out,width,height,slice,splinedeg))
	  {
		mxFree(plhs[0]);
		mexErrMsgTxt("Change of basis failed");
		return;
	  }
	}
}
