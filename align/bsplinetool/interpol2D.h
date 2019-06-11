/*******************************************************************************

	Original source code comes from http://bigwww.epfl.ch/,
	interpolation package written by Philippe Thevenaz, January 3, 2006
	based on the following paper:
        P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
        IEEE Transactions on Medical Imaging,
        vol. 19, no. 7, pp. 739-758, July 2000.

	This is the trivial extension of original code to add gradient value 
	interpolations for mirror and zero boundary condition. 

	Modified by Se Young Chun, Oct 1, 2006, the University of Michigan

*******************************************************************************/

#define                 BASIS_SHIFT \
                                ((double)(0.5 - sqrt(1.0 / 12.0)))


/* point interpolation (2D) using mirror end condition */
extern double	InterpolatedValueMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradXMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradYMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

/* point interpolation (2D) using zero end condition */
extern double	InterpolatedValueZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradXZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradYZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);


