/*******************************************************************************

        Original source code comes from http://bigwww.epfl.ch/,
        interpolation package written by Philippe Thevenaz, January 3, 2006
        based on the following paper:
        P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
        IEEE Transactions on Medical Imaging,
        vol. 19, no. 7, pp. 739-758, July 2000.

        This is the trivial extension of original code to 3D interpolation
	of values, gradients with mirror and zero boundary condition.

        Coded by Se Young Chun, Oct 1, 2006, the University of Michigan

*******************************************************************************/

#define                 BASIS_SHIFT \
                                ((double)(0.5 - sqrt(1.0 / 12.0)))


extern double	InterpolatedValue3DMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradX3DMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradY3DMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradZ3DMirr (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedAll3DMirr (
			double	*gradx,		/* output gradient x */
			double	*grady,		/* output gradient y */
			double	*gradz,		/* output gradient z */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double   InterpolatedThree3DMirr (
                        double  *value1,
                        double  *gradx1,
                        double  *grady1,
                        double  *gradz1,
                        double  *value2,
                        double  *gradx2,
                        double  *grady2,
                        double  *gradz2,
                        double  *value3,
                        double  *gradx3,
                        double  *grady3,
                        double  *gradz3,
                        float   *Bcoeff1,
                        float   *Bcoeff2,
                        float   *Bcoeff3,
                        long    Width,
                        long    Height,
                        long    Slice,
                        double  x,
                        double  y,
                        double  z,
                        long    SplineDegree
                );


extern double	InterpolatedValue3DZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradX3DZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradY3DZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradZ3DZero (
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedAll3DZero (
			double	*gradx,		/* output gradient x */
			double	*grady,		/* output gradient y */
			double	*gradz,		/* output gradient z */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of the image */
			long	Height,		/* height of the image */
			long	Slice,		/* slice of the image */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			double	z,		/* z coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern void   InterpolatedThree3DZero (
                        double  *value1,
                        double  *gradx1,
                        double  *grady1,
                        double  *gradz1,
                        double  *value2,
                        double  *gradx2,
                        double  *grady2,
                        double  *gradz2,
                        double  *value3,
                        double  *gradx3,
                        double  *grady3,
                        double  *gradz3,
                        float   *Bcoeff1,
                        float   *Bcoeff2,
                        float   *Bcoeff3,
                        long    Width,
                        long    Height,
                        long    Slice,
                        double  x,
                        double  y,
                        double  z,
                        long    SplineDegree
                );


