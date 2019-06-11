/*******************************************************************************

        These routines are to calculate the transpose of interpolation matrices

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

extern double	InterpolatedValue3DTranMirr (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradX3DTranMirr (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradY3DTranMirr (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradZ3DTranMirr (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedValue3DTranZero (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradX3DTranZero (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradY3DTranZero (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);

extern double	InterpolatedGradZ3DTranZero (
                        float   *Bout,          /* output + prev values */
                        long    Width,          /* width of Bout */
                        long    Height,         /* height of Bout */
			long	Slice,		/* slice of the image */
                        float   Bcoeff,         /* input B-spline coefficient */
                        double  x,              /* x coord to interpolate */
                        double  y,              /* y coord to interpolate */
			double	z,		/* z coord to interpolate */
                        long    SplineDegree    /* degree of spline model */
		);
