/*******************************************************************************

	These routines are to calculate the transpose of interpolation
	matrices.

	Created by Se Young Chun, Jan 19, 2007, the University of Michigan
	
*******************************************************************************/

extern double	InterpolatedValueTranMirr (
			float	*Bout,		/* output + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradXTranMirr (
			float	*Bout,		/* output  + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradYTranMirr (
			float	*Bout,		/* output  + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedValueTranZero (
			float	*Bout,		/* output + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradXTranZero (
			float	*Bout,		/* output  + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);

extern double	InterpolatedGradYTranZero (
			float	*Bout,		/* output  + prev values */
			long	Width,		/* width of Bout */
			long	Height,		/* height of Bout */
			float	Bcoeff,		/* input B-spline coefficient */
			double	x,		/* x coord to interpolate */
			double	y,		/* y coord to interpolate */
			long	SplineDegree	/* degree of spline model */
		);
