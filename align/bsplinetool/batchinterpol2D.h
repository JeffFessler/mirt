/*******************************************************************************

        Batch interpolation for each location of points (+ deformation)
        using point interpolation functions with mirror and zero end condition

        Created by Se Young Chun, Sep 30, 2006, the University of Michigan

*******************************************************************************/

extern void 	BatchInterpolatedValue2DMirr (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedValue2DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradX2DMirr (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedGradX2DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradY2DMirr (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedGradY2DMirrWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedValue2DZero (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedValue2DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradX2DZero (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedGradX2DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

extern void 	BatchInterpolatedGradY2DZero (
			float	*Bout,		/* output interpolated */
			float	*Bcoeff,	/* input B-spline coeff array */
			long	Width,		/* width of coeff array */
			long	Height,		/* height of coeff array */
			double	nx,		/* width of Bout */
			double	offx,		/* offset of x direction */
			double	mx,		/* magnif factor for x */
			double	ny,		/* height of Bout */
			double	offy,		/* offset of y direction */
			double	my,		/* magnif factor for y */
			long	SplineDegree	/* degree of spline model */
		);

extern void     BatchInterpolatedGradY2DZeroWarp (
                        float   *Bout,          /* output interpolated */
                        float   *Bcoeff,        /* input B-spline coeff array */
                        long    Width,          /* width of coeff array */
                        long    Height,         /* height of coeff array */
                        double  nx,             /* width of Bout */
                        double  offx,           /* offset of x direction */
                        double  mx,             /* magnif factor for x */
                        float   *WarpX,         /* deformed value in x direct */
                        double  ny,             /* height of Bout */
                        double  offy,           /* offset of y direction */
                        double  my,             /* magnif factor for y */
                        float   *WarpY,         /* deformed value in y direct */
                        long    SplineDegree    /* degree of spline model */
                );

