/*******************************************************************************

        Interpolation of grid locations of points using filtering technique
        (kernel) with zero boundary condition

        Created by Se Young Chun, Oct 20, 2006, the University of Michigan

*******************************************************************************/

extern void     FiltInterpolatedValue3DZero (
                        double  *Bout,     /* output interpolated value */
                        double  *Bcoeff,   /* input B-spline coeff array */
                        int     Width,     /* width of coeff array */
                        int     Height,    /* height of coeff array */
                        int     Slice,     /* slice of coeff array */
                        int	nx,        /* width of Bout */
                        int	offx,      /* offset of x direction */
                        int     mx,        /* magnif factor for x direction */
                        int     ny,        /* height of Bout */
                        int     offy,      /* offset of y direction */
                        int     my,        /* magnif factor for y direction */
                        int     nz,        /* slice of Bout */
                        int     offz,      /* offset of z direction */
                        int     mz,        /* magnif factor for z direction */
			double  *kernelx,  /* kernel to be applied, 1D */
			int	nkx,	   /* kernel length */
			double  *kernely,  /* kernel to be applied, 1D */
			int	nky,	   /* kernel length */
			double  *kernelz,  /* kernel to be applied, 1D */
			int	nkz	   /* kernel length */
                );
