/*******************************************************************************

        Transpose interpolation of grid locations of points using filtering 
	technique (kernel) with zero boundary condition

        Created by Se Young Chun, Nov 2, 2006, the University of Michigan

*******************************************************************************/

extern void     FiltInterpolatedValue2DTranZero (
                        double  *Bout,    /* output interpolated value */
                        double  *Bcoeff,  /* input B-spline coeff array */
                        int     Width,    /* width of coeff array */
                        int     Height,   /* height of coeff array */
                        int     nx,       /* width of Bout */
                        int     offx,     /* offset of x direction */
                        int     mx,       /* magnif facor for x direction */
                        int     ny,       /* height of Bout */
                        int     offy,     /* offset of y direction */
                        int     my,       /* magnif facor for y direction */
                        double  *kernelx, /* kernelx applied in x direction */
                        int     nkx,      /* kernelx length */
                        double  *kernely, /* kernely applied in y direction */
                        int     nky       /* kernely length */
                );
