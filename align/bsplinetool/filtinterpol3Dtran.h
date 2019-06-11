/*******************************************************************************

        Transpose interpolation of grid locations of points using filtering 
	technique (kernel) with zero boundary condition

        Created by Se Young Chun, Nov 2, 2006, the University of Michigan

*******************************************************************************/

extern void     FiltInterpolatedValue3DTranZero (
                        double  *Bout,    /* output interpolated value */
                        double  *Bcoeff,  /* input B-spline coeff array */
                        int     Width,    /* width of coeff array */
                        int     Height,   /* height of coeff array */
                        int     Slice,    /* slice of coeff array */
                        int     nx,       /* width of Bout */
                        int     offx,     /* offset of x direction */
                        int     mx,       /* magnif factor for x direction */
                        int     ny,       /* height of Bout */
                        int     offy,     /* offset of y direction */
                        int     my,       /* magnif factor for y direction */
                        int     nz,       /* slice of Bout */
                        int     offz,     /* offset of z direction */
                        int     mz,       /* magnif factor for z direction */
                        double  *kernelx, /* kernelx applied in x direction */
                        int     nkx,      /* kernelx length */
                        double  *kernely, /* kernely applied in y direction */
                        int     nky,      /* kernely length */
                        double  *kernelz, /* kernelz applied in z direction */
                        int     nkz       /* kernelz length */
                );
