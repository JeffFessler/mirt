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

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"interpol3D.h"
 
extern double	InterpolatedValue3DMirr (
			float	*Bcoeff,
			long	Width,	
			long	Height,
			long	Slice,	
			double	x,	
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedValue3DMirr */

	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	interpolated;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree == 1L) {
	   /* perform the shift */
           x -= BASIS_SHIFT;
           y -= BASIS_SHIFT;
           z -= BASIS_SHIFT;
	
           /* compute the interpolation indexes */
           i = (long)floor(x);
           j = (long)floor(y);
           l = (long)floor(z);
           for (k = 0L; k <= 1; k++) {
                xIndex[k] = i++;
                yIndex[k] = j++;
                zIndex[k] = l++;
           }
	}
	else {
	   if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	   else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 1L:
			/* x */
        		w = x - (double)xIndex[0];
        		xWeight[0] = 1.0 - w;
        		xWeight[1] = w;
       	 		/* y */
        		w = y - (double)yIndex[0];
        		yWeight[0] = 1.0 - w;
        		yWeight[1] = w;
       	 		/* z */
        		w = z - (double)zIndex[0];
        		zWeight[0] = 1.0 - w;
        		zWeight[1] = w;
			break;
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;
			break;
		case 6L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= xWeight[0] / 720.0;
			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			xWeight[6] = 1.0 / 2.0 + w;
			xWeight[6] *= xWeight[6] * xWeight[6];
			xWeight[6] *= xWeight[6] / 720.0;
			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[6];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= yWeight[0] / 720.0;
			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			yWeight[6] = 1.0 / 2.0 + w;
			yWeight[6] *= yWeight[6] * yWeight[6];
			yWeight[6] *= yWeight[6] / 720.0;
			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[6];
			/* z */
			w = z - (double)zIndex[3];
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0] * zWeight[0];
			zWeight[0] *= zWeight[0] / 720.0;
			zWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			zWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			zWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			zWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			zWeight[6] = 1.0 / 2.0 + w;
			zWeight[6] *= zWeight[6] * zWeight[6];
			zWeight[6] *= zWeight[6] / 720.0;
			zWeight[5] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[6];
			break;
		case 7L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			xWeight[2] = (397.0 / 7.0 - w*(245.0/3.0 + w*(-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			xWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w * (-3.0 + w)*(-3.0 + w2)))))/144.0;
			xWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			xWeight[7] = w2;
			xWeight[7] *= xWeight[7] * xWeight[7];
			xWeight[7] *= w / 5040.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			yWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			yWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				*(19.0 + w*(-3.0 + w) * (-3.0 + w2))))) / 144.0;
			yWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			yWeight[7] = w2;
			yWeight[7] *= yWeight[7] * yWeight[7];
			yWeight[7] *= w / 5040.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7];
			/* z */
			w = z - (double)zIndex[3];
			zWeight[0] = 1.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] * zWeight[0];
			zWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			zWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			zWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			zWeight[3] = (2416.0/35.0 + w2*(-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			zWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			zWeight[5] = (40.0/7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			zWeight[7] = w2;
			zWeight[7] *= zWeight[7] * zWeight[7];
			zWeight[7] *= w / 5040.0;
			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[2] - zWeight[3]
				- zWeight[4] - zWeight[5] - zWeight[7];
			break;
		case 8L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] / 40320.0;
			w2 = w * w;
			xWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0/16.0 + w * (-15.0/4.0 + w*(9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			xWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 
				+ w * (-1.0 + w)))))))) / 720.0;
			xWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			xWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			xWeight[7] = (39.0 / 16.0 - w * (-6.0 
				+ w * (-9.0 / 2.0 + w2)))
				*(21.0/16.0 + w*(15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			xWeight[8] = 1.0 / 2.0 + w;
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8] / 40320.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[7] - xWeight[8];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] / 40320.0;
			w2 = w * w;
			yWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0/4.0 + w*(9.0/2.0 + w
				* (-3.0 + w)))) / 5040.0;
			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			yWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			yWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			yWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			yWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(15.0/4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			yWeight[8] = 1.0 / 2.0 + w;
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8] / 40320.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[7] - yWeight[8];
			/* z */
			w = z - (double)zIndex[4];
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] / 40320.0;
			w2 = w * w;
			zWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0/2.0 + w2)))
				* (21.0 / 16.0 + w*(-15.0/4.0 + w*(9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			zWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				*(2275.0/16.0 + w*(-487.0/8.0 + w*(-85.0/8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			zWeight[3] = (310661.0 / 1792.0 - w*(14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(53.0 / 4.0 + w * (-8.0 
				+ w * (-1.0 + w)))))))) / 720.0;
			zWeight[4] = (2337507.0/8960.0 + w2*(-2601.0/16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			zWeight[5] = (310661.0 / 1792.0 - w * (-14219.0/64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			zWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w * (15.0/4.0 + w*(9.0/2.0 + w
				* (3.0 + w)))) / 5040.0;
			zWeight[8] = 1.0 / 2.0 + w;
			zWeight[8] *= zWeight[8];
			zWeight[8] *= zWeight[8];
			zWeight[8] *= zWeight[8] / 40320.0;
			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[5] 
				- zWeight[7] - zWeight[8];
			break;
		case 9L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 
				+ w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[6] - yWeight[7] - yWeight[9];
			/* z */
			w = z - (double)zIndex[4];
			zWeight[0] = 1.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] * (1.0 - w) / 362880.0;
			zWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			zWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			zWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			zWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			zWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			zWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			zWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			zWeight[9] = w2 * w2;
			zWeight[9] *= zWeight[9] * w / 362880.0;
			zWeight[8] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[5] 
				- zWeight[6] - zWeight[7] - zWeight[9];
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	if (SplineDegree == 1L) {
	   for (k = 0L; k <= 1; k++) {
                xIndex[k] = (xIndex[k] < 0L) ? (0L)
                        : ((Width <= xIndex[k]) ? (Width - 1L) : (xIndex[k]));
                yIndex[k] = (yIndex[k] < 0L) ? (0L)
                        : ((Height <= yIndex[k]) ? (Height - 1L) : (yIndex[k]));
                zIndex[k] = (zIndex[k] < 0L) ? (0L)
                        : ((Slice <= zIndex[k]) ? (Slice - 1L) : (zIndex[k]));
           }
	}
	else {
	   for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	   }
	}

	/* perform interpolation */
	interpolated = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
	}

	return(interpolated);

} /* end InterpolatedValue3DMirr */

extern double	InterpolatedGradX3DMirr (
			float	*Bcoeff,
			long	Width,	
			long	Height,	
			long	Slice,	
			double	x,	
			double	y,	
			double	z,	
			long	SplineDegree	
		)
{ /* begin InterpolatedGradX3DMirr */

	float	*p;
	double	gradxWeight[10], yWeight[10], zWeight[10];
	double	interpol_gradx;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree == 1L) {
	   /* perform the shift */
           x -= BASIS_SHIFT;
           y -= BASIS_SHIFT;
           z -= BASIS_SHIFT;
	
           /* compute the interpolation indexes */
           i = (long)floor(x);
           j = (long)floor(y);
           l = (long)floor(z);
           for (k = 0L; k <= 1; k++) {
                xIndex[k] = i++;
                yIndex[k] = j++;
                zIndex[k] = l++;
           }
	}
	else {
	   if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	   else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 1L:
			/* x */
        		w = x - (double)xIndex[0];
        		gradxWeight[0] = -1.0;
        		gradxWeight[1] = 1.0;
       	 		/* y */
        		w = y - (double)yIndex[0];
        		yWeight[0] = 1.0 - w;
        		yWeight[1] = w;
       	 		/* z */
        		w = z - (double)zIndex[0];
        		zWeight[0] = 1.0 - w;
        		zWeight[1] = w;
			break;
		case 2L:
			/* gradx */
			w = x - (double)xIndex[1];
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			break;
		case 3L:
			/* gradx */
			w = x - (double)xIndex[1];
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0+gradxWeight[0]-2.0*gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			break;
		case 4L:
			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] 
				- yWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			break;
		case 5L:
			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	if (SplineDegree == 1L) {
	   for (k = 0L; k <= 1; k++) {
                xIndex[k] = (xIndex[k] < 0L) ? (0L)
                        : ((Width <= xIndex[k]) ? (Width - 1L) : (xIndex[k]));
                yIndex[k] = (yIndex[k] < 0L) ? (0L)
                        : ((Height <= yIndex[k]) ? (Height - 1L) : (yIndex[k]));
                zIndex[k] = (zIndex[k] < 0L) ? (0L)
                        : ((Slice <= zIndex[k]) ? (Slice - 1L) : (zIndex[k]));
           }
	}
	else {
	   for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	   }
	}

	/* perform interpolation */
	interpol_gradx = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
		}
		interpol_gradx +=  zWeight[j] * w2;
	}

	return(interpol_gradx);

} /* end InterpolatedGradX3DMirr */

extern double	InterpolatedGradY3DMirr (
			float	*Bcoeff,	
			long	Width,	
			long	Height,		
			long	Slice,		
			double	x,		
			double	y,	
			double	z,	
			long	SplineDegree
		)
{ /* begin InterpolatedGradY3DMirr */

	float	*p;
	double	xWeight[10], gradyWeight[10], zWeight[10];
	double	interpol_grady;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree == 1L) {
	   /* perform the shift */
           x -= BASIS_SHIFT;
           y -= BASIS_SHIFT;
           z -= BASIS_SHIFT;
	
           /* compute the interpolation indexes */
           i = (long)floor(x);
           j = (long)floor(y);
           l = (long)floor(z);
           for (k = 0L; k <= 1; k++) {
                xIndex[k] = i++;
                yIndex[k] = j++;
                zIndex[k] = l++;
           }
	}
	else {
	   if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	   else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 1L:
			/* x */
        		w = x - (double)xIndex[0];
        		xWeight[0] = 1.0 - w;
        		xWeight[1] = w;
       	 		/* y */
        		w = y - (double)yIndex[0];
        		gradyWeight[0] = -1.0;
        		gradyWeight[1] = 1.0;
       	 		/* z */
        		w = z - (double)zIndex[0];
        		zWeight[0] = 1.0 - w;
        		zWeight[1] = w;
			break;
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* grady */
			w = y - (double)yIndex[1];
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* grady */
			w = y - (double)yIndex[1];
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	if (SplineDegree == 1L) {
	   for (k = 0L; k <= 1; k++) {
                xIndex[k] = (xIndex[k] < 0L) ? (0L)
                        : ((Width <= xIndex[k]) ? (Width - 1L) : (xIndex[k]));
                yIndex[k] = (yIndex[k] < 0L) ? (0L)
                        : ((Height <= yIndex[k]) ? (Height - 1L) : (yIndex[k]));
                zIndex[k] = (zIndex[k] < 0L) ? (0L)
                        : ((Slice <= zIndex[k]) ? (Slice - 1L) : (zIndex[k]));
           }
	}
	else {
	   for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	   }
	}

	/* perform interpolation */
	interpol_grady = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
			}
			w2 += gradyWeight[l] * w;
		}
		interpol_grady +=  zWeight[j] * w2;
	}

	return(interpol_grady);

} /* end InterpolatedGradY3DMirr */

extern double	InterpolatedGradZ3DMirr (
			float	*Bcoeff,
			long	Width,		
			long	Height,	
			long	Slice,		
			double	x,	
			double	y,	
			double	z,		
			long	SplineDegree
		)
{ /* begin InterpolatedGradZ3DMirr */

	float	*p;
	double	xWeight[10], yWeight[10], gradzWeight[10];
	double	interpol_gradz;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree == 1L) {
	   /* perform the shift */
           x -= BASIS_SHIFT;
           y -= BASIS_SHIFT;
           z -= BASIS_SHIFT;
	
           /* compute the interpolation indexes */
           i = (long)floor(x);
           j = (long)floor(y);
           l = (long)floor(z);
           for (k = 0L; k <= 1; k++) {
                xIndex[k] = i++;
                yIndex[k] = j++;
                zIndex[k] = l++;
           }
	}
	else {
	   if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	   else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	   }
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 1L:
			/* x */
        		w = x - (double)xIndex[0];
        		xWeight[0] = 1.0 - w;
        		xWeight[1] = w;
       	 		/* y */
        		w = y - (double)yIndex[0];
        		yWeight[0] = 1.0 - w;
        		yWeight[1] = w;
       	 		/* z */
        		w = z - (double)zIndex[0];
        		gradzWeight[0] = -1.0;
        		gradzWeight[1] = 1.0;
			break;
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* gradz */
			w = z - (double)zIndex[1];
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* gradz */
			w = z - (double)zIndex[1];
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* gradz */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	if (SplineDegree == 1L) {
	   for (k = 0L; k <= 1; k++) {
                xIndex[k] = (xIndex[k] < 0L) ? (0L)
                        : ((Width <= xIndex[k]) ? (Width - 1L) : (xIndex[k]));
                yIndex[k] = (yIndex[k] < 0L) ? (0L)
                        : ((Height <= yIndex[k]) ? (Height - 1L) : (yIndex[k]));
                zIndex[k] = (zIndex[k] < 0L) ? (0L)
                        : ((Slice <= zIndex[k]) ? (Slice - 1L) : (zIndex[k]));
           }
	}
	else {
	   for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	   }
	}

	/* perform interpolation */
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
		}
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	return(interpol_gradz);

} /* end InterpolatedGradZ3DMirr */


extern double	InterpolatedAll3DMirr (
			double	*gradx,	
			double	*grady,	
			double	*gradz,	
			float	*Bcoeff,	
			long	Width,		
			long	Height,	
			long	Slice,	
			double	x,
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedAll3DMirr */

	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	gradxWeight[10], gradyWeight[10], gradzWeight[10];
	double	interpolated, interpol_gradx, interpol_grady, interpol_gradz;
	double	w, w2, w4, t, t0, t1, wx, w2x, w2y;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* gradx */
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* grady */
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			/* gradz */
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* gradx */
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* grady */
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			/* gradz */
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] 
				- xWeight[4];

			/* gradx */
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* grady */
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			/* gradz */
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	}

	/* perform interpolation */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*gradx = interpol_gradx;
	*grady = interpol_grady;
	*gradz = interpol_gradz;
	
	return(interpolated);

} /* end InterpolatedAll3DMirr */

extern double	InterpolatedThree3DMirr (
			double	*value1,	
			double	*gradx1,	
			double	*grady1,	
			double	*gradz1,	
			double	*value2,	
			double	*gradx2,	
			double	*grady2,	
			double	*gradz2,	
			double	*value3,	
			double	*gradx3,	
			double	*grady3,	
			double	*gradz3,	
			float	*Bcoeff1,	
			float	*Bcoeff2,	
			float	*Bcoeff3,	
			long	Width,		
			long	Height,	
			long	Slice,	
			double	x,
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedThree3DMirr */

	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	gradxWeight[10], gradyWeight[10], gradzWeight[10];
	double	interpolated, interpol_gradx, interpol_grady, interpol_gradz;
	double	w, w2, w4, t, t0, t1, wx, w2x, w2y;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* gradx */
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* grady */
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			/* gradz */
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* gradx */
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* grady */
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			/* gradz */
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] 
				- xWeight[4];

			/* gradx */
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* grady */
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			/* gradz */
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
		zIndex[k] = (Slice == 1L) ? (0L) : ((zIndex[k] < 0L) ?
			(-zIndex[k] - Slice2 * ((-zIndex[k]) / Slice2))
			: (zIndex[k] - Slice2 * (zIndex[k] / Slice2)));
		if (Slice <= zIndex[k]) {
			zIndex[k] = Slice2 - zIndex[k];
		}
	}

	/* perform interpolation for Bcoeff1 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff1 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value1 = interpolated;
	*gradx1 = interpol_gradx;
	*grady1 = interpol_grady;
	*gradz1 = interpol_gradz;
	
	/* perform interpolation for Bcoeff2 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff2 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value2 = interpolated;
	*gradx2 = interpol_gradx;
	*grady2 = interpol_grady;
	*gradz2 = interpol_gradz;

	/* perform interpolation for Bcoeff3 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff3 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value3 = interpolated;
	*gradx3 = interpol_gradx;
	*grady3 = interpol_grady;
	*gradz3 = interpol_gradz;

	return(interpolated);

} /* end InterpolatedThree3DMirr */


extern double	InterpolatedValue3DZero (
			float	*Bcoeff,
			long	Width,	
			long	Height,
			long	Slice,	
			double	x,	
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedValue3DZero */

	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	interpolated;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
//	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
//		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;
			break;
		case 6L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= xWeight[0] / 720.0;
			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			xWeight[6] = 1.0 / 2.0 + w;
			xWeight[6] *= xWeight[6] * xWeight[6];
			xWeight[6] *= xWeight[6] / 720.0;
			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[6];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= yWeight[0] / 720.0;
			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			yWeight[6] = 1.0 / 2.0 + w;
			yWeight[6] *= yWeight[6] * yWeight[6];
			yWeight[6] *= yWeight[6] / 720.0;
			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[6];
			/* z */
			w = z - (double)zIndex[3];
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0] * zWeight[0];
			zWeight[0] *= zWeight[0] / 720.0;
			zWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			zWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			zWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			zWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			zWeight[6] = 1.0 / 2.0 + w;
			zWeight[6] *= zWeight[6] * zWeight[6];
			zWeight[6] *= zWeight[6] / 720.0;
			zWeight[5] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[6];
			break;
		case 7L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			xWeight[2] = (397.0 / 7.0 - w*(245.0/3.0 + w*(-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			xWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w * (-3.0 + w)*(-3.0 + w2)))))/144.0;
			xWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			xWeight[7] = w2;
			xWeight[7] *= xWeight[7] * xWeight[7];
			xWeight[7] *= w / 5040.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			yWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			yWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				*(19.0 + w*(-3.0 + w) * (-3.0 + w2))))) / 144.0;
			yWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			yWeight[7] = w2;
			yWeight[7] *= yWeight[7] * yWeight[7];
			yWeight[7] *= w / 5040.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7];
			/* z */
			w = z - (double)zIndex[3];
			zWeight[0] = 1.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] * zWeight[0];
			zWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			zWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			zWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			zWeight[3] = (2416.0/35.0 + w2*(-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			zWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			zWeight[5] = (40.0/7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 
				+ w * (-2.0 + w)))))) / 240.0;
			zWeight[7] = w2;
			zWeight[7] *= zWeight[7] * zWeight[7];
			zWeight[7] *= w / 5040.0;
			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[2] - zWeight[3]
				- zWeight[4] - zWeight[5] - zWeight[7];
			break;
		case 8L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] / 40320.0;
			w2 = w * w;
			xWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0/16.0 + w * (-15.0/4.0 + w*(9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			xWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 
				+ w * (-1.0 + w)))))))) / 720.0;
			xWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			xWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			xWeight[7] = (39.0 / 16.0 - w * (-6.0 
				+ w * (-9.0 / 2.0 + w2)))
				*(21.0/16.0 + w*(15.0 / 4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			xWeight[8] = 1.0 / 2.0 + w;
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8] / 40320.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[7] - xWeight[8];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] / 40320.0;
			w2 = w * w;
			yWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0/4.0 + w*(9.0/2.0 + w
				* (-3.0 + w)))) / 5040.0;
			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			yWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			yWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			yWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			yWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(15.0/4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			yWeight[8] = 1.0 / 2.0 + w;
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8] / 40320.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[7] - yWeight[8];
			/* z */
			w = z - (double)zIndex[4];
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] / 40320.0;
			w2 = w * w;
			zWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0/2.0 + w2)))
				* (21.0 / 16.0 + w*(-15.0/4.0 + w*(9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			zWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				*(2275.0/16.0 + w*(-487.0/8.0 + w*(-85.0/8.0 + w
				* (41.0 / 2.0 + w * (-5.0 
				+ w * (-2.0 + w)))))))) / 1440.0;
			zWeight[3] = (310661.0 / 1792.0 - w*(14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(53.0 / 4.0 + w * (-8.0 
				+ w * (-1.0 + w)))))))) / 720.0;
			zWeight[4] = (2337507.0/8960.0 + w2*(-2601.0/16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			zWeight[5] = (310661.0 / 1792.0 - w * (-14219.0/64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w*(-53.0 / 4.0 + w * (-8.0 
				+ w * (1.0 + w)))))))) / 720.0;
			zWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w * (15.0/4.0 + w*(9.0/2.0 + w
				* (3.0 + w)))) / 5040.0;
			zWeight[8] = 1.0 / 2.0 + w;
			zWeight[8] *= zWeight[8];
			zWeight[8] *= zWeight[8];
			zWeight[8] *= zWeight[8] / 40320.0;
			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[5] 
				- zWeight[7] - zWeight[8];
			break;
		case 9L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 
				+ w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[6] - yWeight[7] - yWeight[9];
			/* z */
			w = z - (double)zIndex[4];
			zWeight[0] = 1.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] * (1.0 - w) / 362880.0;
			zWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			zWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			zWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			zWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			zWeight[5] = (44117.0 / 63.0 + w*(809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			zWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			zWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 
				+ w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			zWeight[9] = w2 * w2;
			zWeight[9] *= zWeight[9] * w / 362880.0;
			zWeight[8] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] 
				- zWeight[3] - zWeight[4] - zWeight[5] 
				- zWeight[6] - zWeight[7] - zWeight[9];
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	interpolated = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
           if ((zIndex[j] >= 0) && (zIndex[j] < Slice)) {
              w2 = 0.0;
              for (l = 0L; l <= SplineDegree; l++) {
                 if ((yIndex[l] >= 0) && (yIndex[l] < Height)) {
                    p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width + zIndex[j] 
				* Width * Height);
                    w = 0.0;
                    for (i = 0L; i <= SplineDegree; i++) {
                       if ((xIndex[i] >= 0) && (xIndex[i] < Width)) {
                          w += xWeight[i] * p[xIndex[i]];
                       }
                    }
                    w2 += yWeight[l] * w;
                 }
              }
              interpolated += zWeight[j] * w2;
           }
        }

	return(interpolated);

} /* end InterpolatedValue3DZero */

extern double	InterpolatedGradX3DZero (
			float	*Bcoeff,
			long	Width,	
			long	Height,	
			long	Slice,	
			double	x,	
			double	y,	
			double	z,	
			long	SplineDegree	
		)
{ /* begin InterpolatedGradX3DZero */

	float	*p;
	double	gradxWeight[10], yWeight[10], zWeight[10];
	double	interpol_gradx;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
//	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
//		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* gradx */
			w = x - (double)xIndex[1];
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			break;
		case 3L:
			/* gradx */
			w = x - (double)xIndex[1];
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0+gradxWeight[0]-2.0*gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			break;
		case 4L:
			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] 
				- yWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			break;
		case 5L:
			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	interpol_gradx = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
           if ((zIndex[j] >= 0) && (zIndex[j] < Slice)) {
              w2 = 0.0;
              for (l = 0L; l <= SplineDegree; l++) {
                 if ((yIndex[l] >= 0) && (yIndex[l] < Height)) {
                    p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width + zIndex[j] 
				* Width * Height);
                    w = 0.0;
                    for (i = 0L; i <= SplineDegree; i++) {
                       if ((xIndex[i] >= 0) && (xIndex[i] < Width)) {
                          w += gradxWeight[i] * p[xIndex[i]];
                       }
                    }
                 w2 += yWeight[l] * w;
                 }
              }
              interpol_gradx +=  zWeight[j] * w2;
           }
        }

	return(interpol_gradx);

} /* end InterpolatedGradX3DZero */

extern double	InterpolatedGradY3DZero (
			float	*Bcoeff,	
			long	Width,	
			long	Height,		
			long	Slice,		
			double	x,		
			double	y,	
			double	z,	
			long	SplineDegree
		)
{ /* begin InterpolatedGradY3DZero */

	float	*p;
	double	xWeight[10], gradyWeight[10], zWeight[10];
	double	interpol_grady;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
//	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
//		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* grady */
			w = y - (double)yIndex[1];
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* grady */
			w = y - (double)yIndex[1];
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	interpol_grady = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
           if ((zIndex[j] >= 0) && (zIndex[j] < Slice)) {
              w2 = 0.0;
              for (l = 0L; l <= SplineDegree; l++) {
                 if ((yIndex[l] >= 0) && (yIndex[l] < Height)) {
                    p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width + zIndex[j] 
				* Width * Height);
                    w = 0.0;
                    for (i = 0L; i <= SplineDegree; i++) {
                       if ((xIndex[i] >= 0) && (xIndex[i] < Width)) {
                          w += xWeight[i] * p[xIndex[i]];
                       }
                    }
                    w2 += gradyWeight[l] * w;
                 }
              }
              interpol_grady +=  zWeight[j] * w2;
           }
        }

	return(interpol_grady);

} /* end InterpolatedGradY3DZero */

extern double	InterpolatedGradZ3DZero (
			float	*Bcoeff,
			long	Width,		
			long	Height,	
			long	Slice,		
			double	x,	
			double	y,	
			double	z,		
			long	SplineDegree
		)
{ /* begin InterpolatedGradZ3DZero */

	float	*p;
	double	xWeight[10], yWeight[10], gradzWeight[10];
	double	interpol_gradz;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10], zIndex[10];
//	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
//		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* gradz */
			w = z - (double)zIndex[1];
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* gradz */
			w = z - (double)zIndex[1];
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* gradz */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
           if ((zIndex[j] >= 0) && (zIndex[j] < Slice)) {
              w2 = 0.0;
              for (l = 0L; l <= SplineDegree; l++) {
                 if ((yIndex[l] >= 0) && (yIndex[l] < Height)) {
                    p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width + zIndex[j] 
				* Width * Height);
                    w = 0.0;
                    for (i = 0L; i <= SplineDegree; i++) {
                       if ((xIndex[i] >= 0) && (xIndex[i] < Width)) {
                          w += xWeight[i] * p[xIndex[i]];
                       }
                    }
                    w2 += yWeight[l] * w;
                 }
              }
              interpol_gradz +=  gradzWeight[j] * w2;
           }
        }

	return(interpol_gradz);

} /* end InterpolatedGradZ3DZero */


extern double	InterpolatedAll3DZero (
			double	*gradx,	
			double	*grady,	
			double	*gradz,	
			float	*Bcoeff,	
			long	Width,		
			long	Height,	
			long	Slice,	
			double	x,
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedAll3DZero */

	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	gradxWeight[10], gradyWeight[10], gradzWeight[10];
	double	interpolated, interpol_gradx, interpol_grady, interpol_gradz;
	double	w, w2, w4, t, t0, t1, wx, w2x, w2y;
	long	xIndex[10], yIndex[10], zIndex[10];
//	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L, 
//		Slice2 = 2L * Slice - 2L;
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* gradx */
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* grady */
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			/* gradz */
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* gradx */
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* grady */
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			/* gradz */
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] 
				- xWeight[4];

			/* gradx */
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* grady */
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			/* gradz */
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
           if ((zIndex[j] >= 0) && (zIndex[j] < Slice)) {
	      w2 = 0.0; w2x = 0.0; w2y = 0.0;
              for (l = 0L; l <= SplineDegree; l++) {
                 if ((yIndex[l] >= 0) && (yIndex[l] < Height)) {
                    p = Bcoeff + (ptrdiff_t)(yIndex[l] * Width + zIndex[j] 
				* Width * Height);
		    w = 0.0; wx = 0.0;
                    for (i = 0L; i <= SplineDegree; i++) {
                       if ((xIndex[i] >= 0) && (xIndex[i] < Width)) {
                          w += xWeight[i] * p[xIndex[i]];
			  wx += gradxWeight[i] * p[xIndex[i]];
                       }
                    }
		    w2 += yWeight[l] * w;
		    w2x += yWeight[l] * wx;
		    w2y += gradyWeight[l] * w;
                 }
              }
	      interpolated += zWeight[j] * w2;
	      interpol_gradx += zWeight[j] * w2x;
	      interpol_grady += zWeight[j] * w2y;
	      interpol_gradz += gradzWeight[j] * w2;
           }
        }

	*gradx = interpol_gradx;
	*grady = interpol_grady;
	*gradz = interpol_gradz;
	
	return(interpolated);

} /* end InterpolatedAll3DZero */

extern void	InterpolatedThree3DZero (
			double	*value1,	
			double	*gradx1,	
			double	*grady1,	
			double	*gradz1,	
			double	*value2,	
			double	*gradx2,	
			double	*grady2,	
			double	*gradz2,	
			double	*value3,	
			double	*gradx3,	
			double	*grady3,	
			double	*gradz3,	
			float	*Bcoeff1,	
			float	*Bcoeff2,	
			float	*Bcoeff3,	
			long	Width,		
			long	Height,	
			long	Slice,	
			double	x,
			double	y,
			double	z,
			long	SplineDegree
		)
{ /* begin InterpolatedThree3DZero */

(void) Slice; // jf
	float	*p;
	double	xWeight[10], yWeight[10], zWeight[10];
	double	gradxWeight[10], gradyWeight[10], gradzWeight[10];
	double	interpolated, interpol_gradx, interpol_grady, interpol_gradz;
	double	w, w2, w4, t, t0, t1, wx, w2x, w2y;
	long	xIndex[10], yIndex[10], zIndex[10];
	long	i, j, l, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		l = (long)floor(z) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		l = (long)floor(z + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = l++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];

			/* gradx */
			gradxWeight[1] = - 2.0 * w;
			gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
			gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];

			/* grady */
			gradyWeight[1] = - 2.0 * w;
			gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
			gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[1] = 3.0 / 4.0 - w * w;
			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];

			/* gradz */
			gradzWeight[1] = - 2.0 * w;
			gradzWeight[2] = (1.0 / 2.0) * (1.0 - gradzWeight[1]);
			gradzWeight[0] = - gradzWeight[1] - gradzWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];

			/* gradx */
			gradxWeight[3] = (1.0 / 2.0) * w * w;
			gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
			gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
			gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 
				- gradxWeight[3];

			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];

			/* grady */
			gradyWeight[3] = (1.0 / 2.0) * w * w;
			gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
			gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
			gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];

			/* z */
			w = z - (double)zIndex[1];
			zWeight[3] = (1.0 / 6.0) * w * w * w;
			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- zWeight[3];
			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];

			/* gradz */
			gradzWeight[3] = (1.0 / 2.0) * w * w;
			gradzWeight[0] = w - 1.0 / 2.0 - gradzWeight[3];
			gradzWeight[2] = 1.0 + gradzWeight[0] 
				- 2.0 * gradzWeight[3];
			gradzWeight[1] = - gradzWeight[0] - gradzWeight[2] 
				- gradzWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] 
				- xWeight[4];

			/* gradx */
			gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradxWeight[1] = t1 + t0;
			gradxWeight[3] = t1 - t0;
			gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
			gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];

			/* grady */
			gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradyWeight[1] = t1 + t0;
			gradyWeight[3] = t1 - t0;
			gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
			gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			zWeight[0] = 1.0 / 2.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			zWeight[1] = t1 + t0;
			zWeight[3] = t1 - t0;
			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] 
				- zWeight[3] - zWeight[4];

			/* gradz */
			gradzWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
			t = (1.0 / 2.0) * w;
			t0 = t * w - 11.0 / 24.0;
			t1 = t - (2.0 / 3.0) * w2 * w;
			gradzWeight[1] = t1 + t0;
			gradzWeight[3] = t1 - t0;
			gradzWeight[4] = gradzWeight[0] + t0 + (1.0 / 2.0);
			gradzWeight[2] = - gradzWeight[0] - gradzWeight[1] 
				- gradzWeight[3] - gradzWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;

			/* gradx */
			w = x - (double)xIndex[2];
			w2 = w * w;
			gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradxWeight[2] = t0 + t1;
			gradxWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradxWeight[1] = t0 + t1;
			gradxWeight[4] = t0 - t1;

			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;

			/* grady */
			w = y - (double)yIndex[2];
			w2 = w * w;
			gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradyWeight[2] = t0 + t1;
			gradyWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradyWeight[1] = t0 + t1;
			gradyWeight[4] = t0 - t1;

			/* z */
			w = z - (double)zIndex[2];
			w2 = w * w;
			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- zWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			zWeight[2] = t0 + t1;
			zWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			zWeight[1] = t0 + t1;
			zWeight[4] = t0 - t1;

			/* grady */
			w = z - (double)zIndex[2];
			w2 = w * w;
			gradzWeight[5] = (1.0 / 24.0) * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			gradzWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradzWeight[5];
			t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
			t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
			gradzWeight[2] = t0 + t1;
			gradzWeight[3] = t0 - t1;
			t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
			t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
			gradzWeight[1] = t0 + t1;
			gradzWeight[4] = t0 - t1;

			break;
		default:
			printf("Invalid spline degree\n");
			return;
	}

	/* perform interpolation for Bcoeff1 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff1 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value1 = interpolated;
	*gradx1 = interpol_gradx;
	*grady1 = interpol_grady;
	*gradz1 = interpol_gradz;
	
	/* perform interpolation for Bcoeff2 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff2 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value2 = interpolated;
	*gradx2 = interpol_gradx;
	*grady2 = interpol_grady;
	*gradz2 = interpol_gradz;

	/* perform interpolation for Bcoeff3 */
	interpolated = 0.0;
	interpol_gradx = 0.0;
	interpol_grady = 0.0;
	interpol_gradz = 0.0;

	for (j = 0L; j <= SplineDegree; j++) {
		w2 = 0.0; w2x = 0.0; w2y = 0.0;
		for (l = 0L; l <= SplineDegree; l++) {
			p = Bcoeff3 + (ptrdiff_t)(yIndex[l] * Width 
				+ zIndex[j] * Width * Height);
			w = 0.0; wx = 0.0;
			for (i = 0L; i <= SplineDegree; i++) {
				w += xWeight[i] * p[xIndex[i]];
				wx += gradxWeight[i] * p[xIndex[i]];
			}
			w2 += yWeight[l] * w;
			w2x += yWeight[l] * wx;
			w2y += gradyWeight[l] * w;
		}
		interpolated += zWeight[j] * w2;
		interpol_gradx +=  zWeight[j] * w2x;
		interpol_grady +=  zWeight[j] * w2y;
		interpol_gradz +=  gradzWeight[j] * w2;
	}

	*value3 = interpolated;
	*gradx3 = interpol_gradx;
	*grady3 = interpol_grady;
	*gradz3 = interpol_gradz;

} /* end InterpolatedThree3DZero */

