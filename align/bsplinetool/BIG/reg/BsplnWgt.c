#include		<stdlib.h>
#include		<math.h>

#include		"BsplnWgt.h"

/************************************************************************/
/* FUNCTION: BsplnDDWght												*/
/************************************************************************/
double				BsplnDDWght		(int				degree,
									 long				i,
									 double				x) {

		switch (degree) {
		  case 0:
			return(BsplnDDWght0(i, x));
		  case 1:
			return(BsplnDDWght1(i, x));
		  case 2:
			return(BsplnDDWght2(i, x));
		  case 3:
			return(BsplnDDWght3(i, x));
		  case 4:
			return(BsplnDDWght4(i, x));
		  case 5:
			return(BsplnDDWght5(i, x));
		  case 6:
			return(BsplnDDWght6(i, x));
		  case 7:
			return(BsplnDDWght7(i, x));
		  default:
			return(0.0);
		}
} /* End of BsplnDDWght */

/************************************************************************/
/* FUNCTION: BsplnDDWght0												*/
/************************************************************************/
double				BsplnDDWght0	(long				i,
									 double				x) {

		return((fabs(x - (double)i) == 0.5) ? (HUGE_VAL) : (0.0));
} /* End of BsplnDDWght0 */

/************************************************************************/
/* FUNCTION: BsplnDDWght1												*/
/************************************************************************/
double				BsplnDDWght1	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		return((x == 0.0) ? (-HUGE_VAL) : ((x == 1.0) ? (HUGE_VAL) : (0.0)));
} /* End of BsplnDDWght1 */

/************************************************************************/
/* FUNCTION: BsplnDDWght2												*/
/************************************************************************/
double				BsplnDDWght2	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5)
		  return(-2.0);
		if (x == 0.5)
		  return(-0.5);
		if (x < 1.5)
		  return(1.0);
		if (x == 1.5)
		  return(0.5);
		return(0.0);
} /* End of BsplnDDWght2 */

/************************************************************************/
/* FUNCTION: BsplnDDWght3												*/
/************************************************************************/
double				BsplnDDWght3	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 1.0)
		  return(x * 3.0 - 2.0);
		if (x < 2.0)
		  return(2.0 - x);
		return(0.0);
} /* End of BsplnDDWght3 */

/************************************************************************/
/* FUNCTION: BsplnDDWght4												*/
/************************************************************************/
double				BsplnDDWght4	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5)
		  return(3.0 * x * x - 1.25);
		if (x < 1.5)
		  return(x * (5.0 - x * 2.0) - 2.5);
		if (x < 2.5) {
		  x -= 2.5;
		  return(x * x * 0.5);
		}
		return(0.0);
} /* End of BsplnDDWght4 */

/************************************************************************/
/* FUNCTION: BsplnDDWght5												*/
/************************************************************************/
double				BsplnDDWght5	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 1.0)
		  return(x * x * (3.0 - x * (5.0 / 3.0)) - 1.0);
		if (x < 2.0)
		  return(x * (x * (x * (5.0 / 6.0) - 4.5) + 7.5) - 3.5);
		if (x < 3.0) {
		  x -= 3.0;
		  return(x * x * x * (-1.0 / 6.0));
		}
		return(0.0);
} /* End of BsplnDDWght5 */

/************************************************************************/
/* FUNCTION: BsplnDDWght6												*/
/************************************************************************/
double				BsplnDDWght6	(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5) {
		  x *= x;
		  return(x * (7.0 / 4.0 - x * (5.0 / 6.0)) - 77.0 / 96.0);
		}
		if (x < 1.5)
		  return(x * (x * (x * (x * 0.625 - 35.0 / 12.0) + 3.9375) - 35.0 / 48.0)
			- 0.7109375);
		if (x < 2.5)
		  return(x * (x * (x * (7.0 / 3.0 - x * 0.25) - 7.875) + 133.0 / 12.0)
			- 5.140625);
		if (x < 3.5) {
		  x -= 3.5;
		  x *= x;
		  return(x * x * (1.0 / 24.0));
		}
		return(0.0);
} /* End of BsplnDDWght6 */

/************************************************************************/
/* FUNCTION: BsplnDDWght7												*/
/************************************************************************/
double				BsplnDDWght7	(long				i,
									 double				x) {

		double				g;

		x = fabs(x - (double)i);
		if (x < 1.0) {
		  g = x * x;
		  return(g * (g * (x * (7.0 / 24.0) - 5.0 / 6.0) + 4.0 / 3.0)
			- 2.0 / 3.0);
		}
		if (x < 2.0)
		  return(x * (x * (x * (x * (1.5 - x * 0.175) - 14.0 / 3.0) + 6.0)
			- 7.0 / 3.0) - 0.2);
		if (x < 3.0)
		  return(x * (x * (x * (x * (x * (7.0 / 120.0) - 5.0 / 6.0) + 14.0 / 3.0)
			- 38.0 / 3.0) + 49.0 / 3.0) - 23.0 / 3.0);
		if (x < 4.0) {
		  x -= 4.0;
		  g = x * x;
		  return(x * g * g * (-1.0 / 120.0));
		}
		return(0.0);
} /* End of BsplnDDWght7 */

/************************************************************************/
/* FUNCTION: BsplnDWght													*/
/************************************************************************/
double				BsplnDWght		(int				degree,
									 long				i,
									 double				x) {

		switch (degree) {
		  case 0:
			return(BsplnDWght0(i, x));
		  case 1:
			return(BsplnDWght1(i, x));
		  case 2:
			return(BsplnDWght2(i, x));
		  case 3:
			return(BsplnDWght3(i, x));
		  case 4:
			return(BsplnDWght4(i, x));
		  case 5:
			return(BsplnDWght5(i, x));
		  case 6:
			return(BsplnDWght6(i, x));
		  case 7:
			return(BsplnDWght7(i, x));
		  default:
			return(0.0);
		}
} /* End of BsplnDWght */

/************************************************************************/
/* FUNCTION: BsplnDWght0												*/
/************************************************************************/
double				BsplnDWght0		(long				i,
									 double				x) {

		x -= (double)i;
		return((x == -0.5) ? (HUGE_VAL) : ((x == 0.5) ? (-HUGE_VAL) : (0.0)));
} /* End of BsplnDWght0 */

/************************************************************************/
/* FUNCTION: BsplnDWght1												*/
/************************************************************************/
double				BsplnDWght1		(long				i,
									 double				x) {

		double				f;

		x -= (double)i;
		if (x == 0.0)
		  return(0.0);
		f = fabs(x);
		if (f < 1.0)
		  f = -1.0;
		else if (f == 1.0)
		  f = -0.5;
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght1 */

/************************************************************************/
/* FUNCTION: BsplnDWght2												*/
/************************************************************************/
double				BsplnDWght2		(long				i,
									 double				x) {

		double				f;

		x -= (double)i;
		f = fabs(x);
		if (f < 0.5)
		  f *= -2.0;
		else if (f < 1.5)
		  f -= 1.5;
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght2 */

/************************************************************************/
/* FUNCTION: BsplnDWght3												*/
/************************************************************************/
double				BsplnDWght3		(long				i,
									 double				x) {

		double				f;

		x -= (double)i;
		f = fabs(x);
		if (f < 1.0)
		  f *= f * 1.5 - 2.0;
		else if (f < 2.0) {
		  f = 2.0 - f;
		  f *= -0.5 * f;
		}
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght3 */

/************************************************************************/
/* FUNCTION: BsplnDWght4												*/
/************************************************************************/
double				BsplnDWght4		(long				i,
									 double				x) {

		double				f;

		x -= (double)i;
		f = fabs(x);
		if (f < 0.5)
		  f *= f * f - 1.25;
		else if (f < 1.5)
		  f = f * (f * (2.5 - f * (2.0 / 3.0)) - 2.5) + 5.0 / 24.0;
		else if (f < 2.5) {
		  f -= 2.5;
		  f *= f * f * (1.0 / 6.0);
		}
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght4 */

/************************************************************************/
/* FUNCTION: BsplnDWght5												*/
/************************************************************************/
double				BsplnDWght5		(long				i,
									 double				x) {

		double				f;

		x -= (double)i;
		f = fabs(x);
		if (f < 1.0)
		  f *= f * f * (1.0 - f * (5.0 / 12.0)) - 1.0;
		else if (f < 2.0)
		  f = f * (f * (f * (f * (5.0 / 24.0) - 1.5) + 3.75) - 3.5) + 0.625;
		else if (f < 3.0) {
		  f -= 3.0;
		  f *= f;
		  f *= f * (-1.0 / 24.0);
		}
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght5 */

/************************************************************************/
/* FUNCTION: BsplnDWght6												*/
/************************************************************************/
double				BsplnDWght6		(long				i,
									 double				x) {

		double				f, g;

		x -= (double)i;
		f = fabs(x);
		if (f < 0.5) {
		  g = f * f;
		  f *= g * (7.0 / 12.0 - g * (1.0 / 6.0)) - 77.0 / 96.0;
		}
		else if (f < 1.5)
		  f = f * (f * (f * (f * (f * 0.125 - 35.0 / 48.0) + 1.3125) - 35.0 / 96.0)
			- 0.7109375) - 7.0 / 768.0;
		else if (f < 2.5)
		  f = f * (f * (f * (f * (7.0 / 12.0 - f * 0.05) - 2.625) + 133.0 / 24.0)
			- 5.140625) + 1267.0 / 960.0;
		else if (f < 3.5) {
		  f -= 3.5;
		  g = f * f;
		  f *= g * g * (1.0 / 120.0);
		}
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght6 */

/************************************************************************/
/* FUNCTION: BsplnDWght7												*/
/************************************************************************/
double				BsplnDWght7		(long				i,
									 double				x) {

		double				f, g;

		x -= (double)i;
		f = fabs(x);
		if (f < 1.0) {
		  g = f * f;
		  f *= g * (g * (f * (7.0 / 144.0) - 1.0 / 6.0) + 4.0 / 9.0) - 2.0 / 3.0;
		}
		else if (f < 2.0)
		  f = f * (f * (f * (f * (f * (0.3 - f * (7.0 / 240.0)) - 7.0 / 6.0) + 2.0)
			- 7.0 / 6.0) - 0.2) - 7.0 / 90.0;
		else if (f < 3.0)
		  f = f * (f * (f * (f * (f * (f * (7.0 / 720.0) - 1.0 / 6.0) + 7.0 / 6.0)
			- 38.0 / 9.0) + 49.0 / 6.0) - 23.0 / 3.0) + 217.0 / 90.0;
		else if (f < 4.0) {
		  f -= 4.0;
		  f *= f;
		  f *= f * f * (-1.0 / 720.0);
		}
		else
		  return(0.0);
		return((x > 0.0) ? (f) : (-f));
} /* End of BsplnDWght7 */

/************************************************************************/
/* FUNCTION: BsplnWght													*/
/************************************************************************/
double				BsplnWght		(int				degree,
									 long				i,
									 double				x) {

		switch (degree) {
		  case 0:
			return(BsplnWght0(i, x));
		  case 1:
			return(BsplnWght1(i, x));
		  case 2:
			return(BsplnWght2(i, x));
		  case 3:
			return(BsplnWght3(i, x));
		  case 4:
			return(BsplnWght4(i, x));
		  case 5:
			return(BsplnWght5(i, x));
		  case 6:
			return(BsplnWght6(i, x));
		  case 7:
			return(BsplnWght7(i, x));
		  default:
			return(0.0);
		}
} /* End of BsplnWght */

/************************************************************************/
/* FUNCTION: BsplnWght0													*/
/************************************************************************/
double				BsplnWght0		(long				i,
									 double				x) {

		x -= (double)i;
		return((x >= 0.5) ? (0.0) : ((x < -0.5) ? (0.0) : (1.0)));
} /* End of BsplnWght0 */

/************************************************************************/
/* FUNCTION: BsplnWght1													*/
/************************************************************************/
double				BsplnWght1		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		return((x > 1.0) ? (0.0) : (1.0 - x));
} /* End of BsplnWght1 */

/************************************************************************/
/* FUNCTION: BsplnWght2													*/
/************************************************************************/
double				BsplnWght2		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5)
		  return(0.75 - x * x);
		if (x < 1.5) {
		  x -= 1.5;
		  return(0.5 * x * x);
		}
		return(0.0);
} /* End of BsplnWght2 */

/************************************************************************/
/* FUNCTION: BsplnWght3													*/
/************************************************************************/
double				BsplnWght3		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 1.0)
		  return(0.5 * x * x * (x - 2.0) + 2.0 / 3.0);
		if (x < 2.0) {
		  x -= 2.0;
		  return(x * x * x * (-1.0 / 6.0));
		}
		return(0.0);
} /* End of BsplnWght3 */

/************************************************************************/
/* FUNCTION: BsplnWght4													*/
/************************************************************************/
double				BsplnWght4		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5) {
		  x *= x;
		  return(x * (x * 0.25 - 0.625) + 115.0 / 192.0);
		}
		if (x < 1.5)
		  return(x * (x * (x * (5.0 / 6.0 - x * (1.0 / 6.0)) - 1.25) + 5.0 / 24.0)
			+ 55.0 / 96.0);
		if (x < 2.5) {
		  x -= 2.5;
		  x *= x;
		  return(x * x * (1.0 / 24.0));
		}
		return(0.0);
} /* End of BsplnWght4 */

/************************************************************************/
/* FUNCTION: BsplnWght5													*/
/************************************************************************/
double				BsplnWght5		(long				i,
									 double				x) {

		double				f;

		x = fabs(x - (double)i);
		if (x < 1.0) {
		  f = x * x;
		  return(f * (f * (0.25 - x * (1.0 / 12.0)) - 0.5) + 0.55);
		}
		if (x < 2.0)
		  return(x * (x * (x * (x * (x * (1.0 / 24.0) - 0.375) + 1.25) - 1.75) + 0.625)
			+ 0.425);
		if (x < 3.0) {
		  f = 3.0 - x;
		  x = f * f;
		  return(f * x * x * (1.0 / 120.0));
		}
		return(0.0);
} /* End of BsplnWght5 */

/************************************************************************/
/* FUNCTION: BsplnWght6													*/
/************************************************************************/
double				BsplnWght6		(long				i,
									 double				x) {

		x = fabs(x - (double)i);
		if (x < 0.5) {
		  x *= x;
		  return(x * (x * (7.0 / 48.0 - x * (1.0 / 36.0)) - 77.0 / 192.0) + 5887.0 / 11520.0);
		}
		if (x < 1.5)
		  return(x * (x * (x * (x * (x * (x * (1.0 / 48.0) - 7.0 / 48.0) + 0.328125)
			- 35.0 / 288.0) - 91.0 / 256.0) - 7.0 / 768.0) + 7861.0 / 15360.0);
		if (x < 2.5)
		  return(x * (x * (x * (x * (x * (7.0 / 60.0 - x * (1.0 / 120.0)) - 0.65625)
			+ 133.0 / 72.0) - 2.5703125) + 1267.0 / 960.0) + 1379.0 / 7680.0);
		if (x < 3.5) {
		  x -= 3.5;
		  x *= x * x;
		  return(x * x * (1.0 / 720.0));
		}
		return(0.0);
} /* End of BsplnWght6 */

/************************************************************************/
/* FUNCTION: BsplnWght7													*/
/************************************************************************/
double				BsplnWght7		(long				i,
									 double				x) {

		double				f;

		x = fabs(x - (double)i);
		if (x < 1.0) {
		  f = x * x;
		  return(f * (f * (f * (x * (1.0 / 144.0) - 1.0 / 36.0) + 1.0 / 9.0) - 1.0 / 3.0)
			+ 151.0 / 315.0);
		}
		if (x < 2.0)
		  return(x * (x * (x * (x * (x * (x * (0.05 - x * (1.0 / 240.0)) - 7.0 / 30.0) + 0.5)
			- 7.0 / 18.0) - 0.1) - 7.0 / 90.0) + 103.0 / 210.0);
		if (x < 3.0)
		  return(x * (x * (x * (x * (x * (x * (x * (1.0 / 720.0) - 1.0 / 36.0) + 7.0 / 30.0)
			- 19.0 / 18.0) + 49.0 / 18.0) - 23.0 / 6.0) + 217.0 / 90.0) - 139.0 / 630.0);
		if (x < 4.0) {
		  f = 4.0 - x;
		  x = f * f * f;
		  return(x * x * f * (1.0 / 5040.0));
		}
		return(0.0);
} /* End of BsplnWght7 */

/************************************************************************/
/* FUNCTION: fold														*/
/************************************************************************/
long				fold			(long				i,
									 long				n) {

		ldiv_t				modOp;
		long				n2;

		i = labs(i);
		if (i < n)
		  return(i);
		if (n == 1L)
		  return(0L);
		n2 = (n << 1L) - 2L;
		modOp = ldiv(i, n2);
		return((modOp.rem < n) ? (modOp.rem) : (n2 - modOp.rem));
} /* End of fold */
