#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<string.h>
#include		<math.h>
#include		<float.h>
#include		"phil.h"

#include		"convolve.h"

/************************************************************************/
/* +=============================+										*/
/* | FUNCTION: firConvolveFinite |										*/
/* +=============================+										*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete convolution									*/
/*				outPtr[x] = SUM(k): inPtr[x] * kernel[k - x]			*/
/*																		*/
/* Conventions:	Both the signal an the kernel have a finite support		*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*				The length of the kernel is nk							*/
/*				The kernel origin (hot spot) is at index [(nk - 1) / 2]	*/
/*				The kernel has to be given straight (k[x], *not* k[-x])	*/
/*																		*/
/************************************************************************/
void				firConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk) {

		double				*p, *q;
		double				sum;
		long				i, j;
		long				k, kp, km;

		kernel += (ptrdiff_t)(nk - 1L);
		k = -shift - nk / 2L;
		for (j = -n; (j < 0L); k++, j++) {
		  kp = (k > 0L) ? (k) : (0L);
		  km = k - kp;
		  p = inPtr + (ptrdiff_t)kp;
		  q = kernel + (ptrdiff_t)km;
		  sum = 0.0;
		  for (i = (n <= (k + nk)) ? (kp - n) : (kp - k - nk); (i < 0L); i++)
			sum += *p++ * *q--;
		  *outPtr++ = sum;
		}
} /* End of firConvolveFinite */

/************************************************************************/
/* +=============================+										*/
/* | FUNCTION: firConvolveMirror |										*/
/* +=============================+										*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete convolution									*/
/*				outPtr[x] = SUM(k): inPtr[x] * kernel[k - x]			*/
/*																		*/
/* Conventions:	The kernel has a finite support of size nk				*/
/*				The signal has an infinite, periodic (2 n - 1) support	*/
/*				The input is described by n values; the missing values	*/
/*					are provided by mirror symetry						*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*				The kernel origin (hot spot) is at index [(nk - 1) / 2]	*/
/*				The kernel has to be given straight (k[x], *not* k[-x])	*/
/*																		*/
/************************************************************************/
void				firConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk) {

		double				*p, *q;
		double				sum;
		long				Ok, n2;
		long				i, j, k, l;

		kernel += (ptrdiff_t)(nk - 1L);
		Ok = (nk - 1L) / 2L;
		n2 = 2L * (n - 1L);
		shift = ((nk & 1L) != 0L) ? (shift) : (shift + 1L);
		if (nk <= n) {
		  k = (n == 1L) ? (0L) : (((shift + Ok) > 0L) ? (n2 * ((n2 - 1L - (-Ok - shift)) / n2)
			- Ok - shift) : (-Ok - shift - n2 * ((-Ok - shift) / n2)));
		  for (j = (k < n) ? (n - k) : (0L); (k < n); k++) {
			p = inPtr + (ptrdiff_t)k;
			q = kernel;
			sum = 0.0;
			for (i = k, l = ((nk + k) < (n - 1L)) ? (nk + k) : (n - 1L); (i < l); i++)
			  sum += *p++ * *q--;
			for (l = (n2 < (nk + k)) ? (n2) : (nk + k); (i < l); i++)
			  sum += *p-- * *q--;
			for (l = nk + k; (i < l); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		  for (; ((j < n) && (k < n2)); k++, j++) {
			p = inPtr + (ptrdiff_t)(n2 - k);
			q = kernel;
			sum = 0.0;
			for (i = k, l = ((nk + k) < (n - 1L)) ? (nk + k) : (n - 1L); (i < l); i++)
			  sum += *p++ * *q--;
			for (l = (n2 < (nk + k)) ? (n2) : (nk + k); (i < l); i++)
			  sum += *p-- * *q--;
			for (l = nk + k; (i < l); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		  for (k = 0L; (k < (n - j)); k++) {
			p = inPtr + (ptrdiff_t)k;
			q = kernel;
			sum = 0.0;
			for (i = k, l = ((nk + k) < (n - 1L)) ? (nk + k) : (n - 1L); (i < l); i++)
			  sum += *p++ * *q--;
			for (l = (n2 < (nk + k)) ? (n2) : (nk + k); (i < l); i++)
			  sum += *p-- * *q--;
			for (l = nk + k; (i < l); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		}
		else
		  for (j = -Ok - shift; (j < (n - Ok - shift)); j++) {
			k = (n == 1L) ? (0L) : ((j < 0L) ? (j + n2 * ((n2 - 1L - j) / n2))
			  : (j - n2 * (j / n2)));
			p = (k < n) ? (inPtr + (ptrdiff_t)k) : (inPtr + (ptrdiff_t)(n2 - k));
			q = kernel;
			sum = 0.0;
			for (i = -nk; (i < 0L); i++) {
			  sum += *p * *q--;
			  if (++k < n)
				p++;
			  else
				if (--p <= inPtr) {
				  p = inPtr;
				  k = 0L;
				}
			}
			*outPtr++ = sum;
		  }
} /* End of firConvolveMirror */

/************************************************************************/
/* +===============================+									*/
/* | FUNCTION: firConvolvePeriodic |									*/
/* +===============================+									*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete convolution									*/
/*				outPtr[x] = SUM(k): inPtr[x] * kernel[k - x]			*/
/*																		*/
/* Conventions:	The kernel has a finite support of size nk				*/
/*				The signal has an infinite, periodic (n) support		*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*				The kernel origin (hot spot) is at index [(nk - 1) / 2]	*/
/*				The kernel has to be given straight (k[x], *not* k[-x])	*/
/*																		*/
/************************************************************************/
void				firConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk) {

		double				*p, *q, *in;
		double				sum;
		long				Ok;
		long				i, j, k;

		kernel += (ptrdiff_t)(nk - 1L);
		shift = ((nk & 1L) != 0L) ? (shift) : (shift + 1L);
		shift = (shift < 0L) ? (shift + n * ((n - 1L - shift) / n)) : ((shift >= n) ? (shift - n
		  * (shift / n)) : (shift));
		outPtr += (ptrdiff_t)shift;
		Ok = (nk - 1L) / 2L;
		if (nk <= n) {
		  for (j = -Ok, k = n - Ok; (j < 0L); k++, j++) {
			if (shift++ == n) {
			  shift = 1L;
			  outPtr -= (ptrdiff_t)n;
			}
			p = inPtr + (ptrdiff_t)k;
			q = kernel;
			sum = 0.0;
			for (i = j; (i < 0L); i++)
			  sum += *p++ * *q--;
			p = inPtr;
			for (i = -j - nk; (i < 0L); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		  in = inPtr;
		  for (j = nk - n - 1L; (j < 0L); j++) {
			if (shift++ == n) {
			  shift = 1L;
			  outPtr -= (ptrdiff_t)n;
			}
			p = in++;
			q = kernel;
			sum = 0.0;
			for (i = -nk; (i < 0L); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		  for (j = Ok + 1L - nk, k = n + 1L - nk; (j < 0L); k++, j++) {
			if (shift++ == n) {
			  shift = 1L;
			  outPtr -= (ptrdiff_t)n;
			}
			p = inPtr + (ptrdiff_t)k;
			q = kernel;
			sum = 0.0;
			for (i = j - Ok; (i < 0L); i++)
			  sum += *p++ * *q--;
			p = inPtr;
			for (i = Ok - j - nk; (i < 0L); i++)
			  sum += *p++ * *q--;
			*outPtr++ = sum;
		  }
		}
		else
		  for (j = 0L; (j < n); j++) {
			if (shift++ == n) {
			  shift = 1L;
			  outPtr -= (ptrdiff_t)n;
			}
			k = j - Ok;
			k = (k < 0L) ? (k + n * ((n - 1L - k) / n)) : ((k >= n) ? (k - n * (k / n)) : (k));
			p = inPtr + (ptrdiff_t)k;
			q = kernel;
			sum = 0.0;
			for (i = -nk; (i < 0L); k++, i++) {
			  if (k == n) {
				k = 0L;
				p = inPtr;
			  }
			  sum += *p++ * *q--;
			}
			*outPtr++ = sum;
		  }
} /* End of firConvolvePeriodic */

/************************************************************************/
/* +=============================+										*/
/* | FUNCTION: iirConvolveFinite |										*/
/* +=============================+										*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete recursive convolution							*/
/*																		*/
/* Conventions:	Both the signal an the kernel have a finite support		*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*				The length of the kernel is np							*/
/*				The kernel is given by its (real) poles					*/
/*																		*/
/* Comment:		The output may share data storage with the input		*/
/*																		*/
/************************************************************************/
void				iirConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np) {

		double				*p, *q;
		double				pole;
		long				i, j;

		p = inPtr;
		q = outPtr;
		for (i = 0L; (i < n); i++)
		  *q++ = *p++ * gain;
		if (n == 1L)
		  return;
		for (i = -np; (i < 0L); i++) {
		  pole = *realPoles++;
		  p = q = outPtr;
		  for (p++, j = 1L - n; (j < 0L); j++)
			*p++ += *q++ * pole;
		  p--;
		  while (--q >= outPtr)
			*q += *p-- * pole;
		}
} /* End of iirConvolveFinite */

/************************************************************************/
/* +=============================+										*/
/* | FUNCTION: iirConvolveMirror |										*/
/* +=============================+										*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete recursive convolution							*/
/*																		*/
/* Conventions:	The kernel has a size np								*/
/*				The kernel is given by its (real) poles					*/
/*				The signal has an infinite, periodic (2 n - 1) support	*/
/*				The input is described by n values; the missing values	*/
/*					are provided by mirror symetry						*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*																		*/
/* Comment:		The output may share data storage with the input		*/
/*																		*/
/************************************************************************/
void				iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np) {

		double				*p, *q, *r;
		double				tolerance = log10((double)FLT_EPSILON);
		double				x0, pole;
		long				n2;
		long				i, j, k;

		p = inPtr;
		q = outPtr;
		r = outPtr + (ptrdiff_t)n;
		while (q < r)
		  *q++ = *p++ * gain;
		if (n == 1L)
		  return;
		n2 = 2L * (n - 1L);
		for (i = -np; (i < 0L); i++) {
		  pole = *realPoles++;
		  j = (long)ceil(tolerance / log10(fabs(pole)));
		  k = j - n2 * (j / n2);
		  j -= k;
		  if (k < n) {
			p = outPtr + (ptrdiff_t)k;
			x0 = *p;
		  }
		  else {
			k = n2 - k;
			p = outPtr + (ptrdiff_t)k;
			x0 = *p;
			while (++p < r)
			  x0 = pole * x0 + *p;
			p--;
		  }
		  while (--p >= outPtr)
			x0 = pole * x0 + *p;
		  while (j > 0L) {
			p++;
			while (++p < r)
			  x0 = pole * x0 + *p;
			p--;
			while (--p >= outPtr)
			  x0 = pole * x0 + *p;
			j -= n2;
		  }
		  q = p++;
		  *p++ = x0;
		  x0 = *(q++ + (ptrdiff_t)n);
		  while (p < r)
			*p++ += *q++ * pole;
		  *q = (2.0 * *q - x0) / (1.0 - pole * pole);
		  p--;
		  while (--q >= outPtr)
			*q += *p-- * pole;
		}
} /* End of iirConvolveMirror */

/************************************************************************/
/* +===============================+									*/
/* | FUNCTION: iirConvolvePeriodic |									*/
/* +===============================+									*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Discrete recursive convolution							*/
/*																		*/
/* Conventions:	The kernel has a size np								*/
/*				The kernel is given by its (real) poles					*/
/*				The signal has an infinite, periodic (n) support		*/
/*				The output has same length n than the input				*/
/*				The output has already been allocated					*/
/*																		*/
/* Comment:		The output may share data storage with the input		*/
/*																		*/
/************************************************************************/
void				iirConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np) {

		double				*p, *q;
		double				tolerance = log10((double)FLT_EPSILON);
		double				x0, x1, pole;
		long				i, j, k, k0;

		p = inPtr;
		q = outPtr;
		for (i = 0L; (i < n); i++)
		  *q++ = *p++ * gain;
		if (n == 1L)
		  return;
		for (i = -np; (i < 0L); i++) {
		  pole = *realPoles++;
		  k = (long)ceil(tolerance / log10(fabs(pole)));
		  k0 = k - n * (k / n);
		  p = outPtr + (ptrdiff_t)(n - 1L - k0);
		  q = outPtr + (ptrdiff_t)k0;
		  x0 = *(p++);
		  x1 = *(q--);
		  for (j = -k0; (j < 0L); j++) {
			x0 = pole * x0 + *p++;
			x1 = pole * x1 + *q--;
		  }
		  for (k -= k0; (k > 0L); k -= n) {
			p = outPtr;
			q = outPtr + (ptrdiff_t)(n - 1L);
			for (j = -n; (j < 0L); j++) {
			  x0 = pole * x0 + *p++;
			  x1 = pole * x1 + *q--;
			}
		  }
		  *++q += pole * x0;
		  x0 = *--p;
		  x1 = pole * x1 + x0;
		  p = q++;
		  for (j = 1L - n; (j < 0L); j++)
			*q++ += *p++ * pole;
		  *p = (*p + x1 - x0) / (1.0 - pole * pole);
		  q--;
		  while (--p >= outPtr)
			*p += *q-- * pole;
		}
} /* End of iirConvolvePeriodic */
