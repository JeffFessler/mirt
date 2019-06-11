#include		<stdlib.h>
#include		<stdio.h>
#include		<stddef.h>
#include		<string.h>
#include		<float.h>
#include		<math.h>
#include		<limits.h>
#include		"phil.h"

#include		"BsplnWgt.h"
#include		"quant.h"

/************************************************************************/
static	int			buildFloatImage	(short				*inPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz,
									 float				*floatOutPtr,
									 double				*histogramPtr,
									 short				minLevel);
static	int			buildShortImage	(short				*inPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz,
									 short				*shortOutPtr,
									 double				*histogramPtr,
									 short				minLevel,
									 unsigned			inLevels);
static	int			computeHistSnr	(short				*inPtr,
									 double				*histogramPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz);
static	int			crossParzenHistogram
									(float				*inPtr1,
									 float				*inPtr2,
									 double				*histPtr,
									 double				*mse,
									 double				*snr,
									 int				parzenDegree,
									 double				scale1,
									 double				scale2,
									 float				min1,
									 float				min2,
									 unsigned			bins,
									 int				nx,
									 int				ny,
									 int				nz);
static	void		deleteList		(struct	histList	*head);
static	int			doConvergence	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels,
									 double				*kernelsPtr,
									 double				epsilon);
static	int			equalizedKern	(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 double				*kernelsPtr);
static	int			equalizeHist	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			inLevels,
									 unsigned			outLevels);
static	int			findKernels		(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 double				*kernelsPtr,
									 double				epsilon);
static	int			getFloatHistogram
									(double				**histogramHdl,
									 float				**levelsHdl,
									 long				*levels,
									 float				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz);
static	int			getHistogram	(double				**histogramHdl,
									 unsigned			*inLevelsPtr,
									 short				*minLevel,
									 short				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz);
static	int			insert			(struct	histList	*head,
									 float				value);
static	void		kernToThresh	(double				*kernelsPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels);
static	int			parzenHistogram	(float				*inPtr,
									 double				*histPtr,
									 int				parzenDegree,
									 double				scale,
									 float				min,
									 unsigned			bins,
									 int				nx,
									 int				ny,
									 int				nz);
static	void		rampHistogram	(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 short				oMin,
									 short				oMax);
static	void		setLookUpTable	(double				*histogramPtr,
									 double				*kernelsPtr,
									 short				minLevel,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 int				labeling,
									 short				step);
static	void		sliceHistogram	(double				*histogramPtr,
									 unsigned			inLevels,
									 short				minLevel,
									 short				iMin,
									 short				iMax,
									 float				oMin,
									 float				oMax,
									 float				loBackgrnd,
									 float				hiBackgrnd);
static	void		threshToKern	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels,
									 double				*kernelsPtr);

/************************************************************************/
static	int			buildFloatImage	(short				*inPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz,
									 float				*floatOutPtr,
									 double				*histogramPtr,
									 short				minLevel) {

		int					i, j, k;

		if (computeHistSnr(inPtr, histogramPtr, mse, snr, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to compute signal to noise ratio");
		  return(ERROR);
		}
		for (k = 0; (k < nz); k++)
		  for (j = 0; (j < ny); j++)
			for (i = 0; (i < nx); inPtr++, floatOutPtr++, i++)
			  *floatOutPtr = (float)*(histogramPtr + (ptrdiff_t)*inPtr - (ptrdiff_t)minLevel);
		return(!ERROR);
} /* End of buildFloatImage */

/************************************************************************/
static	int			buildShortImage	(short				*inPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz,
									 short				*shortOutPtr,
									 double				*histogramPtr,
									 short				minLevel,
									 unsigned			inLevels) {

		double				*hstPtr;
		double				r, rounded;
		int					i, j, k;
		unsigned			n;

		hstPtr = histogramPtr;
		for (n = 0U; (n < inLevels); hstPtr++, n++) {
		  r = *hstPtr + 0.5;
		  rounded = (double)((long)r);
		  if ((r < 0.0) && (rounded != r))
			*hstPtr = rounded - 1.0;
		  else
			*hstPtr = rounded;
		}
		if (computeHistSnr(inPtr, histogramPtr, mse, snr, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to compute signal to noise ratio");
		  return(ERROR);
		}
		for (k = 0; (k < nz); k++)
		  for (j = 0; (j < ny); j++)
			for (i = 0; (i < nx); inPtr++, shortOutPtr++, i++)
			  *shortOutPtr = (short)*(histogramPtr + (ptrdiff_t)*inPtr - (ptrdiff_t)minLevel);
		return(!ERROR);
} /* End of buildShortImage */

/************************************************************************/
static	int			computeHistSnr	(short				*inPtr,
									 double				*histogramPtr,
									 double				*mse,
									 double				*snr,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				*oldHistogramPtr, *oHstPtr;
		double				sum;
		double				sumIn;
		short				minLevel;
		short				i;
		unsigned			inLevels, n;

		if (getHistogram(&oldHistogramPtr, &inLevels, &minLevel, inPtr, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to get histogram");
		  return(ERROR);
		}
		sum = sumIn = *mse = 0.0;
		i = minLevel;
		oHstPtr = oldHistogramPtr;
		for (n = 0U; (n < inLevels); oHstPtr++, histogramPtr++, n++, i++) {
		  sumIn += *oHstPtr * (double)i * (double)i;
		  *mse += *oHstPtr * ((double)i - *histogramPtr) * ((double)i - *histogramPtr);
		  sum += *oHstPtr;
		}
		free(oldHistogramPtr);
		if ((*mse == 0.0) || (sumIn == 0.0))
		  *snr = 0.0;
		else
		  *snr = 10.0 * log10(sumIn / *mse);
		if (sum == 0.0)
		  *mse = 0.0;
		else
		  *mse /= sum;
		return(!ERROR);
} /* End of computeHistSnr */

/************************************************************************/
static	int			crossParzenHistogram
									(float				*inPtr1,
									 float				*inPtr2,
									 double				*histPtr,
									 double				*mse,
									 double				*snr,
									 int				parzenDegree,
									 double				scale1,
									 double				scale2,
									 float				min1,
									 float				min2,
									 unsigned			bins,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				*crossHistogram;
		double				*p, *q;
		double				center1, center2;
		double				h, w, scaledW;
		double				r;
		long				rounded;
		long				lBins;
		long				kFrom1, k1, kTo1;
		long				kFrom2, k2, kTo2;
		long				i, j, k;
		int					x, y, z;

		scaledW = (double)(parzenDegree + 1) * 0.5;
		k = (long)(scaledW + 0.25);
		bins = bins + 2U * (unsigned)k;
		lBins = (long)bins;
		crossHistogram = (double *)malloc((size_t)bins * (size_t)bins * sizeof(double));
		if (crossHistogram == (double *)NULL) {
		  message("ERROR - Not enough memory for widening crossHistogram");
		  return(ERROR);
		}
		p = crossHistogram;
		for (i = 0L; (i < lBins); i++)
		  for (j = 0L; (j < lBins); p++, j++)
			*p = 0.0;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr1++, inPtr2++, x++) {
			  *mse += (double)((*inPtr2 - *inPtr1) * (*inPtr2 - *inPtr1));
			  *snr += (double)(*inPtr2 * *inPtr2);
			  center1 = (double)(*inPtr1 - min1) * scale1 + scaledW;
			  r = center1 - scaledW;
			  rounded = (long)r;
			  if ((r > 0.0) && ((double)rounded != r))
				kFrom1 = rounded + 1L;
			  else
				kFrom1 = rounded;
			  kFrom1 = (0L <= kFrom1) ? (kFrom1) : (0L);
			  r = center1 + scaledW;
			  rounded = (long)r;
			  if ((r < 0.0) && ((double)rounded != r))
				kTo1 = rounded - 1L;
			  else
				kTo1 = rounded;
			  kTo1 = (kTo1 < lBins) ? (kTo1) : (lBins - 1L);
			  center2 = (double)(*inPtr2 - min2) * scale2 + scaledW;
			  r = center2 - scaledW;
			  rounded = (long)r;
			  if ((r > 0.0) && ((double)rounded != r))
				kFrom2 = rounded + 1L;
			  else
				kFrom2 = rounded;
			  kFrom2 = (0L <= kFrom2) ? (kFrom2) : (0L);
			  r = center2 + scaledW;
			  rounded = (long)r;
			  if ((r < 0.0) && ((double)rounded != r))
				kTo2 = rounded - 1L;
			  else
				kTo2 = rounded;
			  kTo2 = (kTo2 < lBins) ? (kTo2) : (lBins - 1L);
			  p = crossHistogram + (ptrdiff_t)(kFrom2 * lBins);
			  k2 = kFrom2;
			  while (k2 <= kTo2) {
				k1 = kFrom1;
				q = p + (ptrdiff_t)kFrom1;
				w = BsplnWght(parzenDegree, 0L, center2 - (double)k2++);
				while (k1 <= kTo1) {
				  h = BsplnWght(parzenDegree, 0L, center1 - (double)k1++);
				  *(q++) += h * w;
				}
				p += (ptrdiff_t)lBins;
			  }
			}
		h = 0.0;
		p = crossHistogram;
		for (i = 0L; (i < lBins); i++)
		  for (j = 0L; (j < lBins); p++, j++)
			h += *p;
		if ((h * h) < (double)FLT_EPSILON) {
		  free(crossHistogram);
		  message("ERROR - Empty cross-histogram");
		  return(ERROR);
		}
		p = crossHistogram;
		for (i = 0L; (i < lBins); i++)
		  for (j = 0L; (j < lBins); p++, j++)
			*p /= h;
		p = crossHistogram;
		for (i = 0L; (i < lBins); i++) {
		  h = 0.0;
		  for (j = 0L; (j <= k); j++) {
			h += *p;
			*p++ = 0.0;
		  }
		  *--p = h;
		  p += (ptrdiff_t)(lBins - 2L * k - 1L);
		  h = 0.0;
		  for (j = lBins - (k + 1L); (j < lBins); j++) {
			h += *p;
			*p++ = 0.0;
		  }
		  *(p - (ptrdiff_t)(k + 1L)) = h;
		}
		p = crossHistogram + (ptrdiff_t)k;
		for (j = k; (j < (lBins - k)); j++) {
		  h = 0.0;
		  for (i = 0L; (i <= k); i++) {
			h += *p;
			*p = 0.0;
			p += (ptrdiff_t)lBins;
		  }
		  p -= (ptrdiff_t)lBins;
		  *p = h;
		  p += (ptrdiff_t)(lBins * (lBins - 2L * k - 1L));
		  h = 0.0;
		  for (i = lBins - (k + 1L); (i < lBins); i++) {
			h += *p;
			*p = 0.0;
			p += (ptrdiff_t)lBins;
		  }
		  *(p - (ptrdiff_t)(lBins * (k + 1L))) = h;
		  p -= (ptrdiff_t)(lBins * lBins - 1L);
		}
		p = crossHistogram + (ptrdiff_t)(k * (lBins + 1L));
		q = histPtr;
		bins = bins - 2U * (unsigned)k;
		for (i = 0L; (i < (long)bins); i++) {
		  q = (double *)memcpy(q, p, (size_t)bins * sizeof(double));
		  p += (ptrdiff_t)lBins;
		  q += (ptrdiff_t)bins;
		}
		if ((*mse * *snr) != 0.0)
		  *snr = 10.0 * log10(*snr / *mse);
		else
		  *snr = 0.0;
		*mse /= (double)nx * (double)ny * (double)nz;
		free(crossHistogram);
		return(!ERROR);
} /* End of crossParzenHistogram */

/************************************************************************/
static	void		deleteList		(struct	histList	*head) {

		struct	histList	*p, *q;

		p = head->next;
		while (p != (struct histList *)NULL) {
		  q = p;
		  p = p->next;
		  free(q);
		}
} /* End of deleteList */

/************************************************************************/
static	int			doConvergence	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels,
									 double				*kernelsPtr,
									 double				epsilon) {

		double				*oldThresholdsPtr;
		double				*thrPtr, *oldPtr;
		double				s, sum;
		unsigned			i;

		oldPtr = oldThresholdsPtr = (double *)malloc((size_t)outLevels * sizeof(double));
		if (oldPtr == (double *)NULL) {
		  message("ERROR - Not enough memory for convergence criterion");
		  return(ERROR);
		}
		thrPtr = thresholdsPtr;
		for (i = 0U; (i < outLevels); oldPtr++, thrPtr++, i++)
		  *oldPtr = *thrPtr;
		do {
		  threshToKern(histogramPtr, thresholdsPtr, outLevels, kernelsPtr);
		  kernToThresh(kernelsPtr, thresholdsPtr, outLevels);
		  thrPtr = thresholdsPtr;
		  oldPtr = oldThresholdsPtr;
		  for (sum = s = 0.0, i = 0U; (i < outLevels); thrPtr++, oldPtr++, i++) {
			s += (*thrPtr - *oldPtr) * (*thrPtr - *oldPtr);
			sum += (*oldPtr) * (*oldPtr);
			*oldPtr = *thrPtr;
		  }
		} while ((s / sum) > epsilon);
		free(oldThresholdsPtr);
		return(!ERROR);
} /* End of doConvergence */

/************************************************************************/
static	int			equalizedKern	(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 double				*kernelsPtr) {

		double				*thresholdsPtr;

		thresholdsPtr = (double *)malloc((size_t)outLevels * sizeof(double));
		if (thresholdsPtr == (double *)NULL) {
		  message("ERROR - Not enough memory for holding thresholds");
		  return(ERROR);
		}
		if (equalizeHist(histogramPtr, thresholdsPtr, inLevels, outLevels) == ERROR) {
		  free(thresholdsPtr);
		  message("ERROR - Unable to equalize histogram");
		  return(ERROR);
		}
		threshToKern(histogramPtr, thresholdsPtr, outLevels, kernelsPtr);
		free(thresholdsPtr);
		return(!ERROR);
} /* End of equalizedKern */

/************************************************************************/
static	int			equalizeHist	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			inLevels,
									 unsigned			outLevels) {

		double				*accumulatedHistPtr;
		double				*hstPtr, *aHstPtr;
		double				*thrPtr;
		double				s, sum;
		unsigned			i, j;

		aHstPtr = accumulatedHistPtr = (double *)malloc((size_t)inLevels * sizeof(double));
		if (aHstPtr == (double *)NULL) {
		  message("ERROR - Not enough memory for holding image repartition");
		  return(ERROR);
		}
		hstPtr = histogramPtr;
		thrPtr = thresholdsPtr;
		*aHstPtr = *hstPtr;
		hstPtr++;
		aHstPtr++;
		for (i = 1U; (i < inLevels); hstPtr++, aHstPtr++, i++)
		  *aHstPtr = *(aHstPtr - (ptrdiff_t)1) + *hstPtr;
		aHstPtr--;
		sum = *aHstPtr;
		if (sum * sum == 0.0) {
		  free(accumulatedHistPtr);
		  message("ERROR - Histogram is empty");
		  return(ERROR);
		}
		for (i = 1U; (i <= outLevels); thrPtr++, i++)
		  *thrPtr = sum * (double)i / (double)outLevels;
		hstPtr = histogramPtr;
		aHstPtr = accumulatedHistPtr;
		thrPtr = thresholdsPtr;
		i = 0U;
		j = 1U;
		s = 0.0;
		while (j < outLevels) {
		  while ((s <= *thrPtr) && (*thrPtr < *aHstPtr)) {
			 *thrPtr -= s;
			 *thrPtr /= *hstPtr;
			 *thrPtr -= 0.5;
			 *thrPtr += (double)i;
			j++;
			thrPtr++;
		  }
		  s = *aHstPtr;
		  i++;
		  hstPtr++;
		  aHstPtr++;
		}
		*thrPtr = (double)inLevels - 0.5;
		free(accumulatedHistPtr);
		return(!ERROR);
} /* End of equalizeHist */

/************************************************************************/
static	int			findKernels		(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 double				*kernelsPtr,
									 double				epsilon) {

		double				*thresholdsPtr;

		thresholdsPtr = (double *)malloc((size_t)outLevels * sizeof(double));
		if (thresholdsPtr == (double *)NULL) {
		  message("ERROR - Not enough memory for holding thresholds");
		  return(ERROR);
		}
		if (equalizeHist(histogramPtr, thresholdsPtr, inLevels, outLevels) == ERROR) {
		  free(thresholdsPtr);
		  message("ERROR - Unable to equalize histogram");
		  return(ERROR);
		}
		if (doConvergence(histogramPtr, thresholdsPtr, outLevels, kernelsPtr, epsilon) == ERROR) {
		  free(thresholdsPtr);
		  message("ERROR - Unable to converge");
		  return(ERROR);
		}
		free(thresholdsPtr);
		return(!ERROR);
} /* End of findKernels */

/************************************************************************/
static	int			getFloatHistogram
									(double				**histogramHdl,
									 float				**levelsHdl,
									 long				*levels,
									 float				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz) {

		struct	histList	head;
		struct	histList	*p, *q;
		double				*h;
		float				*l;
		double				sum;
		int					x, y, z;

		p = (struct histList *)malloc(sizeof(struct histList));
		if (p == (struct histList *)NULL) {
		  message("ERROR - Not enough memory for holding histogram tail");
		  return(ERROR);
		}
		p->next = (struct histList *)NULL;
		p->occurences = 0.0;
		p->value = *inPtr;
		head.next = p;
		head.occurences = 0.0;
		head.value = 0.0F;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, x++)
			  if (insert(&head, *inPtr) == ERROR) {
				deleteList(&head);
				message("ERROR - Unable to get float histogram");
				return(ERROR);
			  }
		*levels = 0L;
		p = head.next;
		sum = 0.0;
		while (p != (struct histList *)NULL) {
		  (*levels)++;
		  sum += p->occurences;
		  p = p->next;
		}
		*histogramHdl = (double *)malloc((size_t)*levels * sizeof(double));
		if (*histogramHdl == (double *)NULL) {
		  deleteList(&head);
		  message("ERROR - Not enough memory for holding histogram");
		  return(ERROR);
		}
		*levelsHdl = (float *)malloc((size_t)*levels * sizeof(float));
		if (*levelsHdl == (float *)NULL) {
		  free(*histogramHdl);
		  deleteList(&head);
		  message("ERROR - Not enough memory for holding levels");
		  return(ERROR);
		}
		h = *histogramHdl;
		l = *levelsHdl;
		p = head.next;
		while (p != (struct histList *)NULL) {
		  *(h++) = p->occurences / sum;
		  *(l++) = p->value;
		  q = p;
		  p = p->next;
		  free(q);
		}
		return(!ERROR);
} /* End of getFloatHistogram */

/************************************************************************/
static	int			getHistogram	(double				**histogramHdl,
									 unsigned			*inLevelsPtr,
									 short				*minLevelPtr,
									 short				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				*hstPtr;
		short				*imgPtr;
		int					i, j, k;
		short				minLevel, maxLevel;
		unsigned			n;

		imgPtr = inPtr;
		minLevel = maxLevel = *imgPtr;
		for (k = 0; (k < nz); k++)
		  for (j = 0; (j < ny); j++)
			for (i = 0; (i < nx); imgPtr++, i++)
			  if (maxLevel < *imgPtr)
				maxLevel = *imgPtr;
			  else if (*imgPtr < minLevel)
				minLevel = *imgPtr;
		*minLevelPtr = minLevel;
		if (((long)maxLevel - (long)minLevel) >= (long)USHRT_MAX) {
		  message("ERROR - Too many levels");
		  return(ERROR);
		}
		*inLevelsPtr = (unsigned)((long)maxLevel - (long)minLevel + 1L);
		hstPtr = *histogramHdl = (double *)malloc((size_t)*inLevelsPtr * sizeof(double));
		if (hstPtr == (double *)NULL) {
		  message("ERROR - Not enough memory for holding histogram");
		  return(ERROR);
		}
		for (n = 0U; (n < *inLevelsPtr); hstPtr++, n++)
		  *hstPtr = 0.0;
		imgPtr = inPtr;
		for (k = 0; (k < nz); k++)
		  for (j = 0; (j < ny); j++)
			for (i = 0; (i < nx); imgPtr++, i++)
			  *(*histogramHdl + (ptrdiff_t)*imgPtr - (ptrdiff_t)minLevel) += 1.0;
		return(!ERROR);
} /* End of getHistogram */

/************************************************************************/
static	int			insert			(struct	histList	*head,
									 float				value) {

		struct	histList	*p, *q, *r;

		r = head;
		p = head->next;
		while (p != (struct histList *)NULL) {
		  if (p->value > value) {
			q = (struct histList *)malloc(sizeof(struct histList));
			if (q == (struct histList *)NULL) {
			  message("ERROR - Not enough memory for holding histogram entry");
			  return(ERROR);
			}
			q->next = p;
			q->occurences = 1.0;
			q->value = value;
			r->next = q;
			return(!ERROR);
		  }
		  if (p->value == value) {
			p->occurences += 1.0;
			return(!ERROR);
		  }
		  r = p;
		  p = p->next;
		}
		q = (struct histList *)malloc(sizeof(struct histList));
		if (q == (struct histList *)NULL) {
		  message("ERROR - Not enough memory for holding histogram entry");
		  return(ERROR);
		}
		q->next = (struct histList *)NULL;
		q->occurences = 1.0;
		q->value = value;
		r->next = q;
		return(!ERROR);
} /* End of insert */

/************************************************************************/
static	void		kernToThresh	(double				*kernelsPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels) {

		double				s;
		unsigned			i;

		s = *kernelsPtr;
		kernelsPtr++;
		for(i = 1U; (i < outLevels); kernelsPtr++, thresholdsPtr++, i++) {
		  *thresholdsPtr = (s + *kernelsPtr) / 2.0;
		  s = *kernelsPtr;
		}
} /* End of kernToThresh */

/************************************************************************/
static	int			parzenHistogram	(float				*inPtr,
									 double				*histPtr,
									 int				parzenDegree,
									 double				scale,
									 float				min,
									 unsigned			bins,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				*p;
		double				center;
		double				h, scaledW;
		double				r;
		long				rounded;
		long				lBins;
		long				kFrom, k, kTo;
		int					x, y, z;
		unsigned			n;

		p = histPtr;
		for (n = 0U; (n < bins); p++, n++)
		  *p = 0.0;
		lBins = (long)bins;
		scaledW = (double)(parzenDegree + 1) * 0.5;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, x++) {
			  center = (double)(*inPtr - min) * scale;
			  r = center - scaledW;
			  rounded = (long)r;
			  if ((r > 0.0) && ((double)rounded != r))
				kFrom = rounded + 1L;
			  else
				kFrom = rounded;
			  kFrom = (0L <= kFrom) ? (kFrom) : (0L);
			  r = center + scaledW;
			  rounded = (long)r;
			  if ((r < 0.0) && ((double)rounded != r))
				kTo = rounded - 1L;
			  else
				kTo = rounded;
			  kTo = (kTo < lBins) ? (kTo) : (lBins - 1L);
			  for (k = kFrom; (k <= kTo); k++)
				*(histPtr + (ptrdiff_t)k) += BsplnWght(parzenDegree, 0L, center - (double)k);
			}
		h = 0.0;
		p = histPtr;
		for (n = 0U; (n < bins); p++, n++)
		  h += *p;
		if ((h * h) < (double)FLT_EPSILON) {
		  message("ERROR - Empty histogram");
		  return(ERROR);
		}
		p = histPtr;
		for (n = 0U; (n < bins); p++, n++)
		  *p /= h;
		return(!ERROR);
} /* End of parzenHistogram */

/************************************************************************/
static	void		rampHistogram	(double				*histogramPtr,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 short				oMin,
									 short				oMax) {

		double				*hstPtr;
		double				r, rounded;
		unsigned			n;

		hstPtr = histogramPtr;
		for (n = 0U; (n < inLevels); hstPtr++, n++)
		  *hstPtr = (double)n * ((double)outLevels - 1.0) / ((double)inLevels - 1.0);
		hstPtr = histogramPtr;
		for (n = 0U; (n < inLevels); hstPtr++, n++) {
		  r = *hstPtr + 0.5;
		  rounded = (double)((long)r);
		  if ((r < 0.0) && (rounded != r))
			*hstPtr = (double)oMin + ((double)oMax - (double)oMin) * (rounded - 1.0)
			  / ((double)outLevels - 1.0);
		  else
			*hstPtr = (double)oMin + ((double)oMax - (double)oMin) * rounded
			  / ((double)outLevels - 1.0);
		}
} /* End of rampHistogram */

/************************************************************************/
static	void		setLookUpTable	(double				*histogramPtr,
									 double				*kernelsPtr,
									 short				minLevel,
									 unsigned			inLevels,
									 unsigned			outLevels,
									 int				labeling,
									 short				step) {

		unsigned			i, j;
		double				thr;

		i = 0U;
		for (j = 1U; (j < outLevels); kernelsPtr++, j++) {
		  thr = 0.5 * (*kernelsPtr + *(kernelsPtr + 1));
		  while ((double)i < thr) {
			if (labeling)
			  *histogramPtr = (double)step * ((double)j - 1.0) + (double)minLevel;
			else
			  *histogramPtr = *kernelsPtr + (double)minLevel;
			histogramPtr++;
			i++;
		  }
		}
		while (i < inLevels) {
		  if (labeling)
			*histogramPtr = (double)step * ((double)j - 1.0) + (double)minLevel;
		  else
			*histogramPtr = *kernelsPtr + (double)minLevel;
		  histogramPtr++;
		  i++;
		}
} /* End of setLookUpTable */

/************************************************************************/
static	void		sliceHistogram	(double				*histogramPtr,
									 unsigned			inLevels,
									 short				minLevel,
									 short				iMin,
									 short				iMax,
									 float				oMin,
									 float				oMax,
									 float				loBackgrnd,
									 float				hiBackgrnd) {

		short				i, maxLevel;
		unsigned			m, n, sliceLen;

		maxLevel = (short)((long)minLevel + (long)inLevels - 1L);
		sliceLen = (unsigned)((long)iMax - (long)iMin + 1L);
		if (maxLevel < iMin) {
		  for (m = 0U; (m < inLevels); histogramPtr++, m++)
			*histogramPtr = (double)loBackgrnd;
		  return;
		}
		if (iMax < minLevel) {
		  for (m = 0U; (m < inLevels); histogramPtr++, m++)
			*histogramPtr = (double)hiBackgrnd;
		  return;
		}
		if (iMin >= minLevel) {
		  for (m = 0U, i = minLevel; (i < iMin); histogramPtr++,
			m++, i++)
			*histogramPtr = (double)loBackgrnd;
		  if (iMax > maxLevel)
			sliceLen -= (unsigned)((long)iMax - (long)maxLevel);
		  for (n = 0U; (n < sliceLen); histogramPtr++, m++, n++)
			if (iMin == iMax)
			  *histogramPtr = (double)oMin;
			else
			  *histogramPtr = (double)oMin + ((double)oMax - (double)oMin) * (double)n
				/ ((double)iMax - (double)iMin);
		  while (m < inLevels) {
			*histogramPtr = (double)hiBackgrnd;
			histogramPtr++;
			m++;
		  }
		}
		else {
		  n = (unsigned)((long)minLevel - (long)iMin);
		  sliceLen -= n;
		  if (iMax > maxLevel)
			sliceLen -= (unsigned)((long)iMax - (long)maxLevel);
		  for (m = 0U; (m < sliceLen); histogramPtr++, n++, m++)
			if (iMin == iMax)
			  *histogramPtr = (double)oMin;
			else
			  *histogramPtr = (double)oMin + ((double)oMax - (double)oMin) * (double)n
				/ ((double)iMax - (double)iMin);
		  while (m < inLevels) {
			*histogramPtr = (double)hiBackgrnd;
			histogramPtr++;
			m++;
		  }
		}
} /* End of sliceHistogram */

/************************************************************************/
static	void		threshToKern	(double				*histogramPtr,
									 double				*thresholdsPtr,
									 unsigned			outLevels,
									 double				*kernelsPtr) {

		double				*thrPtr;
		double				*hstPtr;
		double				*kPtr;
		double				s, norm;
		unsigned			i, j;

		kPtr = kernelsPtr;
		hstPtr = histogramPtr;
		thrPtr = thresholdsPtr;
		for (i = j = 0U, s = -0.5; (j < outLevels); thrPtr++, kPtr++, j++) {
		  if (*thrPtr <= ((double)i + 0.5)) {
			*kPtr = *hstPtr * ((*thrPtr) * (*thrPtr) - s * s) / 2.0;
			norm = *hstPtr * (*thrPtr - s);
		  }
		  else {
			*kPtr = *hstPtr * (((double)i - s) * ((double)i + s) + (double)i + 0.25) / 2.0;
			norm = *hstPtr * ((double)i + 0.5 - s);
			hstPtr++;
			i++;
			while (*thrPtr > (double)i + 0.5) {
			  *kPtr += *hstPtr * (double)i;
			  norm += *hstPtr;
			  hstPtr++;
			  i++;
			}
			*kPtr += *hstPtr*((*thrPtr - (double)i) * (*thrPtr + (double)i) + (double)i - 0.25)
			  / 2.0;
			norm += *hstPtr * (*thrPtr - (double)i + 0.5);
		  }
		  if (norm == 0.0)
			*kPtr = 0.5 * (s + *thrPtr);
		  else
			*kPtr /= norm;
		  s = *thrPtr;
		}
} /* End of threshToKern */

/************************************************************************/
/* FUNCTION: histEntropy												*/
/************************************************************************/
int					histEntropy		(struct qParam		*entDataPtr) {

		double				*histogramPtr1, *histogramPtr2;
		double				*crossHistogram;
		double				*h1, *h2, *ch;
		short				*inPtr1, *inPtr2;
		struct qParam		quantData;
		double				H, mse, snr;
		double				nxyz;
		int					nx, ny, nz;
		int					x, y, z;
		short				minLevel1, minLevel2;
		unsigned			levels1, levels2;
		unsigned			i, j;

		if (entDataPtr->floatInput) {
		  message("ERROR - Input has to be of type (short)");
		  return(ERROR);
		}
		inPtr1 = entDataPtr->inPtr.s;
		inPtr2 = entDataPtr->outPtr.s;
		nx = entDataPtr->nx;
		ny = entDataPtr->ny;
		nz = entDataPtr->nz;
		nxyz = (double)nx * (double)ny * (double)nz;
		quantData = *entDataPtr;
		quantData.inPtr.s = inPtr1;
		quantData.outPtr.s = inPtr1;
		quantData.floatOutput = FALSE;
		quantData.first = (short)SHRT_MIN;
		quantData.step = (short)1;
		if (histQuant(&quantData) == ERROR) {
		  message("ERROR - Unable to adjust the number of levels for 1st volume");
		  return(ERROR);
		}
		if (getHistogram(&histogramPtr1, &levels1, &minLevel1, inPtr1, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to get 1st histogram");
		  return(ERROR);
		}
		quantData.inPtr.s = inPtr2;
		quantData.outPtr.s = inPtr2;
		if (histQuant(&quantData) == ERROR) {
		  free(histogramPtr1);
		  message("ERROR - Unable to adjust the number of levels for 2nd volume");
		  return(ERROR);
		}
		if (getHistogram(&histogramPtr2, &levels2, &minLevel2, inPtr2, nx, ny, nz) == ERROR) {
		  free(histogramPtr1);
		  message("ERROR - Unable to get 2nd histogram");
		  return(ERROR);
		}
		crossHistogram = (double *)malloc((size_t)levels1 * (size_t)levels2 * sizeof(double));
		if (crossHistogram == (double *)NULL) {
		  free(histogramPtr2);
		  free(histogramPtr1);
		  message("ERROR - Not enough memory for holding crossHistogram");
		  return(ERROR);
		}
		ch = crossHistogram;
		for (j = 0U; (j < levels2); j++)
		  for (i = 0U; (i < levels1); ch++, i++)
			*ch = 0.0;
		mse = 0.0;
		snr = 0.0;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr1++, inPtr2++, x++) {
			  *(crossHistogram + (ptrdiff_t)(((long)*inPtr2 - (long)minLevel2) * (long)levels1
				+ (long)*inPtr1 - (long)minLevel1)) += 1.0;
			  mse += ((double)*inPtr2 - (double)*inPtr1) * ((double)*inPtr2 - (double)*inPtr1);
			  snr += (double)*inPtr2 * (double)*inPtr2;
			}
		if ((mse * snr) != 0.0)
		  snr = 10.0 * log10(snr / mse);
		else
		  snr = 0.0;
		mse /= nxyz;
		H = 0.0;
		h2 = histogramPtr2;
		ch = crossHistogram;
		for (j = 0U; (j < levels2); h2++, j++)
		  if (*h2 > 0.0) {
			h1 = histogramPtr1;
			for (i = 0U; (i < levels1); ch++, h1++, i++)
			  if ((*ch > 0.0) && (*h1 > 0.0))
				H += *ch * log(*ch * nxyz / (*h1 * *h2));
		  }
		  else
			ch += (ptrdiff_t)levels1;
		H /= log(2.0) * nxyz;
		entDataPtr->entropy = H;
		entDataPtr->mse = mse;
		entDataPtr->snr = snr;
		free(crossHistogram);
		free(histogramPtr2);
		free(histogramPtr1);
		return(!ERROR);
} /* End of histEntropy */

/************************************************************************/
/* FUNCTION: histogram													*/
/*----------------------------------------------------------------------*/
/*																		*/
/* By convention, we consider a value within a volume of type short as  */
/*   a representation of a bin of size 1, extending half to the left	*/
/*   and half to the right												*/
/*																		*/
/* -------|-------|-------|-------|-------|-------						*/
/*			  \_______/ <-- extended value								*/
/*																		*/
/* This means that the correspondence between a histogram of M levels	*/
/*   and one of N levels is done thus (M = 5, N = 3)					*/
/*																		*/
/* -------|-------|-------|-------|-------|-------						*/
/*	  \_______._______._______._______._______/							*/
/*																		*/
/* ----------|------------|------------|----------						*/
/*	  \____________._____________.____________/							*/
/*																		*/
/* In summary, we align to the extended value for type short			*/
/*																		*/
/*======================================================================*/
/*																		*/
/* By convention, we consider a value within a volume of type float as  */
/*   a representation of a bin of size 0. The correspondence between a	*/
/*   histogram of M levels and one of N levels becomes (M = 5, N = 3)	*/
/*																		*/
/* -------|-------|-------|-------|-------|-------						*/
/* -------|---------------|---------------|-------						*/
/*																		*/
/* In summary, we align to the exact value for type float				*/
/*																		*/
/************************************************************************/
int					histogram		(struct qParam		*histDataPtr) {

		FILE				*fp;
		double				*histogramPtr, *p;
		double				*outHistPtr;
		float				*floatInPtr, *f;
		float				*levelsPtr;
		short				*shortInPtr;
		double				scaledW;
		double				h, j, k;
		double				r, rounded;
		float				min, max;
		long				levels, l;
		int					x, y, z;
		short				minLevel;
		short				i;
		unsigned			inLevels, outLevels;
		unsigned			n;

		floatInPtr = histDataPtr->inPtr.f;
		shortInPtr = histDataPtr->inPtr.s;
		outLevels = histDataPtr->bins;
		fp = fopen(histDataPtr->fileName, "w");
		if (fp == (FILE *)NULL) {
		  message("ERROR - Unable to open histogram file");
		  return(ERROR);
		}
		if (fprintf(fp, "index intensity value\n") == EOF) {
		  if (fclose(fp) == EOF)
			message("ERROR - Unable to close histogram file");
		  message("ERROR - Unable to write histogram headings");
		  return(ERROR);
		}
		outHistPtr = (double *)malloc((size_t)outLevels * sizeof(double));
		if (outHistPtr == (double *)NULL) {
		  if (fclose(fp) == EOF)
			message("ERROR - Unable to close histogram file");
		  message("ERROR - Not enough memory for holding output histogram");
		  return(ERROR);
		}
		p = outHistPtr;
		for (n = 0U; (n < outLevels); p++, n++)
		  *p = 0.0;
		switch (histDataPtr->parzenDegree) {
		  case -1:
			message("WARNING - Ignoring the width and levels directives");
			if (histDataPtr->floatInput) {
			  free(outHistPtr);
			  if (getFloatHistogram(&outHistPtr, &levelsPtr, &levels, floatInPtr,
				histDataPtr->nx, histDataPtr->ny, histDataPtr->nz) == ERROR) {
				if (fclose(fp) == EOF)
				  message("ERROR - Unable to close histogram file");
				message("ERROR - Unable to get float histogram");
				return(ERROR);
			  }
			  p = outHistPtr;
			  f = levelsPtr;
			  for (l = 0L; (l < levels); p++, f++, l++)
				if (fprintf(fp, "%ld %e %e\n", l, (double)*f, *p) == EOF) {
				  free(outHistPtr);
				  free(levelsPtr);
				  if (fclose(fp) == EOF)
					message("ERROR - Unable to close histogram file");
				  message("ERROR - Unable to write histogram");
				  return(ERROR);
				}
			  free(levelsPtr);
			}
			else {
			  if (getHistogram(&histogramPtr, &inLevels, &minLevel, shortInPtr,
				histDataPtr->nx, histDataPtr->ny, histDataPtr->nz) == ERROR) {
				free(outHistPtr);
				if (fclose(fp) == EOF)
				  message("ERROR - Unable to close histogram file");
				message("ERROR - Unable to get histogram");
				return(ERROR);
			  }
			  h = 0.0;
			  p = histogramPtr;
			  for (n = 0U; (n < inLevels); p++, n++)
				h += *p;
			  i = minLevel;
			  p = histogramPtr;
			  for (n = 0U; (n < inLevels); p++, n++, i++)
				if (fprintf(fp, "%u %hd %e\n", n, i, *p / h) == EOF) {
				  free(histogramPtr);
				  free(outHistPtr);
				  if (fclose(fp) == EOF)
					message("ERROR - Unable to close histogram file");
				  message("ERROR - Unable to write histogram");
				  return(ERROR);
				}
			  free(histogramPtr);
			}
			break;
		  case 0:
		  case 1:
		  case 2:
		  case 3:
		  case 4:
		  case 5:
		  case 6:
		  case 7:
			if (histDataPtr->floatInput) {
			  f = floatInPtr;
			  min = *f;
			  max = *f;
			  for (z = 0; (z < histDataPtr->nz); z++)
				for (y = 0; (y < histDataPtr->ny); y++)
				  for (x = 0; (x < histDataPtr->nx); f++, x++) {
					min = (*f < min) ? (*f) : (min);
					max = (max < *f) ? (*f) : (max);
				  }
			  if (min == max) {
				h = 1.0 / (double)(outLevels);
				j = (double)min;
				for (n = 0U; (n < outLevels); n++)
				  if (fprintf(fp, "%u %e %e\n", n, j, h) == EOF) {
					free(outHistPtr);
					if (fclose(fp) == EOF)
					  message("ERROR - Unable to close histogram file");
					message("ERROR - Unable to write histogram");
					return(ERROR);
				  }
			  }
			  else {
				if (parzenHistogram(floatInPtr, outHistPtr, histDataPtr->parzenDegree,
				  (double)(outLevels - 1U) / ((double)max - (double)min), min, outLevels,
				  histDataPtr->nx, histDataPtr->ny, histDataPtr->nz) == ERROR) {
				  free(outHistPtr);
				  if (fclose(fp) == EOF)
					message("ERROR - Unable to close histogram file");
				  message("ERROR - Unable to perform parzenHistogram");
				  return(ERROR);
				}
				p = outHistPtr;
				j = (double)min;
				for (n = 0U; (n < outLevels); n++, p++, j += ((double)max - (double)min)
				  / (double)(outLevels - 1U))
				  if (fprintf(fp, "%u %e %e\n", n, j, *p) == EOF) {
					free(outHistPtr);
					if (fclose(fp) == EOF)
					  message("ERROR - Unable to close histogram file");
					message("ERROR - Unable to write histogram");
					return(ERROR);
				  }
			  }
			}
			else {
			  if (getHistogram(&histogramPtr, &inLevels, &minLevel, shortInPtr,
				histDataPtr->nx, histDataPtr->ny, histDataPtr->nz) == ERROR) {
				free(outHistPtr);
				if (fclose(fp) == EOF)
				  message("ERROR - Unable to close histogram file");
				message("ERROR - Unable to get histogram");
				return(ERROR);
			  }
			  if (inLevels < 2U) {
				h = 1.0 / (double)(outLevels);
				j = (double)minLevel;
				for (n = 0U; (n < outLevels); n++)
				  if (fprintf(fp, "%u %e %e\n", n, j, h) == EOF) {
					free(histogramPtr);
					free(outHistPtr);
					if (fclose(fp) == EOF)
					  message("ERROR - Unable to close histogram file");
					message("ERROR - Unable to write histogram");
					return(ERROR);
				  }
			  }
			  else {
				scaledW = (double)(histDataPtr->parzenDegree + 1) * 0.5;
				p = histogramPtr;
				for (n = 0U; (n < inLevels); p++, n++) {
				  j = (double)n * ((double)outLevels - 1.0) / ((double)inLevels - 1.0);
				  r = j - scaledW;
				  rounded = (double)((long)r);
				  if ((r > 0.0) && (rounded != r))
					k = rounded + 1.0;
				  else
					k = rounded;
				  k = (0L <= k) ? (k) : (0L);
				  r = j + scaledW;
				  rounded = (double)((long)r);
				  if ((r < 0.0) && (rounded != r))
					h = rounded - 1.0;
				  else
					h = rounded;
				  h = (h < (double)outLevels) ? (h) : ((double)outLevels - 1.0);
				  while (k < h) {
					*(outHistPtr + (ptrdiff_t)k) += *p * BsplnWght(histDataPtr->parzenDegree,
					  0L, j - k);
					k += 1.0;
				  }
				}
				h = 0.0;
				p = outHistPtr;
				for (n = 0U; (n < outLevels); p++, n++)
				  h += *p;
				p = outHistPtr;
				j = (double)minLevel - 0.5 + 0.5 * (double)inLevels / (double)outLevels;
				for (n = 0U; (n < outLevels); p++, n++, j += (double)inLevels
				  / (double)outLevels)
				  if (fprintf(fp, "%u %e %e\n", n, j, *p / h) == EOF) {
					free(histogramPtr);
					free(outHistPtr);
					if (fclose(fp) == EOF)
					  message("ERROR - Unable to close histogram file");
					message("ERROR - Unable to write histogram");
					return(ERROR);
				  }
			  }
			  free(histogramPtr);
			}
			break;
		  default:
			free(outHistPtr);
			if (fclose(fp) == EOF)
			  message("ERROR - Unable to close histogram file");
			message("ERROR - Unknown type of Parzen window");
			return(ERROR);
		}
		free(outHistPtr);
		if (fclose(fp) == EOF) {
		  message("ERROR - Unable to close histogram file");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of histogram */

/************************************************************************/
/* FUNCTION: histQuant													*/
/*----------------------------------------------------------------------*/
/* linear																*/
/*......................................................................*/
/*																		*/
/* Objective:	Cut histogram in equal pieces with respect to L levels	*/
/*																		*/
/*				x	 in [xMin, xMax]		is level					*/
/*				Q(x) in [0, L[				is quantization index k		*/
/*				R(k) in [xMin, xMax]		is reconstructed value		*/
/*																		*/
/*				k = Q(x) = round((L - 1) (x - xMin) / (xMax - xMin))	*/
/*				R(k) = xMin + (xMax - xMin) * k / (L - 1)				*/
/*																		*/
/*======================================================================*/
/* equalize																*/
/*......................................................................*/
/*																		*/
/* Objective:	Cut histogram in L equal density pieces					*/
/*																		*/
/*						 / x(k+1)										*/
/*				1 / L = |		 p(z) dz		for all k in [0, L]		*/
/*						/ x(k)											*/
/*																		*/
/*				x	is level											*/
/*				Q(x) is quantization index								*/
/*				R(k) is reconstructed value								*/
/*				p(x) is probability density function					*/
/*				x(k) are the (L + 1) thresholds							*/
/*																		*/
/* Solution:	Let y(x) be the repartition function of p(x)			*/
/*																		*/
/*							{ -inf				k == 0					*/
/*				then x(k) = { inv(y)(k / L)		k in ]0, L[				*/
/*							{ +inf				k == L					*/
/*																		*/
/*								 / x(k+1)								*/
/*				with R(k) = L * |		  z p(z) dz		k in [0, L[		*/
/*								/ x(k)									*/
/*																		*/
/* Model:		Piecewise constant p(x) given by histogram				*/
/*				Support for p(x) is [xMin - 0.5, xMax + 0.5 [			*/
/*																		*/
/*======================================================================*/
/* LloydMax																*/
/*......................................................................*/
/*																		*/
/* Objective:	Minimization of quadratic quantization error e			*/
/*																		*/
/*					 / inf					2							*/
/*				e = |		   (x - R(Q(x)))  p(x) dx					*/
/*					/ -inf												*/
/*																		*/
/*				x	is level											*/
/*				Q(x) is quantization index								*/
/*				R(k) is reconstructed value								*/
/*				p(x) is probability density function					*/
/*																		*/
/*								   2			   2					*/
/* Solution:		Provided that d (Ln(p(x)) / d x  < 0				*/
/*																		*/
/*						 / x(k+1)				   / x(k+1)				*/
/* (1)		then R(k) = |			  z p(z) dz / |			p(z) dz		*/
/*						/ x(k)					  / x(k)				*/
/*																		*/
/*						{ -inf						k == 0				*/
/* (2)		with x(k) = { 0.5 * (y(k-1) + y(k))		k in ]0, L[			*/
/*						{ +inf						k == L				*/
/*																		*/
/*				L is the number of reconstruction levels				*/
/*				x(k) are the (L + 1) thresholds							*/
/*																		*/
/* Model:		Piecewise constant p(x) given by histogram				*/
/*				Support for p(x) is [xMin - 0.5, xMax + 0.5 [			*/
/*																		*/
/* Comment:		This solution is referenced as Lloyd-Max quantization	*/
/*				For our model, the first condition is not satisfied		*/
/*				  however, you know... it is almost satisfied!			*/
/*																		*/
/* Algorithm:	Compute histogram										*/
/*				Find a good initial solution for x(k)					*/
/*				Iterate (1) and (2) until convergence					*/
/*				Produce the output image with levels = y(k)				*/
/*																		*/
/* Details:		Our initial condition is the set of thresholds which	*/
/*				  produce equi-probable intervals in p(x)				*/
/*				  (technique known as histogram equalization)			*/
/*				Our convergence criterion is quasi-absolute				*/
/*				  (parameter epsilon tolerates some round-off problems)	*/
/*																		*/
/************************************************************************/
int					histQuant		(struct qParam		*quantDataPtr) {

		double				*histogramPtr;
		double				*kernelsPtr;
		float				*floatOutPtr;
		short				*inPtr, *shortOutPtr;
		double				epsilon;
		int					nx, ny, nz;
		short				minLevel;
		unsigned			inLevels, outLevels;

		if (quantDataPtr->floatInput) {
		  message("ERROR - Input has to be of type (short)");
		  return(ERROR);
		}
		if (quantDataPtr->floatOutput && (quantDataPtr->labels != none)) {
		  message("ERROR - Inconsistent parameters");
		  return(ERROR);
		}
		inPtr = quantDataPtr->inPtr.s;
		nx = quantDataPtr->nx;
		ny = quantDataPtr->ny;
		nz = quantDataPtr->nz;
		outLevels = quantDataPtr->bins;
		epsilon = quantDataPtr->epsilon;
		floatOutPtr = quantDataPtr->outPtr.f;
		shortOutPtr = quantDataPtr->outPtr.s;
		if (getHistogram(&histogramPtr, &inLevels, &minLevel, inPtr, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to get histogram");
		  return(ERROR);
		}
		switch (quantDataPtr->mode) {
		  case linear:
			if (quantDataPtr->floatOutput) {
			  rampHistogram(histogramPtr, inLevels, outLevels, minLevel, (short)((long)minLevel
				+ (long)inLevels - 1L));
			  if (buildFloatImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				floatOutPtr, histogramPtr, minLevel) == ERROR) {
				free(histogramPtr);
				message("ERROR - Unable to build (float) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == none) {
			  rampHistogram(histogramPtr, inLevels, outLevels, minLevel, (short)((long)minLevel
				+ (long)inLevels - 1L));
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == from) {
			  if ((long)((short)((long)quantDataPtr->first + (long)quantDataPtr->step
				* ((long)outLevels - 1L))) != ((long)quantDataPtr->first
				+ (long)quantDataPtr->step * ((long)outLevels - 1L))) {
				free(histogramPtr);
				message("ERROR - Too large label step or too many levels");
				return(ERROR);
			  }
			  rampHistogram(histogramPtr, inLevels, outLevels, quantDataPtr->first,
				(short)((long)quantDataPtr->first + (long)quantDataPtr->step * ((long)outLevels
				- 1L)));
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else {
			  free(histogramPtr);
			  message("ERROR - Unknown label tag");
			  return(ERROR);
			}
			break;
		  case equalize:
			kernelsPtr = (double *)malloc((size_t)outLevels * sizeof(double));
			if (kernelsPtr == (double *)NULL) {
			  free(histogramPtr);
			  message("ERROR - Not enough memory for holding kernels");
			  return(ERROR);
			}
			if (equalizedKern(histogramPtr, inLevels, outLevels, kernelsPtr) == ERROR) {
			  free(kernelsPtr);
			  free(histogramPtr);
			  message("ERROR - Unable to find kernels for histogram equalization");
			  return(ERROR);
			}
			if (quantDataPtr->floatOutput) {
			  setLookUpTable(histogramPtr, kernelsPtr, minLevel, inLevels, outLevels, FALSE,
				(short)0);
			  if (buildFloatImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				floatOutPtr, histogramPtr, minLevel) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == none) {
			  setLookUpTable(histogramPtr, kernelsPtr, minLevel, inLevels, outLevels, FALSE,
				(short)0);
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == from) {
			  setLookUpTable(histogramPtr, kernelsPtr, quantDataPtr->first, inLevels, outLevels,
				TRUE, quantDataPtr->step);
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else {
			  free(kernelsPtr);
			  free(histogramPtr);
			  message("ERROR - Unknown label tag");
			  return(ERROR);
			}
			free(kernelsPtr);
			break;
		  case LloydMax:
			kernelsPtr = (double *)malloc((size_t)outLevels * sizeof(double));
			if (kernelsPtr == (double *)NULL) {
			  free(histogramPtr);
			  message("ERROR - Not enough memory for holding kernels");
			  return(ERROR);
			}
			if (findKernels(histogramPtr, inLevels, outLevels, kernelsPtr, epsilon) == ERROR) {
			  free(kernelsPtr);
			  free(histogramPtr);
			  message("ERROR - Unable to find kernels for histogram posterization");
			  return(ERROR);
			}
			if (quantDataPtr->floatOutput) {
			  setLookUpTable(histogramPtr, kernelsPtr, minLevel, inLevels, outLevels, FALSE,
				(short)0);
			  if (buildFloatImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				floatOutPtr, histogramPtr, minLevel) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (float) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == none) {
			  setLookUpTable(histogramPtr, kernelsPtr, minLevel, inLevels, outLevels, FALSE,
				(short)0);
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else if (quantDataPtr->labels == from) {
			  setLookUpTable(histogramPtr, kernelsPtr, quantDataPtr->first, inLevels, outLevels,
				TRUE, quantDataPtr->step);
			  if (buildShortImage(inPtr, &(quantDataPtr->mse), &(quantDataPtr->snr), nx, ny, nz,
				shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
				free(kernelsPtr);
				free(histogramPtr);
				message("ERROR - Unable to build (short) output image");
				return(ERROR);
			  }
			}
			else {
			  free(kernelsPtr);
			  free(histogramPtr);
			  message("ERROR - Unknown label tag");
			  return(ERROR);
			}
			free(kernelsPtr);
			break;
		  default:
			free(histogramPtr);
			message("ERROR - Unknown mode tag");
			return(ERROR);
		}
		free(histogramPtr);
		return(!ERROR);
} /* End of histQuant */

/************************************************************************/
/* FUNCTION: parzenEntropy												*/
/************************************************************************/
int					parzenEntropy	(struct qParam		*entDataPtr) {

		double				*histogram1, *histogram2;
		double				*crossHistogram;
		double				*h1, *h2, *ch;
		float				*inPtr1, *inPtr2;
		float				*p, *q;
		double				scale1, scale2;
		double				mse, snr, H;
		double				nxyz;
		float				min1, min2;
		float				max1, max2;
		int					nx, ny, nz;
		int					x, y, z;
		unsigned			bins;
		unsigned			i, j;

		inPtr1 = entDataPtr->inPtr.f;
		inPtr2 = entDataPtr->outPtr.f;
		nx = entDataPtr->nx;
		ny = entDataPtr->ny;
		nz = entDataPtr->nz;
		nxyz = (double)nx * (double)ny * (double)nz;
		bins = entDataPtr->bins;
		p = inPtr1;
		q = inPtr2;
		min1 = *p;
		max1 = *p;
		min2 = *q;
		max2 = *q;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); p++, q++, x++) {
			  min1 = (*p < min1) ? (*p) : (min1);
			  max1 = (max1 < *p) ? (*p) : (max1);
			  min2 = (*q < min2) ? (*q) : (min2);
			  max2 = (max2 < *q) ? (*q) : (max2);
			}
		if ((min1 == max1) || (min2 == max2)) {
		  message("ERROR - At least one of the input volumes is blank");
		  return(ERROR);
		}
		scale1 = (double)(bins - 1U) / (double)(max1 - min1);
		scale2 = (double)(bins - 1U) / (double)(max2 - min2);
		crossHistogram = (double *)malloc((size_t)bins * (size_t)bins * sizeof(double));
		if (crossHistogram == (double *)NULL) {
		  message("ERROR - Not enough memory for holding crossHistogram");
		  return(ERROR);
		}
		if (crossParzenHistogram(inPtr1, inPtr2, crossHistogram, &mse, &snr,
		  entDataPtr->parzenDegree, scale1, scale2, min1, min2, bins, nx, ny, nz) == ERROR) {
		  free(crossHistogram);
		  message("ERROR - Unable to compute the cross-histogram");
		  return(ERROR);
		}
		histogram1 = (double *)malloc((size_t)bins * sizeof(double));
		if (histogram1 == (double *)NULL) {
		  free(crossHistogram);
		  message("ERROR - Not enough memory for holding output histogram1");
		  return(ERROR);
		}
		histogram2 = (double *)malloc((size_t)bins * sizeof(double));
		if (histogram2 == (double *)NULL) {
		  free(crossHistogram);
		  free(histogram1);
		  message("ERROR - Not enough memory for holding output histogram2");
		  return(ERROR);
		}
		h1 = histogram1;
		for (i = 0U; (i < bins); h1++, i++)
		  *h1 = 0.0;
		ch = crossHistogram;
		h2 = histogram2;
		for (j = 0U; (j < bins); h2++, j++) {
		  *h2 = 0.0;
		  h1 = histogram1;
		  for (i = 0U; (i < bins); ch++, h1++, i++) {
			*h1 += *ch;
			*h2 += *ch;
		  }
		}
		H = 0.0;
		h1 = histogram1;
		h2 = histogram2;
		ch = crossHistogram;
		for (j = 0U; (j < bins); h2++, j++)
		  if (*h2 > (double)FLT_EPSILON) {
			h1 = histogram1;
			for (i = 0U; (i < bins); ch++, h1++, i++)
			  if ((*ch > (double)FLT_EPSILON) && (*h1 > (double)FLT_EPSILON))
				H += *ch * log(*ch / (*h1 * *h2));
		  }
		  else
			ch += (ptrdiff_t)bins;
		H /= log(2.0);
		entDataPtr->entropy = H;
		entDataPtr->mse = mse;
		entDataPtr->snr = snr;
		free(crossHistogram);
		free(histogram2);
		free(histogram1);
		return(!ERROR);
} /* End of parzenEntropy */

/************************************************************************/
/* FUNCTION: quantizeFromBins											*/
/************************************************************************/
void				quantizeFromBins(float				*inPtr,
									 float				*binPtr,
									 short				*outPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 unsigned			bins) {

		float				*leq, *heq;
		int					x, y, z;
		short				b;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, outPtr++, x++) {
			  leq = binPtr + (ptrdiff_t)bins;
			  heq = binPtr;
			  b = (short)0;
			  while (heq < leq)
				if (*heq < *inPtr) {
				  heq++;
				  b++;
				}
				else {
				  if (heq == binPtr)
					leq = heq;
				  else {
					leq = heq - (ptrdiff_t)1;
					b = ((*inPtr + *inPtr) < (*leq + *heq)) ? (b - (short)1) : (b);
				  }
				  break;
				}
			  if (leq == (binPtr + (ptrdiff_t)bins))
				*outPtr = (short)bins;
			  else
				*outPtr = b;
			}
} /* End of quantizeFromBins */

/************************************************************************/
/* FUNCTION: sliceQuant													*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	Contrast enhancement by stretching histogram			*/
/*				  (within a given range)								*/
/*																		*/
/* Solution:	Levels within a range are linearly remapped x -> f(x)	*/
/*				The range is given by endpoints		[iMin, iMax]		*/
/*				The linear operation is given by endpoints mapping		*/
/*				Any value out of range is mapped to some constant		*/
/*																		*/
/*				iMin -> oMin			low endpoint					*/
/*				iMax -> oMax			high endpoint					*/
/*				x in [iMin, iMax]		level							*/
/*																		*/
/*				f(x) = oMin + (oMax - oMin) (x - iMin) / (iMax - iMin)	*/
/*				f(x) = backgrnd		x not in [iMin, iMax]				*/
/*																		*/
/************************************************************************/
int					sliceQuant		(struct qParam		*sliceDataPtr) {

		double				*histogramPtr;
		float				*floatOutPtr;
		short				*inPtr, *shortOutPtr;
		float				oMin, oMax;
		float				loBackgrnd, hiBackgrnd;
		int					nx, ny, nz;
		short				minLevel;
		short				iMin, iMax;
		unsigned			inLevels;

		if (sliceDataPtr->floatInput) {
		  message("ERROR - Input has to be of type (short)");
		  return(ERROR);
		}
		iMin = sliceDataPtr->iMin;
		iMax = sliceDataPtr->iMax;
		oMin = sliceDataPtr->oMin;
		oMax = sliceDataPtr->oMax;
		if (iMin > iMax) {
		  message("ERROR - Parameter iMin should not exceed iMax");
		  return(ERROR);
		}
		if ((iMin == iMax) && (oMin != oMax)) {
		  message("ERROR - Conflicting pair (oMin, oMax)");
		  return(ERROR);
		}
		inPtr = sliceDataPtr->inPtr.s;
		nx = sliceDataPtr->nx;
		ny = sliceDataPtr->ny;
		nz = sliceDataPtr->nz;
		floatOutPtr = sliceDataPtr->outPtr.f;
		shortOutPtr = sliceDataPtr->outPtr.s;
		loBackgrnd = sliceDataPtr->backgrnd;
		hiBackgrnd = sliceDataPtr->foregrnd;
		if (getHistogram(&histogramPtr, &inLevels, &minLevel, inPtr, nx, ny, nz) == ERROR) {
		  message("ERROR - Unable to get histogram");
		  return(ERROR);
		}
		sliceHistogram(histogramPtr, inLevels, minLevel, iMin, iMax, oMin, oMax, loBackgrnd,
		  hiBackgrnd);
		if (sliceDataPtr->floatOutput)
		  if (buildFloatImage(inPtr, &(sliceDataPtr->mse), &(sliceDataPtr->snr), nx, ny, nz,
			floatOutPtr, histogramPtr, minLevel) == ERROR) {
			message("ERROR - Unable to build (float) output image");
			return(ERROR);
		  }
		  else
			;
		else
		  if (buildShortImage(inPtr, &(sliceDataPtr->mse), &(sliceDataPtr->snr), nx, ny, nz,
			shortOutPtr, histogramPtr, minLevel, inLevels) == ERROR) {
			message("ERROR - Unable to build (short) output image");
			return(ERROR);
		  }
		free(histogramPtr);
		return(!ERROR);
} /* End of sliceQuant */

/************************************************************************/
/* FUNCTION: thrQuant													*/
/*----------------------------------------------------------------------*/
/*																		*/
/* Objective:	A priori image segmentation								*/
/*																		*/
/* Solution:	All levels lower or equal to some given threshold		*/
/*				  are replaced by a given backgrnd value				*/
/*				All levels strictly greater than the previous threshold	*/
/*				  are replaced by a given foregrnd value				*/
/*																		*/
/*				x in R				is level							*/
/*				Q(x) in [0, 1]		is quantization index				*/
/*				R(i)				is reconstructed value				*/
/*				T					is threshold						*/
/*																		*/
/*					   { 0		x <= T									*/
/*				Q(x) = {												*/
/*					   { 1		x > T									*/
/*																		*/
/*					   { backgrnd		i == 0							*/
/*				R(i) = {												*/
/*					   { foregrnd		i == 1							*/
/*																		*/
/************************************************************************/
int					thrQuant		(struct qParam		*thrDataPtr) {

		double				*histogramPtr, *hstPtr;
		float				*floatInPtr, *floatOutPtr;
		short				*shortInPtr, *shortOutPtr;
		double				sum;
		double				sumIn, sumNoise, snr;
		double				r, rounded;
		int					x, y, z;
		short				minLevel;
		short				i;
		short				foregrnd, backgrnd;
		unsigned			inLevels, n;

		floatInPtr = thrDataPtr->inPtr.f;
		shortInPtr = thrDataPtr->inPtr.s;
		floatOutPtr = thrDataPtr->outPtr.f;
		shortOutPtr = thrDataPtr->outPtr.s;
		r = (double)thrDataPtr->foregrnd + 0.5;
		rounded = (double)((long)r);
		if ((r < 0.0) && (rounded != r))
		  rounded--;
		if (rounded < (double)SHRT_MIN)
		  foregrnd = (short)SHRT_MIN;
		else if ((double)SHRT_MAX < rounded)
		  foregrnd = (short)SHRT_MAX;
		else
		  foregrnd = (short)rounded;
		r = (double)thrDataPtr->backgrnd + 0.5;
		rounded = (double)((long)r);
		if ((r < 0.0) && (rounded != r))
		  rounded--;
		if (rounded < (double)SHRT_MIN)
		  backgrnd = (short)SHRT_MIN;
		else if ((double)SHRT_MAX < rounded)
		  backgrnd = (short)SHRT_MAX;
		else
		  backgrnd = (short)rounded;
		if (thrDataPtr->floatInput) {
		  sum = sumIn = sumNoise = 0.0;
		  for (z = 0; (z < thrDataPtr->nz); z++)
			for (y = 0; (y < thrDataPtr->ny); y++)
			  for (x = 0; (x < thrDataPtr->nx); floatInPtr++, x++) {
				if (thrDataPtr->floatOutput) {
				  if (*floatInPtr <= thrDataPtr->threshold)
					*floatOutPtr = thrDataPtr->backgrnd;
				  else
					*floatOutPtr = thrDataPtr->foregrnd;
				  sumNoise += ((double)*floatInPtr - (double)*floatOutPtr) * ((double)*floatInPtr
					- (double)*floatOutPtr);
				  floatOutPtr++;
				}
				else {
				  if (*floatInPtr <= thrDataPtr->threshold)
					*shortOutPtr = backgrnd;
				  else
					*shortOutPtr = foregrnd;
				  sumNoise += ((double)*floatInPtr - (double)*shortOutPtr) * ((double)*floatInPtr
					- (double)*shortOutPtr);
				  shortOutPtr++;
				}
				sumIn += (double)*floatInPtr * (double)*floatInPtr;
				sum += 1.0;
			  }
		  if ((sumNoise == 0.0) || (sumIn == 0.0))
			snr = 0.0;
		  else
			snr = 10.0 * log10(sumIn / sumNoise);
		  if (sum == 0.0)
			sumNoise = 0.0;
		  else
			sumNoise /= sum;
		  thrDataPtr->mse = sumNoise;
		  thrDataPtr->snr = snr;
		}
		else {
		  if (getHistogram(&histogramPtr, &inLevels, &minLevel, shortInPtr,
			thrDataPtr->nx, thrDataPtr->ny, thrDataPtr->nz) == ERROR) {
			message("ERROR - Unable to get histogram");
			return(ERROR);
		  }
		  i = minLevel;
		  hstPtr = histogramPtr;
		  for (n = 0U; (n < inLevels); hstPtr++, n++, i++) {
			r = (double)thrDataPtr->threshold + 0.5;
			rounded = (double)((long)r);
			if ((r < 0.0) && (rounded != r))
			  rounded--;
			if (i <= (short)rounded)
			  *hstPtr = (double)thrDataPtr->backgrnd;
			else
			  *hstPtr = (double)thrDataPtr->foregrnd;
		  }
		  if (thrDataPtr->floatOutput)
			if (buildFloatImage(shortInPtr, &(thrDataPtr->mse), &(thrDataPtr->snr),
			  thrDataPtr->nx, thrDataPtr->ny, thrDataPtr->nz, floatOutPtr, histogramPtr, minLevel)
			  == ERROR) {
			  message("ERROR - Unable to build output (float) image");
			  return(ERROR);
			}
			else
			  ;
		  else
			if (buildShortImage(shortInPtr, &(thrDataPtr->mse), &(thrDataPtr->snr),
			  thrDataPtr->nx, thrDataPtr->ny, thrDataPtr->nz, shortOutPtr, histogramPtr, minLevel,
			  inLevels) == ERROR) {
			  message("ERROR - Unable to build output (short) image");
			  return(ERROR);
			}
		  free(histogramPtr);
		}
		return(!ERROR);
} /* End of thrQuant */
