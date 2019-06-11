#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<math.h>
#include		<float.h>
#include		"phil.h"

#include		"register.h"
#include		"reg1.h"

#include		"quant.h"
#include		"reg0.h"

/************************************************************************/
static	void		lineReduce		(float				xIn[],
									 long				nxIn,
									 float				xOut[],
									 float				downFilter[],
									 short				downTapsNb);
static	int			reduce2D		(float				*inPtr,
									 long				nxIn,
									 long				nyIn,
									 long				nzIn,
									 float				*outPtr,
									 long				nxOut,
									 long				nyOut,
									 long				nzOut,
									 float				downFilter[],
									 short				downTapsNb);
static	int			reduce3D		(float				*inPtr,
									 long				nxIn,
									 long				nyIn,
									 long				nzIn,
									 float				*outPtr,
									 long				nxOut,
									 long				nyOut,
									 long				nzOut,
									 float				downFilter[],
									 short				downTapsNb);

/************************************************************************/
static	void		lineReduce		(float				xIn[],
									 long				nxIn,
									 float				xOut[],
									 float				downFilter[],
									 short				downTapsNb) {

		float				*d;
		float				x;
		long				k, i, i1, i2, kk, kn;

		if (nxIn == 1L) {
		  xOut[0] = xIn[0];
		  return;
		}
		kn = 2L * nxIn - 2L;
		if (downTapsNb == (short)1)
		  for (kk = 0L; ((kk + kk) < nxIn); kk++) {
			k = kk + kk;
			i2 = k + 1L;
			if (i2 >= nxIn)
			  i2 = kn - i2;
			xOut[kk] = (xIn[k] + xIn[i2]) * 0.5F;
		  }
		else
		  for (kk = 0L; ((kk + kk) < nxIn); kk++) {
			k= kk + kk;
			d = downFilter;
			x = *d * xIn[k];
			for (d++, i = 1L; (i < (long)downTapsNb); d++, i++) {
			  i1 = k - i;
			  i2 = k + i;
			  if (i1 < 0L) {
				i1 = (-i1) - kn * ((-i1) / kn);
				if (i1 >= nxIn)
				  i1 = kn - i1;
			  }
			  if (i2 >= nxIn) {
				i2 -= kn * (i2 / kn);
				if (i2 >= nxIn)
				  i2 = kn - i2;
			  }
			  x += *d * (xIn[i1] + xIn[i2]);
			}
			xOut[kk] = x;
		  }
} /* End of lineReduce */

/************************************************************************/
static	int			reduce2D		(float				*inPtr,
									 long				nxIn,
									 long				nyIn,
									 long				nzIn,
									 float				*outPtr,
									 long				nxOut,
									 long				nyOut,
									 long				nzOut,
									 float				downFilter[],
									 short				downTapsNb) {

		float				*tmp, *p;
		float				*pIn, *pOut;
		float				*yIn, *yOut;
		long				i, j;

		if (nzIn != nzOut) {
		  message("ERROR - Incoherent 2D pyramid reduction parameters");
		  return(ERROR);
		}
		tmp = (float *)malloc((size_t)nxOut * (size_t)nyIn * sizeof(float));
		if (tmp == (float *)NULL) {
		  message("ERROR - Not enough memory for 2-D pyramid workspace");
		  return(ERROR);
		}
		yIn = (float *)malloc((size_t)nyIn * sizeof(float));
		if (yIn == (float *)NULL) {
		  free(tmp);
		  message("ERROR - Not enough memory for extracting column");
		  return(ERROR);
		}
		yOut = (float *)malloc((size_t)nyOut * sizeof(float));
		if (yOut == (float *)NULL) {
		  free(yIn);
		  free(tmp);
		  message("ERROR - Not enough memory for saving column");
		  return(ERROR);
		}
		for (nzIn = 0L; (nzIn < nzOut); nzIn++) {
		  pIn = inPtr + (ptrdiff_t)nxIn * (ptrdiff_t)nyIn * (ptrdiff_t)nzIn;
		  pOut = tmp;
		  for (j = 0L; (j < nyIn); pIn += (ptrdiff_t)nxIn, pOut += (ptrdiff_t)nxOut, j++)
			lineReduce(pIn, nxIn, pOut, downFilter, downTapsNb);
		  pIn = tmp;
		  pOut = outPtr + (ptrdiff_t)nxOut * (ptrdiff_t)nyOut * (ptrdiff_t)nzIn;
		  for (i = 0L; (i < nxOut); pIn++, pOut++, i++) {
			p = pIn;
			for (j = 0L; (j < nyIn); p += (ptrdiff_t)nxOut, j++)
			  yIn[j] = *p;
			lineReduce(yIn, nyIn, yOut, downFilter, downTapsNb);
			p = pOut;
			for (j = 0L; (j < nyOut); p += (ptrdiff_t)nxOut, j++)
			  *p = yOut[j];
		  }
		}
		free(yOut);
		free(yIn);
		free(tmp);
		return(!ERROR);
} /* End of reduce2D */

/************************************************************************/
static	int			reduce3D		(float				*inPtr,
									 long				nxIn,
									 long				nyIn,
									 long				nzIn,
									 float				*outPtr,
									 long				nxOut,
									 long				nyOut,
									 long				nzOut,
									 float				downFilter[],
									 short				downTapsNb) {

		float				*tmp;
		float				*p, *q;
		float				*pIn, *pOut;
		float				*yIn, *yOut, *zIn, *zOut;
		long				yInxOut = nyIn * nxOut;
		long				yOutxOut = nyOut * nxOut;
		long				i, j, k;

		tmp = (float *)malloc((size_t)yInxOut * (size_t)nzIn * sizeof(float));
		if (tmp == (float *)NULL) {
		  message("ERROR - Not enough memory for 3-D pyramid workspace");
		  return(ERROR);
		}
		yIn = (float *)malloc((size_t)nyIn * sizeof(float));
		if (yIn == (float *)NULL) {
		  free(tmp);
		  message("ERROR - Not enough memory for extracting column");
		  return(ERROR);
		}
		yOut = (float *)malloc((size_t)nyOut * sizeof(float));
		if (yOut == (float *)NULL) {
		  free(yIn);
		  free(tmp);
		  message("ERROR - Not enough memory for saving column");
		  return(ERROR);
		}
		zIn = (float *)malloc((size_t)nzIn * sizeof(float));
		if (zIn == (float *)NULL) {
		  free(yOut);
		  free(yIn);
		  free(tmp);
		  message("ERROR - Not enough memory for extracting tower");
		  return(ERROR);
		}
		zOut = (float *)malloc((size_t)nzOut * sizeof(float));
		if (zOut == (float *)NULL) {
		  free(zIn);
		  free(yOut);
		  free(yIn);
		  free(tmp);
		  message("ERROR - Not enough memory for saving tower");
		  return(ERROR);
		}
		pIn = inPtr;
		pOut = tmp;
		for (k = 0L; (k < nzIn); k++)
		  for (j = 0L; (j < nyIn); pIn += (ptrdiff_t)nxIn, pOut += (ptrdiff_t)nxOut, j++)
			lineReduce(pIn, nxIn, pOut, downFilter, downTapsNb);
		pIn = tmp;
		for (k = 0L; (k < nzIn); pIn += (ptrdiff_t)yInxOut, k++) {
		  p = pIn;
		  for (i = 0L; (i < nxOut); p++, i++) {
			q = p;
			for (j = 0L; (j < nyIn); q += (ptrdiff_t)nxOut, j++)
			  yIn[j] = *q;
			lineReduce(yIn, nyIn, yOut, downFilter, downTapsNb);
			pOut = p;
			for (j = 0L; (j < nyOut); pOut += (ptrdiff_t)nxOut, j++)
			  *pOut = yOut[j];
		  }
		}
		pIn = tmp;
		pOut = outPtr;
		for (j = 0L; (j < nyOut); j++)
		  for (i = 0L; (i < nxOut); pIn++, pOut++, i++) {
			p = pIn;
			for (k = 0L; (k < nzIn); p += (ptrdiff_t)yInxOut, k++)
			  zIn[k] = *p;
			lineReduce(zIn, nzIn, zOut, downFilter, downTapsNb);
			q = pOut;
			for (k = 0L; (k < nzOut); q += (ptrdiff_t)yOutxOut, k++)
			  *q = zOut[k];
		  }
		free(zOut);
		free(zIn);
		free(yOut);
		free(yIn);
		free(tmp);
		return(!ERROR);
} /* End of reduce3D */

/************************************************************************/
/* FUNCTION: convertOrigin												*/
/************************************************************************/
void				convertOrigin	(struct	fitRec		*fit,
									 double				origin[]) {

		double				dx[3], oldOrigin[3];
		int					i, j;

		for (i = 0; (i < 3); i++) {
		  oldOrigin[i] = fit->origin[i];
		  dx[i] = oldOrigin[i];
		}
		for (i = 0; (i < 3); i++) {
		  for (j = 0; (j < 3); j++)
			fit->dx[i] += fit->skew[i][j] * dx[j];
		  fit->dx[i] -= dx[i];
		}
		for (i = 0; (i < 3); i++) {
		  fit->origin[i] = origin[i];
		  dx[i] = fit->origin[i];
		}
		for (i = 0; (i < 3); i++) {
		  for (j = 0; (j < 3); j++)
			fit->dx[i] -= fit->skew[i][j] * dx[j];
		  fit->dx[i] += dx[i];
		}
		for (i = 0; (i < 3); i++)
		  origin[i] = oldOrigin[i];
} /* End of convertOrigin */

/************************************************************************/
/* FUNCTION: downscaleFit												*/
/************************************************************************/
void				downscaleFit	(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze) {

		if (nxRange[fromLevel + 1] > 1L) {
		  fit->dx[0] /= 2.0;
		  fit->origin[0] /= 2.0;
		}
		else {
		  fit->skew[1][0] /= 2.0;
		  fit->skew[2][0] /= 2.0;
		  fit->skew[0][1] *= 2.0;
		  fit->skew[0][2] *= 2.0;
		}
		if (nyRange[fromLevel + 1] > 1L) {
		  fit->dx[1] /= 2.0;
		  fit->origin[1] /= 2.0;
		}
		else {
		  fit->skew[0][1] /= 2.0;
		  fit->skew[2][1] /= 2.0;
		  fit->skew[1][0] *= 2.0;
		  fit->skew[1][2] *= 2.0;
		}
		if (zSqueeze && (nzRange[fromLevel + 1] > 1L)) {
		  fit->dx[2] /= 2.0;
		  fit->origin[2] /= 2.0;
		}
		else {
		  fit->skew[0][2] /= 2.0;
		  fit->skew[1][2] /= 2.0;
		  fit->skew[2][0] *= 2.0;
		  fit->skew[2][1] *= 2.0;
		}
} /* End of downscaleFit */

/************************************************************************/
/* FUNCTION: freePyramids												*/
/************************************************************************/
void				freePyramids	(float				*pyrHdl1[],
									 float				*pyrHdl2[],
									 float				*mskHdl1[],
									 float				*mskHdl2[],
									 int				levels) {

		int					l;

		for (l = 0; (l < levels); l++) {
		  if (mskHdl2[l] != (float *)NULL)
			free(mskHdl2[l]);
		  mskHdl2[l] = (float *)NULL;
		  if (mskHdl1[l] != (float *)NULL)
			free(mskHdl1[l]);
		  mskHdl1[l] = (float *)NULL;
		  if (pyrHdl2[l] != (float *)NULL)
			free(pyrHdl2[l]);
		  pyrHdl2[l] = (float *)NULL;
		  if (pyrHdl1[l] != (float *)NULL)
			free(pyrHdl1[l]);
		  pyrHdl1[l] = (float *)NULL;
		}
} /* End of freePyramids */

/************************************************************************/
/* FUNCTION: pyramid													*/
/************************************************************************/
int					pyramid			(float				*inPtr,
									 int				levels,
									 float				*pyrHdl[],
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				nx,
									 int				ny,
									 int				nz,
									 int				zSqueeze,
									 enum	iDegree		interpolation) {

		float				downFilter[maxTapsNumber], upFilter[maxTapsNumber];
		short				downTapsNb, upTapsNb;
		int					l;
		int					error;

		switch (interpolation) {
		  case zero:
			pyr_filters(downFilter, &downTapsNb, upFilter, &upTapsNb, (short)0);
			break;
		  case one:
			pyr_filters(downFilter, &downTapsNb, upFilter, &upTapsNb, (short)1);
			break;
		  case three:
			pyr_filters(downFilter, &downTapsNb, upFilter, &upTapsNb, (short)3);
			break;
		  default:
			message("ERROR - Invalid interpolation degree for pyramid computation");
			return(ERROR);
		}
		levels = pyr_getsize(nx, ny, nz, levels, nxRange, nyRange, nzRange);
		if (!zSqueeze)
		  for (l = 0; (l < levels); l++)
			nzRange[l] = (long)nz;
		for (l = 0; (l < levels); l++) {
		  pyrHdl[l] = (float *)malloc((size_t)nxRange[l] * (size_t)nyRange[l]
			* (size_t)nzRange[l] * sizeof(float));
		  if (pyrHdl[l] == (float *)NULL) {
			message("ERROR - Not enough memory for holding pyramid");
			return(ERROR);
		  }
		}
		copyClip(inPtr, nx, ny, nz, pyrHdl[0], (int)nxRange[0], (int)nyRange[0], (int)nzRange[0]);
		error = !ERROR;
		if (zSqueeze)
		  for (l = 1; (l < levels); l++)
			error = error || reduce3D(pyrHdl[l - 1], nxRange[l - 1], nyRange[l - 1],
			  nzRange[l - 1], pyrHdl[l], nxRange[l], nyRange[l], nzRange[l],
			  downFilter, downTapsNb);
		else
		  for (l = 1; (l < levels); l++)
			error = error || reduce2D(pyrHdl[l - 1], nxRange[l - 1], nyRange[l - 1],
			  nzRange[l - 1], pyrHdl[l], nxRange[l], nyRange[l], nzRange[l],
			  downFilter, downTapsNb);
		return(error);
} /* End of pyramid */

/************************************************************************/
/* FUNCTION: upscaleFit													*/
/************************************************************************/
void				upscaleFit		(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze) {

		if (nxRange[fromLevel - 1] > 1L) {
		  fit->dx[0] *= 2.0;
		  fit->origin[0] *= 2.0;
		}
		else {
		  fit->skew[1][0] *= 2.0;
		  fit->skew[2][0] *= 2.0;
		  fit->skew[0][1] /= 2.0;
		  fit->skew[0][2] /= 2.0;
		}
		if (nyRange[fromLevel - 1] > 1L) {
		  fit->dx[1] *= 2.0;
		  fit->origin[1] *= 2.0;
		}
		else {
		  fit->skew[0][1] *= 2.0;
		  fit->skew[2][1] *= 2.0;
		  fit->skew[1][0] /= 2.0;
		  fit->skew[1][2] /= 2.0;
		}
		if (zSqueeze && (nzRange[fromLevel - 1] > 1L)) {
		  fit->dx[2] *= 2.0;
		  fit->origin[2] *= 2.0;
		}
		else {
		  fit->skew[0][2] *= 2.0;
		  fit->skew[1][2] *= 2.0;
		  fit->skew[2][0] /= 2.0;
		  fit->skew[2][1] /= 2.0;
		}
} /* End of upscaleFit */
