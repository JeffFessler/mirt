#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<string.h>
#include		<math.h>
#include		<float.h>
#include		"phil.h"

#include		"quant.h"
#include		"register.h"
#include		"regFlt3d.h"

#include		"reg0.h"
#include		"reg1.h"
#include		"reg2.h"

/************************************************************************/
static	double		averageMskData	(float				*inPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz);
static	void		bilevelFormat	(float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz);
static	int			getCenter		(float				*inPtr,
									 float				*inMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				*Gx,
									 double				*Gy,
									 double				*Gz);
static	void		logicalFormat	(float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz);
static	void		zapMean			(float				*inPtr,
									 double				mean,
									 int				nx,
									 int				ny,
									 int				nz);

/************************************************************************/
static	double		averageMskData	(float				*inPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				n = 0.0;
		double				mean = 0.0;
		int					x, y, z;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, mskPtr++, x++)
			  if ((int)*mskPtr) {
				mean += (double)*inPtr;
				n += 1.0;
			  }
		if (n != 0.0)
		  mean /= n;
		return(mean);
} /* End of averageMskData */

/************************************************************************/
static	void		bilevelFormat	(float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz) {

		int					x, y, z;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); mskPtr++, x++)
			  *mskPtr = (*mskPtr != 0.0F) ? ((float)ONlevel) : ((float)OFFlevel);
} /* End of bilevelFormat */

/************************************************************************/
static	int			getCenter		(float				*inPtr,
									 float				*inMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				*Gx,
									 double				*Gy,
									 double				*Gz) {

		double				w;
		double				xO, yO, zO;
		double				x, y, z;

		w = 0.0;
		*Gx = *Gy = *Gz = 0.0;
		xO = (double)((nx - 1) / 2);
		yO = (double)((ny - 1) / 2);
		zO = (double)((nz - 1) / 2);
		for (z = -zO; (z < (nz - zO)); z++)
		  for (y = -yO; (y < (ny - yO)); y++)
			for (x = -xO; (x < (nx - xO)); inPtr++, inMsk++, x++)
			  if ((int)*inMsk) {
				*Gx += x * (double)*inPtr;
				*Gy += y * (double)*inPtr;
				*Gz += z * (double)*inPtr;
				w += (double)*inPtr;
			  }
		if ((w * w) < (double)FLT_EPSILON) {
		  message("ERROR - Weightless image");
		  return(ERROR);
		}
		*Gx /= w;
		*Gy /= w;
		*Gz /= w;
		return(!ERROR);
} /* End of getCenter */

/************************************************************************/
static	void		logicalFormat	(float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz) {

		int					x, y, z;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); mskPtr++, x++)
			  *mskPtr = (float)(*mskPtr >= (((float)OFFlevel + (float)ONlevel) / 2.0F));
} /* End of logicalFormat */

/************************************************************************/
static	void		zapMean			(float				*inPtr,
									 double				mean,
									 int				nx,
									 int				ny,
									 int				nz) {

		int					x, y, z;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, x++)
			  *inPtr -= (float)mean;
} /* End of zapMean */

/************************************************************************/
/* FUNCTION: regFloat3d													*/
/************************************************************************/
int					regFloat3d		(struct	rParam		*regDataPtr) {

		struct	fitRec		fit, invFit;
		float				*pyrHdl1[maxPyrLevels], *pyrHdl2[maxPyrLevels];
		float				*mskHdl1[maxPyrLevels], *mskHdl2[maxPyrLevels];
		float				*inPtr1, *inPtr2;
		float				*outPtr;
		float				*mskPtr1, *mskPtr2;
		float				*mskPtr;
		float				*u;
		double				origin[3];
		double				epsilon[maxPyrLevels];
		long				nxRange[maxPyrLevels];
		long				nyRange[maxPyrLevels];
		long				nzRange[maxPyrLevels];
		int					ia[hessianSize];
		double				mse, snr;
		double				lambda;
		double				x1, y1, z1;
		double				x2, y2, z2;
		double				mean;
		int					nx, ny, nz;
		int					k, l;
		int					x, y, z;

		switch (regDataPtr->directives.convergence) {
		  case gravity:
			if (regDataPtr->directives.isoScaling || regDataPtr->directives.xRot
			  || regDataPtr->directives.yRot || regDataPtr->directives.zRot
			  || regDataPtr->directives.xSkew || regDataPtr->directives.ySkew
			  || regDataPtr->directives.zSkew || regDataPtr->directives.matchGrey) {
			  message("ERROR - Gravity center can only be translated");
			  return(ERROR);
			}
			else if ((regDataPtr->directives.importFit != 0) && (regDataPtr->directives.xTrans
			  || regDataPtr->directives.yTrans || regDataPtr->directives.zTrans))
			  message("WARNING - Ignoring 'importFit' parameter");
			break;
		  case Marquardt:
			break;
		  default:
			message("ERROR - Unknown convergence specification");
			return(ERROR);
		}
		switch (regDataPtr->directives.interpolation) {
		  case zero:
			if (regDataPtr->directives.convergence != gravity) {
			  message("ERROR - Nearest neighbor interpolation valid only for gravity");
			  return(ERROR);
			}
		  case one:
		  case three:
			break;
		  default:
			message("ERROR - Unknown interpolation specification");
			return(ERROR);
		}
		if (regDataPtr->directives.zSqueeze && (regDataPtr->nz == 1)) {
		  message("ERROR - Inconsistent parameters: unable to shrink Z-axis in 2-D");
		  return(ERROR);
		}
		if ((!regDataPtr->directives.zSqueeze) && (regDataPtr->nz != 1)
		  && regDataPtr->directives.isoScaling) {
		  message("ERROR - Inconsistent parameters: zSqueeze needed for isoScaling in 3-D");
		  return(ERROR);
		}
		if (regDataPtr->directives.isoScaling && (regDataPtr->directives.xSkew
		  || regDataPtr->directives.ySkew || regDataPtr->directives.zSkew)) {
		  message("ERROR - Inconsistent parameters: (x,y,z)Skew and isoScaling");
		  return(ERROR);
		}
		if ((regDataPtr->directives.xRot || regDataPtr->directives.yRot
		  || regDataPtr->directives.zRot) && (regDataPtr->directives.xSkew
		  || regDataPtr->directives.ySkew || regDataPtr->directives.zSkew)) {
		  message("ERROR - Inconsistent parameters: (x,y,z)Rot and (x,y,z)Skew");
		  return(ERROR);
		}
		if (regDataPtr->directives.xRot && (!regDataPtr->directives.zSqueeze)) {
		  message("ERROR - Inconsistent parameters: xRot needs zSqueeze");
		  return(ERROR);
		}
		if (regDataPtr->directives.yRot && (!regDataPtr->directives.zSqueeze)) {
		  message("ERROR - Inconsistent parameters: yRot needs zSqueeze");
		  return(ERROR);
		}
		message("STATUS - Initializing...");
		inPtr1 = regDataPtr->inPtr1;
		inPtr2 = regDataPtr->inPtr2;
		mskPtr1 = regDataPtr->inMsk1;
		mskPtr2 = regDataPtr->inMsk2;
		outPtr = regDataPtr->outPtr;
		mskPtr = regDataPtr->mskPtr;
		nx = regDataPtr->nx;
		ny = regDataPtr->ny;
		nz = regDataPtr->nz;
		initialEstimate(&fit, nx, ny, nz);
		if ((regDataPtr->directives.importFit != 0)
		  && (regDataPtr->directives.convergence != gravity)) {
		  if (importFit(&fit, regDataPtr->inFit) == ERROR) {
			message("ERROR - Importation of fit variables failed");
			return(ERROR);
		  }
		  switch (regDataPtr->directives.importFit) {
			case -1:
			  if (invertFit(&fit, &invFit) == ERROR) {
				message("ERROR - Unable to compute the backward transformation");
				return(ERROR);
			  }
			  fit = invFit;
			  break;
			case 1:
			  break;
			default:
			  message("ERROR - Unknown importFit specification");
			  return(ERROR);
		  }
		}
		switch (regDataPtr->directives.testMask) {
		  case blank:
			u = mskPtr1;
			for (z = 0; (z < nz); z++)
			  for (y = 0; (y < ny); y++)
				for (x = 0; (x < nx); u++, x++)
				  *u = 1.0F;
			break;
		  case provided:
			break;
		  case computed:
			if (maskFromData(inPtr1, mskPtr1, nx, ny, nz, regDataPtr->sx, regDataPtr->sy,
			  regDataPtr->sz) == ERROR) {
			  message("ERROR - Unable to extract mask from test data");
			  return(ERROR);
			}
			break;
		  default:
			message("ERROR - Unknown test mask specification");
			return(ERROR);
		}
		switch (regDataPtr->directives.referenceMask) {
		  case blank:
			u = mskPtr2;
			for (z = 0; (z < nz); z++)
			  for (y = 0; (y < ny); y++)
				for (x = 0; (x < nx); u++, x++)
				  *u = 1.0F;
			break;
		  case provided:
			break;
		  case computed:
			if (maskFromData(inPtr2, mskPtr2, nx, ny, nz, regDataPtr->sx, regDataPtr->sy,
			  regDataPtr->sz) == ERROR) {
			  message("ERROR - Unable to extract mask from reference data");
			  return(ERROR);
			}
			break;
		  default:
			message("ERROR - Unknown reference mask specification");
			return(ERROR);
		}
		if (regDataPtr->directives.zapMean) {
		  mean = averageMskData(inPtr1, mskPtr1, nx, ny, nz);
		  zapMean(inPtr1, mean, nx, ny, nz);
		  mean = averageMskData(inPtr2, mskPtr2, nx, ny, nz);
		  zapMean(inPtr2, mean, nx, ny, nz);
		}
		if ((!regDataPtr->directives.isoScaling) && (!regDataPtr->directives.xTrans)
		  && (!regDataPtr->directives.yTrans) && (!regDataPtr->directives.zTrans)
		  && (!regDataPtr->directives.xRot) && (!regDataPtr->directives.yRot)
		  && (!regDataPtr->directives.zRot) && (!regDataPtr->directives.xSkew)
		  && (!regDataPtr->directives.ySkew) && (!regDataPtr->directives.zSkew)
		  && (!regDataPtr->directives.matchGrey)) {
		  message("WARNING - Nothing to optimize");
		  message("STATUS - Creating output...");
		  if (dirMskTransform(&fit, inPtr1, outPtr, mskPtr1, mskPtr, nx, ny, nz,
			regDataPtr->directives.greyRendering, regDataPtr->directives.interpolation)
			== ERROR) {
			message("ERROR - Final transformation of test image failed");
			return(ERROR);
		  }
		  adornImage(outPtr, mskPtr, nx, ny, nz, nx, ny, nz, regDataPtr->directives.clipping,
			regDataPtr->backgrnd);
		  computeMskSnr(&mse, &snr, inPtr2, outPtr, mskPtr2, mskPtr,
			regDataPtr->directives.maskCombine, nx, ny, nz);
		  if (exportFit(&fit) == ERROR ) {
			message("ERROR - Assignment of transformation variables failed");
			return(ERROR);
		  }
		  if (exportSummary(&fit, mse, snr) == ERROR ) {
			message("ERROR - Assignment of summary variables failed");
			return(ERROR);
		  }
		  if (regDataPtr->directives.exportFit)
			if (fExportFit(&fit, regDataPtr->outFit) == ERROR ) {
			  message("ERROR - Exportation of transformation variables failed");
			  return(ERROR);
			}
		  return(!ERROR);
		}
		bilevelFormat(mskPtr1, nx, ny, nz);
		bilevelFormat(mskPtr2, nx, ny, nz);
		if (regDataPtr->directives.convergence == gravity) {
		  logicalFormat(mskPtr1, nx, ny, nz);
		  logicalFormat(mskPtr2, nx, ny, nz);
		  message("STATUS - Aligning centers of gravity...");
		  if (getCenter(inPtr1, mskPtr1, nx, ny, nz, &x1, &y1, &z1) == ERROR) {
			message("ERROR - Unable to compute 1st center of gravity");
			return(ERROR);
		  }
		  if (getCenter(inPtr2, mskPtr2, nx, ny, nz, &x2, &y2, &z2) == ERROR) {
			message("ERROR - Unable to compute 2nd center of gravity");
			return(ERROR);
		  }
		  fit.dx[0] = x2 - x1;
		  fit.dx[1] = y2 - y1;
		  fit.dx[2] = z2 - z1;
		  message("STATUS - Creating output...");
		  if (dirMskTransform(&fit, inPtr1, outPtr, mskPtr1, mskPtr, nx, ny, nz,
			regDataPtr->directives.greyRendering, regDataPtr->directives.interpolation)
			== ERROR) {
			message("ERROR - Final transformation of test image failed");
			return(ERROR);
		  }
		  adornImage(outPtr, mskPtr, nx, ny, nz, nx, ny, nz, regDataPtr->directives.clipping,
			regDataPtr->backgrnd);
		  computeMskSnr(&mse, &snr, inPtr2, outPtr, mskPtr2, mskPtr,
			regDataPtr->directives.maskCombine, nx, ny, nz);
		  if (exportFit(&fit) == ERROR ) {
			message("ERROR - Assignment of transformation variables failed");
			return(ERROR);
		  }
		  if (exportSummary(&fit, mse, snr) == ERROR ) {
			message("ERROR - Assignment of summary variables failed");
			return(ERROR);
		  }
		  if (regDataPtr->directives.exportFit)
			if (fExportFit(&fit, regDataPtr->outFit) == ERROR ) {
			  message("ERROR - Exportation of transformation variables failed");
			  return(ERROR);
			}
		  return(!ERROR);
		}
		for (l = 0; (l < regDataPtr->levels); l++) {
		  pyrHdl1[l] = (float *)NULL;
		  pyrHdl2[l] = (float *)NULL;
		  mskHdl1[l] = (float *)NULL;
		  mskHdl2[l] = (float *)NULL;
		}
		if (pyramid(mskPtr1, regDataPtr->levels, mskHdl1, nxRange, nyRange, nzRange, nx, ny, nz,
		  regDataPtr->directives.zSqueeze, regDataPtr->directives.interpolation) == ERROR) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Not enough memory for 1st mask pyramid computation");
		  return(ERROR);
		}
		if (pyramid(mskPtr2, regDataPtr->levels, mskHdl2, nxRange, nyRange, nzRange, nx, ny, nz,
		  regDataPtr->directives.zSqueeze, regDataPtr->directives.interpolation) == ERROR) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Not enough memory for 2nd mask pyramid computation");
		  return(ERROR);
		}
		logicalFormat(mskPtr1, nx, ny, nz);
		logicalFormat(mskPtr2, nx, ny, nz);
		for (l = 0; (l < regDataPtr->levels); l++) {
		  logicalFormat(mskHdl1[l], (int)nxRange[l], (int)nyRange[l], (int)nzRange[l]);
		  logicalFormat(mskHdl2[l], (int)nxRange[l], (int)nyRange[l], (int)nzRange[l]);
		}
		if (pyramid(inPtr1, regDataPtr->levels, pyrHdl1, nxRange, nyRange, nzRange, nx, ny, nz,
		  regDataPtr->directives.zSqueeze, regDataPtr->directives.interpolation) == ERROR) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Not enough memory for 1st pyramid computation");
		  return(ERROR);
		}
		if (pyramid(inPtr2, regDataPtr->levels, pyrHdl2, nxRange, nyRange, nzRange, nx, ny, nz,
		  regDataPtr->directives.zSqueeze, regDataPtr->directives.interpolation) == ERROR ) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Not enough memory for 2nd pyramid computation");
		  return(ERROR);
		}
		epsilon[0] = regDataPtr->epsilon;
		switch (regDataPtr->directives.convergence) {
		  case Marquardt:
			for (l = 1; (l < regDataPtr->levels); l++)
			  epsilon[l] = epsilon[l - 1];
			break;
		  default:
			freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
			message("ERROR - Unknown convergence");
			return(ERROR);
		}
		message("STATUS - Optimizing...");
		origin[0] = (double)((nxRange[0] - 1L) / 2L);
		origin[1] = (double)((nyRange[0] - 1L) / 2L);
		origin[2] = (double)((nzRange[0] - 1L) / 2L);
		convertOrigin(&fit, origin);
		for (l = 0; (l < (regDataPtr->levels - 1)); l++)
		  downscaleFit(&fit, nxRange, nyRange, nzRange, l, regDataPtr->directives.zSqueeze);
		lambda = regDataPtr->firstLambda;
		for (l = regDataPtr->levels - 1; (l >= 0); l--) {
		  if ((l == (regDataPtr->levels - 1)) || (l >= (regDataPtr->lastLevel - 1))) {
			ia[0] = regDataPtr->directives.xTrans;							/* dx[0] */
			ia[1] = regDataPtr->directives.yTrans && (nyRange[l] != 1L);	/* dx[1] */
			ia[2] = regDataPtr->directives.zTrans && (nzRange[l] != 1L);	/* dx[2] */
			for (k = 3; (k < 12); k++)
			  ia[k] = FALSE;
			ia[12] = regDataPtr->directives.matchGrey;						/* gamma */
			if (regDataPtr->directives.xSkew) {
			  ia[3] = TRUE;													/* skew[0][0] */
			  ia[6] = (nyRange[l] != 1L);									/* skew[1][0] */
			  ia[9] = (nzRange[l] != 1L);									/* skew[2][0] */
			}
			if (regDataPtr->directives.ySkew) {
			  ia[4] = (nyRange[l] != 1L);									/* skew[0][1] */
			  ia[7] = (nyRange[l] != 1L);									/* skew[1][1] */
			  ia[10] = (nzRange[l] != 1L);									/* skew[2][1] */
			}
			if (regDataPtr->directives.zSkew) {
			  ia[5] = (nzRange[l] != 1L);									/* skew[0][2] */
			  ia[8] = (nzRange[l] != 1L);									/* skew[1][2] */
			  ia[11] = (nzRange[l] != 1L);									/* skew[2][2] */
			}
			if (regDataPtr->directives.xRot || regDataPtr->directives.yRot
			  || regDataPtr->directives.zRot || regDataPtr->directives.isoScaling) {
			  ia[3] = regDataPtr->directives.xRot && (nyRange[l] != 1L);	/* phi */
			  ia[4] = regDataPtr->directives.yRot && (nzRange[l] != 1L);	/* theta */
			  ia[5] = regDataPtr->directives.zRot && (nyRange[l] != 1L);	/* psi */
			  ia[6] = regDataPtr->directives.isoScaling;					/* lambda */
			  ia[7] = regDataPtr->directives.matchGrey;						/* gamma */
			  ia[8] = FALSE;
			  ia[9] = FALSE;
			  ia[10] = FALSE;
			  ia[11] = FALSE;
			  ia[12] = FALSE;
			}
			switch (regDataPtr->directives.convergence) {
			  case Marquardt:
				if (optimize(&fit, &(regDataPtr->directives), ia, pyrHdl1[l], pyrHdl2[l],
				  mskHdl1[l], mskHdl2[l], outPtr, mskPtr, &lambda, regDataPtr->firstLambda,
				  regDataPtr->lambdaScale, epsilon[l], regDataPtr->minGain,
				  (int)nxRange[l], (int)nyRange[l], (int)nzRange[l]) == ERROR) {
				  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
				  message("ERROR - Optimization failed (Marquardt)");
				  return(ERROR);
				}
				lambda = regDataPtr->firstLambda;
				break;
			  default:
				freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
				message("ERROR - Unknown convergence");
				return(ERROR);
			}
		  }
		  if (l != 0)
			upscaleFit(&fit, nxRange, nyRange, nzRange, l, regDataPtr->directives.zSqueeze);
		}
		convertOrigin(&fit, origin);
		message("STATUS - Creating output...");
		if (dirMskTransform(&fit, inPtr1, outPtr, mskPtr1, mskPtr, nx, ny, nz,
		  regDataPtr->directives.greyRendering, regDataPtr->directives.interpolation) == ERROR) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Final transformation of test image failed");
		  return(ERROR);
		}
		adornImage(outPtr, mskPtr, nx, ny, nz, (int)nxRange[0], (int)nyRange[0],
		  (int)nzRange[0], regDataPtr->directives.clipping, regDataPtr->backgrnd);
		computeMskSnr(&mse, &snr, inPtr2, outPtr, mskPtr2, mskPtr,
		  regDataPtr->directives.maskCombine, nx, ny, nz);
		if (exportFit(&fit) == ERROR ) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Assignment of transformation variables failed");
		  return(ERROR);
		}
		if (exportSummary(&fit, mse, snr) == ERROR ) {
		  freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		  message("ERROR - Assignment of summary variables failed");
		  return(ERROR);
		}
		if (regDataPtr->directives.exportFit)
		  if (fExportFit(&fit, regDataPtr->outFit) == ERROR ) {
			freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
			message("ERROR - Exportation of transformation variables failed");
			return(ERROR);
		  }
		freePyramids(pyrHdl1, pyrHdl2, mskHdl1, mskHdl2, regDataPtr->levels);
		return(!ERROR);
} /* End of regFloat3d */
