#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<string.h>
#include		<float.h>
#include		<math.h>
#include		"phil.h"

#include		"register.h"
#include		"reg2.h"

#include		"BsplnTrf.h"
#include		"convolve.h"
#include		"getPut.h"
#include		"quant.h"
#include		"reg0.h"
#include		"reg3.h"

/************************************************************************/
static	void		computeRotHessian
									(double				hessian[][hessianSize],
									 float				*inPtr,
									 float				*gradxPtr,
									 float				*gradyPtr,
									 float				*gradzPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				xO,
									 double				yO,
									 double				zO);
static	void		freeTmp			(float				*gradHdl[],
									 double				*splinePtr);

/************************************************************************/
static	void		computeRotHessian
									(double				hessian[][hessianSize],
									 float				*inPtr,
									 float				*gradxPtr,
									 float				*gradyPtr,
									 float				*gradzPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				xO,
									 double				yO,
									 double				zO) {

		double				agc;
		double				gxx, gxy, gxz;
		double				gyy, gyz;
		double				gzz;
		double				xgy, xgz;
		double				ygx, ygz;
		double				zgx, zgy;
		double				xgxygyzgz;
		double				xgxx, xgxy, xgxz, xgyy, xgyz, xgzz;
		double				ygxx, ygxy, ygxz, ygyy, ygyz, ygzz;
		double				zgxx, zgxy, zgxz, zgyy, zgyz, zgzz;
		double				normalization = 0.0;
		double				x, y, z;
		int					i, j;

		for (i = 0; (i < hessianSize); i++)
		  for (j = 0; (j < hessianSize); j++)
			hessian[i][j] = 0.0;
		for (z = -zO; (z < ((double)nz - zO)); z += 1.0) {
		  for (y = -yO; (y < ((double)ny - yO)); y += 1.0) {
			for (x = -xO; (x < ((double)nx - xO)); inPtr++, mskPtr++,
			  gradxPtr++, gradyPtr++, gradzPtr++, x += 1.0)
			  if ((int)*mskPtr) {
				agc = (double)*inPtr;
				gxx = (double)*gradxPtr * (double)*gradxPtr;
				gxy = (double)*gradxPtr * (double)*gradyPtr;
				gxz = (double)*gradxPtr * (double)*gradzPtr;
				gyy = (double)*gradyPtr * (double)*gradyPtr;
				gyz = (double)*gradyPtr * (double)*gradzPtr;
				gzz = (double)*gradzPtr * (double)*gradzPtr;
				xgy = x * (double)*gradyPtr;
				xgz = x * (double)*gradzPtr;
				ygx = y * (double)*gradxPtr;
				ygz = y * (double)*gradzPtr;
				zgx = z * (double)*gradxPtr;
				zgy = z * (double)*gradyPtr;
				xgxygyzgz = x * (double)*gradxPtr + y * (double)*gradyPtr + z * (double)*gradzPtr;
				xgxx = x * gxx;
				xgxy = x * gxy;
				xgxz = x * gxz;
				xgyy = x * gyy;
				xgyz = x * gyz;
				xgzz = x * gzz;
				ygxx = y * gxx;
				ygxy = y * gxy;
				ygxz = y * gxz;
				ygyy = y * gyy;
				ygyz = y * gyz;
				ygzz = y * gzz;
				zgxx = z * gxx;
				zgxy = z * gxy;
				zgxz = z * gxz;
				zgyy = z * gyy;
				zgyz = z * gyz;
				zgzz = z * gzz;
				hessian[0][0] += gxx;
				hessian[0][1] += gxy;
				hessian[0][2] += gxz;
				hessian[0][3] += zgxy - ygxz;
				hessian[0][4] += -zgxx + xgxz;
				hessian[0][5] += ygxx - xgxy;
				hessian[0][6] += -xgxx - ygxy - zgxz;
				hessian[0][7] -= (double)*gradxPtr * agc;
				hessian[1][1] += gyy;
				hessian[1][2] += gyz;
				hessian[1][3] += zgyy - ygyz;
				hessian[1][4] += -zgxy + xgyz;
				hessian[1][5] += ygxy - xgyy;
				hessian[1][6] += -xgxy - ygyy - zgyz;
				hessian[1][7] -= (double)*gradyPtr * agc;
				hessian[2][2] += gzz;
				hessian[2][3] += zgyz - ygzz;
				hessian[2][4] += -zgxz + xgzz;
				hessian[2][5] += ygxz - xgyz;
				hessian[2][6] += -xgxz - ygyz - zgzz;
				hessian[2][7] -= (double)*gradzPtr * agc;
				hessian[3][3] += (-zgy + ygz) * (-zgy + ygz);
				hessian[3][4] += (-zgy + ygz) * (zgx - xgz);
				hessian[3][5] += (-zgy + ygz) * (-ygx + xgy);
				hessian[3][6] += (-zgy + ygz) * xgxygyzgz;
				hessian[3][7] += (-zgy + ygz) * agc;
				hessian[4][4] += (zgx - xgz) * (zgx - xgz);
				hessian[4][5] += (zgx - xgz) * (-ygx + xgy);
				hessian[4][6] += (zgx - xgz) * xgxygyzgz;
				hessian[4][7] += (zgx - xgz) * agc;
				hessian[5][5] += (-ygx + xgy) * (-ygx + xgy);
				hessian[5][6] += (-ygx + xgy) * xgxygyzgz;
				hessian[5][7] += (-ygx + xgy) * agc;
				hessian[6][6] += xgxygyzgz * xgxygyzgz;
				hessian[6][7] += xgxygyzgz * agc;
				hessian[7][7] += agc * agc;
				normalization += 1.0;
			  }
		  }
		}
		if (normalization != 0.0)
		  for (i = 0; (i < hessianSize); i++) {
			for (j = 0; (j < i); j++) {
			  hessian[j][i] /= normalization;
			  hessian[i][j] = hessian[j][i];
			}
			hessian[i][i] /= normalization;
		  }
} /* End of computeRotHessian */

/************************************************************************/
static	void		freeTmp			(float				*gradHdl[],
									 double				*splinePtr) {

		if (gradHdl[2] != (float *)NULL)
		  free(gradHdl[2]);
		gradHdl[2] = (float *)NULL;
		if (gradHdl[1] != (float *)NULL)
		  free(gradHdl[1]);
		gradHdl[1] = (float *)NULL;
		if (gradHdl[0] != (float *)NULL)
		  free(gradHdl[0]);
		gradHdl[0] = (float *)NULL;
		if (splinePtr != (double *)NULL)
		  free(splinePtr);
} /* End of freeTmp */

/************************************************************************/
/* FUNCTION: computeSkewHessian											*/
/************************************************************************/
void				computeSkewHessian
									(double				hessian[][hessianSize],
									 float				*inPtr,
									 float				*gradxPtr,
									 float				*gradyPtr,
									 float				*gradzPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				xO,
									 double				yO,
									 double				zO) {

		double				agc;
		double				xx, xy, xz;
		double				yy, yz;
		double				zz;
		double				gxx, gxy, gxz;
		double				gyy, gyz;
		double				gzz;
		double				normalization = 0.0;
		double				x, y, z;
		int					i, j;

		for (i = 0; (i < hessianSize); i++)
		  for (j = 0; (j < hessianSize); j++)
			hessian[i][j] = 0.0;
		for (z = -zO; (z < ((double)nz - zO)); z += 1.0) {
		  zz = z * z;
		  for (y = -yO; (y < ((double)ny - yO)); y += 1.0) {
			yy = y * y;
			yz = y * z;
			for (x = -xO; (x < ((double)nx - xO)); inPtr++, mskPtr++,
			  gradxPtr++, gradyPtr++, gradzPtr++, x += 1.0)
			  if ((int)*mskPtr) {
				agc = (double)*inPtr;
				xx = x * x;
				xy = x * y;
				xz = x * z;
				gxx = (double)*gradxPtr * (double)*gradxPtr;
				gxy = (double)*gradxPtr * (double)*gradyPtr;
				gxz = (double)*gradxPtr * (double)*gradzPtr;
				gyy = (double)*gradyPtr * (double)*gradyPtr;
				gyz = (double)*gradyPtr * (double)*gradzPtr;
				gzz = (double)*gradzPtr * (double)*gradzPtr;
				hessian[0][0] += gxx;
				hessian[0][1] += gxy;
				hessian[0][2] += gxz;
				hessian[0][3] -= x * gxx;
				hessian[0][4] -= y * gxx;
				hessian[0][5] -= z * gxx;
				hessian[0][6] -= x * gxy;
				hessian[0][7] -= y * gxy;
				hessian[0][8] -= z * gxy;
				hessian[0][9] -= x * gxz;
				hessian[0][10] -= y * gxz;
				hessian[0][11] -= z * gxz;
				hessian[0][12] -= (double)*gradxPtr * agc;
				hessian[1][1] += gyy;
				hessian[1][2] += gyz;
				hessian[1][6] -= x * gyy;
				hessian[1][7] -= y * gyy;
				hessian[1][8] -= z * gyy;
				hessian[1][9] -= x * gyz;
				hessian[1][10] -= y * gyz;
				hessian[1][11] -= z * gyz;
				hessian[1][12] -= (double)*gradyPtr * agc;
				hessian[2][2] += gzz;
				hessian[2][9] -= x * gzz;
				hessian[2][10] -= y * gzz;
				hessian[2][11] -= z * gzz;
				hessian[2][12] -= (double)*gradzPtr * agc;
				hessian[3][3] += xx * gxx;
				hessian[3][4] += xy * gxx;
				hessian[3][5] += xz * gxx;
				hessian[3][6] += xx * gxy;
				hessian[3][7] += xy * gxy;
				hessian[3][8] += xz * gxy;
				hessian[3][9] += xx * gxz;
				hessian[3][10] += xy * gxz;
				hessian[3][11] += xz * gxz;
				hessian[3][12] += x * (double)*gradxPtr * agc;
				hessian[4][4] += yy * gxx;
				hessian[4][5] += yz * gxx;
				hessian[4][7] += yy * gxy;
				hessian[4][8] += yz * gxy;
				hessian[4][10] += yy * gxz;
				hessian[4][11] += yz * gxz;
				hessian[4][12] += y * (double)*gradxPtr * agc;
				hessian[5][5] += zz * gxx;
				hessian[5][8] += zz * gxy;
				hessian[5][11] += zz * gxz;
				hessian[5][12] += z * (double)*gradxPtr * agc;
				hessian[6][6] += xx * gyy;
				hessian[6][7] += xy * gyy;
				hessian[6][8] += xz * gyy;
				hessian[6][9] += xx * gyz;
				hessian[6][10] += xy * gyz;
				hessian[6][11] += xz * gyz;
				hessian[6][12] += x * (double)*gradyPtr * agc;
				hessian[7][7] += yy * gyy;
				hessian[7][8] += yz * gyy;
				hessian[7][10] += yy * gyz;
				hessian[7][11] += yz * gyz;
				hessian[7][12] += y * (double)*gradyPtr * agc;
				hessian[8][8] += zz * gyy;
				hessian[8][11] += zz * gyz;
				hessian[8][12] += z * (double)*gradyPtr * agc;
				hessian[9][9] += xx * gzz;
				hessian[9][10] += xy * gzz;
				hessian[9][11] += xz * gzz;
				hessian[9][12] += x * (double)*gradzPtr * agc;
				hessian[10][10] += yy * gzz;
				hessian[10][11] += yz * gzz;
				hessian[10][12] += y * (double)*gradzPtr * agc;
				hessian[11][11] += zz * gzz;
				hessian[11][12] += z * (double)*gradzPtr * agc;
				hessian[12][12] += agc * agc;
				normalization += 1.0;
			  }
		  }
		}
		hessian[1][3] = hessian[0][6];
		hessian[1][4] = hessian[0][7];
		hessian[1][5] = hessian[0][8];
		hessian[2][3] = hessian[0][9];
		hessian[2][4] = hessian[0][10];
		hessian[2][5] = hessian[0][11];
		hessian[2][6] = hessian[1][9];
		hessian[2][7] = hessian[1][10];
		hessian[2][8] = hessian[1][11];
		hessian[4][6] = hessian[3][7];
		hessian[5][6] = hessian[3][8];
		hessian[4][9] = hessian[3][10];
		hessian[5][9] = hessian[3][11];
		hessian[5][7] = hessian[4][8];
		hessian[5][10] = hessian[4][11];
		hessian[7][9] = hessian[6][10];
		hessian[8][9] = hessian[6][11];
		hessian[8][10] = hessian[7][11];
		if (normalization != 0.0)
		  for (i = 0; (i < hessianSize); i++) {
			for (j = 0; (j < i); j++) {
			  hessian[j][i] /= normalization;
			  hessian[i][j] = hessian[j][i];
			}
			hessian[i][i] /= normalization;
		  }
} /* End of computeSkewHessian */

/************************************************************************/
/* FUNCTION: optimize													*/
/************************************************************************/
int					optimize		(struct fitRec		*fit,
									 struct	rOper		*directives,
									 int				ia[],
									 float				*inPtr1,
									 float				*inPtr2,
									 float				*mskPtr1,
									 float				*mskPtr2,
									 float				*outPtr,
									 float				*mskPtr,
									 double				*lambda,
									 double				firstLambda,
									 double				lambdaScale,
									 double				epsilon,
									 double				minGain,
									 int				nx,
									 int				ny,
									 int				nz) {

		float				*gradHdl[3];
		double				*splinePtr, *p;
		double				mse = 0.0;
		double				hessian[hessianSize][hessianSize];
		double				covar[hessianSize][hessianSize];
		double				beta[hessianSize];
		double				chisq, fChisq;
		int					successful;
		int					i, j;

		splinePtr = (double *)NULL;
		gradHdl[0] = gradHdl[1] = gradHdl[2] = (float *)NULL;
		splinePtr = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
		if (splinePtr == (double *)NULL) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Not enough memory for B-spline coefficients");
		  return(ERROR);
		}
		p = splinePtr;
		for (j = 0; (j < nz); j++)
		  for (i = 0; (i < ny); p += (ptrdiff_t)nx, i++)
			getxF2D(inPtr2, 0L, (long)i, (long)j, (long)nx, (long)ny, p, (long)nx);
		switch (directives->interpolation) {
		  case one:
			if (directBsplineMirror(splinePtr, splinePtr, (long)nx, (long)ny, (long)nz, 1)
			  == ERROR) {
			  freeTmp(gradHdl, splinePtr);
			  message("ERROR - B(1) spline computation failed");
			  return(ERROR);
			}
			break;
		  case three:
			if (directBsplineMirror(splinePtr, splinePtr, (long)nx, (long)ny, (long)nz, 3)
			  == ERROR) {
			  freeTmp(gradHdl, splinePtr);
			  message("ERROR - B(3) spline computation failed");
			  return(ERROR);
			}
			break;
		  default:
			message("ERROR - Invalid interpolation degree");
			return(ERROR);
		}
		gradHdl[0] = (float *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(float));
		if (gradHdl[0] == (float *)NULL) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Not enough memory for x gradient coefficients");
		  return(ERROR);
		}
		gradHdl[1] = (float *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(float));
		if (gradHdl[1] == (float *)NULL) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Not enough memory for y gradient coefficients");
		  return(ERROR);
		}
		gradHdl[2] = (float *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(float));
		if (gradHdl[2] == (float *)NULL) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Not enough memory for z gradient coefficients");
		  return(ERROR);
		}
		if (spatialGradient(gradHdl, inPtr1, nx, ny, nz, directives->interpolation) == ERROR) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Unable to perform spatial gradient computation");
		  return(ERROR);
		}
		if (directives->xRot || directives->yRot || directives->zRot || directives->isoScaling) {
		  eulerRotation(fit);
		  for (i = 0; (i < 3); i++)
			for (j = 0; (j < 3); j++)
			  fit->skew[i][j] *= exp(fit->lambda);
		}
		if (invTransform(fit, inPtr2, splinePtr, outPtr, mskPtr2, mskPtr, nx, ny, nz,
		  directives->interpolation) == ERROR) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Unable to perform inverse transformation");
		  return(ERROR);
		}
		combineMasks(mskPtr1, mskPtr, directives->maskCombine, nx, ny, nz);
		if (directives->xSkew || directives->ySkew || directives->zSkew) {
		  computeSkewHessian(hessian, inPtr1, gradHdl[0], gradHdl[1], gradHdl[2], mskPtr,
			nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		  computeSkewChisq(beta, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			&fChisq, nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		}
		else if (directives->xRot || directives->yRot || directives->zRot
		  || directives->isoScaling) {
		  computeRotHessian(hessian, inPtr1, gradHdl[0], gradHdl[1], gradHdl[2], mskPtr,
			nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		  computeRotChisq(beta, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			&fChisq, nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		}
		else {
		  computeSkewHessian(hessian, inPtr1, gradHdl[0], gradHdl[1], gradHdl[2], mskPtr,
			nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		  computeSkewChisq(beta, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			&fChisq, nx, ny, nz, fit->origin[0], fit->origin[1], fit->origin[2]);
		}
		fChisq *= exp(2.0 * fit->gamma) / absDeterminant(fit->skew);
		chisq = fChisq;
		successful = TRUE;
		do {
		  if (marquardt(directives, fit, &chisq, lambda, inPtr1, inPtr2, splinePtr,
			mskPtr1, mskPtr2, mskPtr, outPtr, gradHdl, &successful, ia, covar, hessian, beta,
			firstLambda, lambdaScale, epsilon, minGain, nx, ny, nz) == ERROR) {
			freeTmp(gradHdl, splinePtr);
			message("ERROR - Optimization aborted");
			return(ERROR);
		  }
		} while (chisq > 0.0);
		chisq = fabs(chisq);
		*lambda = 0.0;
		if (marquardt(directives, fit, &chisq, lambda, inPtr1, inPtr2, splinePtr,
		  mskPtr1, mskPtr2, mskPtr, outPtr, gradHdl, &successful, ia, covar, hessian, beta,
		  firstLambda, lambdaScale, epsilon, minGain, nx, ny, nz) == ERROR) {
		  freeTmp(gradHdl, splinePtr);
		  message("ERROR - Optimization aborted");
		  return(ERROR);
		}
		freeTmp(gradHdl, splinePtr);
		return(!ERROR);
} /* End of optimize */

/************************************************************************/
/* FUNCTION: spatialGradient											*/
/************************************************************************/
int					spatialGradient	(float				*gradHdl[],
									 float				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 enum	iDegree		interpolation) {

		double				*splinePtr;
		double				*inData, *outData;
		double				*p;
		double				*gx, *gy, *gz;
		double				derivative[3], reconstruction[3];
		long				dLength, rLength;
		long				i, j, k;

		splinePtr = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
		if (splinePtr == (double *)NULL) {
		  message("ERROR - Not enough memory for B-spline gradient coefficients");
		  return(ERROR);
		}
		gx = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
		if (gx == (double *)NULL) {
		  free(splinePtr);
		  message("ERROR - Not enough memory for temporary x gradient");
		  return(ERROR);
		}
		gy = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
		if (gy == (double *)NULL) {
		  free(gx);
		  free(splinePtr);
		  message("ERROR - Not enough memory for temporary y gradient");
		  return(ERROR);
		}
		gz = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
		if (gz == (double *)NULL) {
		  free(gy);
		  free(gx);
		  free(splinePtr);
		  message("ERROR - Not enough memory for temporary z gradient");
		  return(ERROR);
		}
		p = splinePtr;
		for (j = 0L; (j < (long)nz); j++)
		  for (i = 0L; (i < (long)ny); p += (ptrdiff_t)nx, i++)
			getxF2D(inPtr, 0L, i, j, (long)nx, (long)ny, p, (long)nx);
		switch (interpolation) {
		  case one:
			dLength = 3L;
			derivative[0] = 0.5;
			derivative[1] = 0.0;
			derivative[2] = -0.5;
			rLength = 1L;
			reconstruction[0] = 1.0;
			if (directBsplineMirror(splinePtr, splinePtr, (long)nx, (long)ny, (long)nz, 1)
			  == ERROR) {
			  free(gz);
			  free(gy);
			  free(gx);
			  free(splinePtr);
			  message("ERROR - B(1) spline computation for spatial gradient failed");
			  return(ERROR);
			}
			break;
		  case three:
			dLength = 3L;
			derivative[0] = 0.5;
			derivative[1] = 0.0;
			derivative[2] = -0.5;
			rLength = 3L;
			reconstruction[0] = 1.0 / 6.0;
			reconstruction[1] = 4.0 / 6.0;
			reconstruction[2] = 1.0 / 6.0;
			if (directBsplineMirror(splinePtr, splinePtr, (long)nx, (long)ny, (long)nz, 3)
			  == ERROR) {
			  free(gz);
			  free(gy);
			  free(gx);
			  free(splinePtr);
			  message("ERROR - B(3) spline computation for spatial gradient failed");
			  return(ERROR);
			}
			break;
		  default:
			message("ERROR - Invalid interpolation degree");
			return(ERROR);
		}
		p = splinePtr;
		for (j = 0L; (j < (long)nz); j++)
		  for (i = 0L; (i < (long)ny); p += (ptrdiff_t)nx,
			gx += (ptrdiff_t)nx, gy += (ptrdiff_t)nx, gz += (ptrdiff_t)nx, i++) {
			firConvolveMirror(p, gx, (long)nx, 0L, derivative, dLength);
			if (ny > 1)
			  firConvolveMirror(p, gy, (long)nx, 0L, reconstruction, rLength);
			if (nz > 1)
			  firConvolveMirror(p, gz, (long)nx, 0L, reconstruction, rLength);
		  }
		free(splinePtr);
		gx -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		gy -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		gz -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		if (ny > 1) {
		  inData = (double *)malloc((size_t)ny * sizeof(double));
		  if (inData == (double *)NULL) {
			free(gz);
			free(gy);
			free(gx);
			message("ERROR - Not enough memory for (y) inData");
			return(ERROR);
		  }
		  outData = (double *)malloc((size_t)ny * sizeof(double));
		  if (outData == (double *)NULL) {
			free(inData);
			free(gz);
			free(gy);
			free(gx);
			message("ERROR - Not enough memory for (y) outData");
			return(ERROR);
		  }
		  for (j = 0L; (j < (long)nz); j++)
			for (i = 0L; (i < (long)nx); i++) {
			  getyD2D(gx, i, 0L, j, (long)nx, (long)ny, inData, (long)ny);
			  firConvolveMirror(inData, outData, (long)ny, 0L, reconstruction, rLength);
			  putyD2D(gx, i, 0L, j, (long)nx, (long)ny, outData, (long)ny);
			  getyD2D(gy, i, 0L, j, (long)nx, (long)ny, inData, (long)ny);
			  firConvolveMirror(inData, outData, (long)ny, 0L, derivative, dLength);
			  putyD2D(gy, i, 0L, j, (long)nx, (long)ny, outData, (long)ny);
			  if (nz > 1) {
				getyD2D(gz, i, 0L, j, (long)nx, (long)ny, inData, (long)ny);
				firConvolveMirror(inData, outData, (long)ny, 0L, reconstruction, rLength);
				putyD2D(gz, i, 0L, j, (long)nx, (long)ny, outData, (long)ny);
			  }
			}
		  free(outData);
		  free(inData);
		}
		else {
		  p = gy;
		  for (k = 0L; (k < (long)nz); k++)
			for (j = 0L; (j < (long)ny); j++)
			  for (i = 0L; (i < (long)nx); i++)
				*p++ = 0.0;
		}
		if (nz > 1) {
		  inData = (double *)malloc((size_t)nz * sizeof(double));
		  if (inData == (double *)NULL) {
			free(gz);
			free(gy);
			free(gx);
			message("ERROR - Not enough memory for (z) inData");
			return(ERROR);
		  }
		  outData = (double *)malloc((size_t)nz * sizeof(double));
		  if (outData == (double *)NULL) {
			free(inData);
			free(gz);
			free(gy);
			free(gx);
			message("ERROR - Not enough memory for (z) outData");
			return(ERROR);
		  }
		  for (j = 0L; (j < (long)ny); j++)
			for (i = 0L; (i < (long)nx); i++) {
			  getzD2D(gx, i, j, 0L, (long)nx, (long)ny, inData, (long)nz);
			  firConvolveMirror(inData, outData, (long)nz, 0L, reconstruction, rLength);
			  putzD2D(gx, i, j, 0L, (long)nx, (long)ny, outData, (long)nz);
			  getzD2D(gy, i, j, 0L, (long)nx, (long)ny, inData, (long)nz);
			  firConvolveMirror(inData, outData, (long)nz, 0L, reconstruction, rLength);
			  putzD2D(gy, i, j, 0L, (long)nx, (long)ny, outData, (long)nz);
			  getzD2D(gz, i, j, 0L, (long)nx, (long)ny, inData, (long)nz);
			  firConvolveMirror(inData, outData, (long)nz, 0L, derivative, dLength);
			  putzD2D(gz, i, j, 0L, (long)nx, (long)ny, outData, (long)nz);
			}
		  free(outData);
		  free(inData);
		}
		else {
		  p = gz;
		  for (k = 0L; (k < (long)nz); k++)
			for (j = 0L; (j < (long)ny); j++)
			  for (i = 0L; (i < (long)nx); i++)
				*p++ = 0.0;
		}
		for (j = 0L; (j < (long)nz); j++)
		  for (i = 0L; (i < (long)ny); gx += (ptrdiff_t)nx, gy += (ptrdiff_t)nx,
			gz += (ptrdiff_t)nx, i++) {
			putxD2F(gradHdl[0], 0L, i, j, (long)nx, (long)ny, gx, (long)nx);
			putxD2F(gradHdl[1], 0L, i, j, (long)nx, (long)ny, gy, (long)nx);
			putxD2F(gradHdl[2], 0L, i, j, (long)nx, (long)ny, gz, (long)nx);
		  }
		gx -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		gy -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		gz -= (ptrdiff_t)((long)nx * (long)ny * (long)nz);
		free(gz);
		free(gy);
		free(gx);
		return(!ERROR);
} /* End of spatialGradient */
