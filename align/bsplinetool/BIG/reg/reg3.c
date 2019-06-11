#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<float.h>
#include		<math.h>
#include		"phil.h"

#include		"register.h"
#include		"reg3.h"

#include		"quant.h"
#include		"reg0.h"
#include		"svdcmp.h"

/************************************************************************/
/* FUNCTION: marquardt													*/
/************************************************************************/
int					marquardt		(struct	rOper		*directives,
									 struct	fitRec		*fit,
									 double				*chisq,
									 double				*lambda,
									 float				*inPtr1,
									 float				*inPtr2,
									 double				*splinePtr,
									 float				*mskPtr1,
									 float				*mskPtr2,
									 float				*mskPtr,
									 float				*outPtr,
									 float				*gradHdl[],
									 int				*successful,
									 int				ia[],
									 double				covar[][hessianSize],
									 double				hessian[][hessianSize],
									 double				beta[],
									 double				firstLambda,
									 double				lambdaScale,
									 double				epsilon,
									 double				minGain,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				u[hessianSize][hessianSize], v[hessianSize][hessianSize];
		double				da[hessianSize], w[hessianSize];
		struct	fitRec		fitTry, deltaFit;
		double				oChisq;
		double				wmax, thresh;
		double				more;
		int					i, j, k;

		oChisq = *chisq;
		for (i = 0; (i < hessianSize); i++) {
		  for (j = 0; (j < hessianSize); j++)
			if (ia[i] && ia[j])
			  u[i][j] = hessian[i][j];
			else
			  u[i][j] = 0.0;
		  if (ia[i])
			u[i][i] *= 1.0 + *lambda;
		}
		if (svdcmp(&(u[0][0]), (long)hessianSize, (long)hessianSize, w, &(v[0][0])) == ERROR) {
		  message("ERROR - Unable to perform svdcmp in marquardt");
		  return(ERROR);
		}
		wmax = 0.0;
		for (k = 0; (k < hessianSize); k++)
		  if (w[k] > wmax)
			wmax = w[k];
		thresh = epsilon * wmax;
		for (k = 0; (k < hessianSize); k++) {
		  if (w[k] < thresh) {
			w[k] = 0.0;
			for (i = 0; (i < hessianSize); i++) {
			  u[i][k] = 0.0;
			  v[i][k] = 0.0;
			}
		  }
		  if (!ia[k])
			for (i = 0; (i < hessianSize); i++)
			  v[k][i] = 0.0;
		}
		for (i = 0; (i < hessianSize); i++)
		  for (j = i; (j < hessianSize); j++)
			for (k = 0, covar[i][j] = 0.0; (k < hessianSize); k++)
			  if (w[k] != 0.0)
				covar[i][j] += v[i][k] * v[j][k] / (w[k] * w[k]);
		for (i = 1; (i < hessianSize); i++)
		  for (j = 0; (j < i); j++)
			covar[i][j] = covar[j][i];
		if (svbksb(&(u[0][0]), w, &(v[0][0]), (long)hessianSize, (long)hessianSize, beta, da)
		  == ERROR) {
		  message("ERROR - Unable to perform svbksb in marquardt");
		  return(ERROR);
		}
		more = 0.0;
		k = 0;
		for (i = 0; (i < 3); k++, i++)
		  if (ia[k]) {
			if ((fit->dx[i] * fit->dx[i]) > epsilon)
			  more += da[i] * da[i] / (fit->dx[i] * fit->dx[i]);
			else if ((beta[i] * beta[i]) > epsilon)
			  more += (da[i] * da[i]) * (hessian[i][i] * hessian[i][i]) / (beta[i] * beta[i]);
			deltaFit.dx[i] = da[i];
		  }
		  else
			deltaFit.dx[i] = 0.0;
		if (directives->xSkew || directives->ySkew || directives->zSkew) {
		  for (i = 0; (i < 3); i++)
			for (j = 0; (j < 3); k++, j++)
			  if (ia[k]) {
				if ((fit->skew[i][j] * fit->skew[i][j]) > epsilon)
				  more += da[k] * da[k] / (fit->skew[i][j] * fit->skew[i][j]);
				else if ((beta[k] * beta[k]) > epsilon)
				  more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
				deltaFit.skew[i][j] = da[k];
			  }
			  else
				deltaFit.skew[i][j] = 0.0;
		}
		else if (directives->xRot || directives->yRot || directives->zRot
		  || directives->isoScaling) {
		  if (ia[k]) {
			if ((fit->phi * fit->phi) > epsilon)
			  more += da[k] * da[k] / (fit->phi * fit->phi);
			else if ((beta[k] * beta[k]) > epsilon)
			  more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
			deltaFit.phi = da[k];
		  }
		  else
			deltaFit.phi = 0.0;
		  k++;
		  if (ia[k]) {
			if ((fit->theta * fit->theta) > epsilon)
			  more += da[k] * da[k] / (fit->theta * fit->theta);
			else if ((beta[k] * beta[k]) > epsilon)
			  more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
			deltaFit.theta = da[k];
		  }
		  else
			deltaFit.theta = 0.0;
		  k++;
		  if (ia[k]) {
			if ((fit->psi * fit->psi) > epsilon)
			  more += da[k] * da[k] / (fit->psi * fit->psi);
			else if ((beta[k] * beta[k]) > epsilon)
			  more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
			deltaFit.psi = da[k];
		  }
		  else
			deltaFit.psi = 0.0;
		  k++;
		  if (ia[k]) {
			if ((fit->lambda * fit->lambda) > epsilon)
			  more += (da[k] * da[k]) / (fit->lambda * fit->lambda);
			else if ((beta[k] * beta[k]) > epsilon)
			  more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
			deltaFit.lambda = da[k];
		  }
		  else
			deltaFit.lambda = 0.0;
		  k++;
		}
		else
		  k = 12;
		if (ia[k]) {
		  if ((fit->gamma * fit->gamma) > epsilon)
			more += da[k] * da[k] / (fit->gamma * fit->gamma);
		  else if ((beta[k] * beta[k]) > epsilon)
			more += (da[k] * da[k]) * (hessian[k][k] * hessian[k][k]) / (beta[k] * beta[k]);
		  deltaFit.gamma = da[k];
		}
		else
		  deltaFit.gamma = 0.0;
		fitTry = *fit;
		if (directives->xSkew || directives->ySkew || directives->zSkew) {
		  for (i = 0; (i < 3); i++)
			deltaFit.skew[i][i] += 1.0;
		  for (i = 0; (i < 3); i++)
			for (j = 0; (j < 3); j++)
			  for (k = 0, fitTry.skew[i][j] = 0.0; (k < 3); k++)
				fitTry.skew[i][j] += fit->skew[i][k] * deltaFit.skew[k][j];
		  k = 3;
		  for (i = 0; (i < 3); i++)
			for (j = 0; (j < 3); k++, j++)
			  if (!ia[k])
				fitTry.skew[i][j] = fit->skew[i][j];
		  if (ia[k])
			fitTry.gamma += deltaFit.gamma;
		}
		else if (directives->xRot || directives->yRot || directives->zRot
		  || directives->isoScaling) {
		  k = 4;
		  if (ia[k])
			fitTry.theta = asin((cos(deltaFit.phi) * sin(fit->theta) + sin(deltaFit.phi)
			  * cos(fit->theta) * sin(fit->psi)) * cos(deltaFit.theta) + sin(deltaFit.theta)
			  * cos(fit->theta) * cos(fit->psi));
		  k = 3;
		  if (ia[k])
			fitTry.phi = asin((cos(deltaFit.theta) * (cos(deltaFit.phi) * sin(fit->phi)
			  * cos(fit->theta) + (cos(fit->phi) * cos(fit->psi) - sin(fit->phi) * sin(fit->theta)
			  * sin(fit->psi)) * sin(deltaFit.phi)) - (cos(fit->phi) * sin(fit->psi) + sin(fit->phi)
			  * sin(fit->theta) * cos(fit->psi)) * sin(deltaFit.theta)) / cos(fitTry.theta));
		  k = 5;
		  if (ia[k])
			fitTry.psi = asin((sin(deltaFit.psi) * (cos(deltaFit.theta) * cos(fit->theta)
			  * cos(fit->psi) - (cos(deltaFit.phi) * sin(fit->theta) + sin(deltaFit.phi)
			  * cos(fit->theta) * sin(fit->psi)) * sin(deltaFit.theta)) - (sin(deltaFit.phi)
			  * sin(fit->theta) - cos(deltaFit.phi) * cos(fit->theta) * sin(fit->psi))
			  * cos(deltaFit.psi)) / cos(fitTry.theta));
		  k++;
		  if (ia[k])
			fitTry.lambda += deltaFit.lambda;
		  k++;
		  if (ia[k])
			fitTry.gamma += deltaFit.gamma;
		  eulerRotation(&deltaFit);
		  eulerRotation(fit);
		  eulerRotation(&fitTry);
		  for (i = 0; (i < 3); i++)
			for (j = 0; (j < 3); j++) {
			  fit->skew[i][j] *= exp(fit->lambda);
			  fitTry.skew[i][j] *= exp(fitTry.lambda);
			}
		}
		else if (ia[12])
		  fitTry.gamma += deltaFit.gamma;
		for (i = 0; (i < 3); i++)
		  for (j = 0; (j < 3); j++)
			fitTry.dx[i] += fitTry.skew[i][j] * deltaFit.dx[j];
		wmax = 0.0;
		for (i = 0; (i < hessianSize); i++)
		  if (ia[i])
			wmax += 1.0;
		if (((more / (wmax * wmax)) > epsilon) || (*lambda == 0.0)) {
		  if (invTransform(&fitTry, inPtr2, splinePtr, outPtr, mskPtr2, mskPtr, nx, ny, nz,
			directives->interpolation) == ERROR) {
			message("ERROR - Unable to perform inverse transformation");
			return(ERROR);
		  }
		  combineMasks(mskPtr1, mskPtr, directives->maskCombine, nx, ny, nz);
		  if (directives->xSkew || directives->ySkew || directives->zSkew) {
			computeSkewChisq(da, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			  chisq, nx, ny, nz, fitTry.origin[0], fitTry.origin[1], fitTry.origin[2]);
		  }
		  else if (directives->xRot || directives->yRot || directives->zRot
			|| directives->isoScaling) {
			computeRotChisq(da, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			  chisq, nx, ny, nz, fitTry.origin[0], fitTry.origin[1], fitTry.origin[2]);
		  }
		  else {
			computeSkewChisq(da, inPtr1, mskPtr, outPtr, gradHdl[0], gradHdl[1], gradHdl[2],
			  chisq, nx, ny, nz, fitTry.origin[0], fitTry.origin[1], fitTry.origin[2]);
		  }
		  *chisq *= exp(2.0 * fitTry.gamma) / absDeterminant(fitTry.skew);
		  if (*chisq < oChisq) {
			*successful = TRUE;
			*lambda /= lambdaScale;
			*fit = fitTry;
			for (i = 0; (i < hessianSize); i++)
			  beta[i] = da[i];
			if (*chisq > ((1.0 - minGain) * oChisq))
			  *chisq = -*chisq;
		  }
		  else {
			*successful = FALSE;
			if (*lambda < firstLambda)
			  *lambda = firstLambda;
			else
			  if (*lambda < (double)FLT_MAX)
				*lambda *= lambdaScale;
			*chisq = oChisq;
		  }
		} else *chisq = -oChisq;
		return(!ERROR);
} /* End of marquardt */
