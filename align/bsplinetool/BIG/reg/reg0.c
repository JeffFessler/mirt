#include		<stddef.h>
#include		<stdlib.h>
#include		<stdio.h>
#include		<math.h>
#include		<limits.h>
#include		<float.h>
#include		"phil.h"

#include		"quant.h"
#include		"register.h"
#include		"reg0.h"

#include		"BsplnTrf.h"
#include		"BsplnWgt.h"
#include		"getPut.h"
#include		"reg1.h"

/************************************************************************/
static	int			anisoGauss3d	(float				*inPtr,
									 float				*outPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				sx,
									 double				sy,
									 double				sz);
static	int			gaussj			(double				a[][3],
									 double				b[][3]);

/************************************************************************/
static	int			anisoGauss3d	(float				*inPtr,
									 float				*outPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				sx,
									 double				sy,
									 double				sz) {

		double				*data;
		double				xKernel[3], yKernel[3], zKernel[3];
		double				xAlpha, xLambda, xC0;
		double				yAlpha, yLambda, yC0;
		double				zAlpha, zLambda, zC0;
		long				x, y, z;
		int					i;

		xLambda = (sx * sx) / 6.0;
		yLambda = (sy * sy) / 6.0;
		zLambda = (sz * sz) / 6.0;
		xAlpha = (2.0 + (1.0 / xLambda) - sqrt(pow(2.0 + 1.0 / xLambda, 2.0) - 4.0)) * 0.5;
		yAlpha = (2.0 + (1.0 / yLambda) - sqrt(pow(2.0 + 1.0 / yLambda, 2.0) - 4.0)) * 0.5;
		zAlpha = (2.0 + (1.0 / zLambda) - sqrt(pow(2.0 + 1.0 / zLambda, 2.0) - 4.0)) * 0.5;
		xC0 = 1.0;
		yC0 = 1.0;
		zC0 = 1.0;
		for (i = 0; (i < 3); i++) {
		  xC0 *= (1.0 - xAlpha) * (1.0 - xAlpha);
		  yC0 *= (1.0 - yAlpha) * (1.0 - yAlpha);
		  zC0 *= (1.0 - zAlpha) * (1.0 - zAlpha);
		  xKernel[i] = xAlpha;
		  yKernel[i] = yAlpha;
		  zKernel[i] = zAlpha;
		}
		data = (double *)malloc((size_t)nx * sizeof(double));
		if (data == (double *)NULL) {
		  message("ERROR - Not enough memory for (x) data");
		  return(ERROR);
		}
		for (z = 0L; (z < nz); z++)
		  for (y = 0L; (y < ny); y++) {
			getxF2D(inPtr, 0L, y, z, nx, ny, data, nx);
			iirConvolveMirror(data, data, nx, xC0, xKernel, 3L);
			putxD2F(outPtr, 0L, y, z, nx, ny, data, nx);
		}
		free(data);
		if (ny > 1L) {
		  data = (double *)malloc((size_t)ny * sizeof(double));
		  if (data == (double *)NULL) {
			message("ERROR - Not enough memory for (y) data");
			return(ERROR);
		  }
		  for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  getyF2D(outPtr, x, 0L, z, nx, ny, data, ny);
			  iirConvolveMirror(data, data, ny, yC0, yKernel, 3L);
			  putyD2F(outPtr, x, 0L, z, nx, ny, data, ny);
			}
		  free(data);
		}
		if (nz > 1L) {
		  data = (double *)malloc((size_t)nz * sizeof(double));
		  if (data == (double *)NULL) {
			message("ERROR - Not enough memory for (z) data");
			return(ERROR);
		  }
		  for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  getzF2D(outPtr, x, y, 0L, nx, ny, data, nz);
			  iirConvolveMirror(data, data, nz, zC0, zKernel, 3L);
			  putzD2F(outPtr, x, y, 0L, nx, ny, data, nz);
			}
		  free(data);
		}
		return(!ERROR);
} /* End of anisoGauss3d */

/************************************************************************/
static	int			gaussj			(double				a[][3],
									 double				b[][3]) {

		int					indxc[3], indxr[3];
		int					ipiv[3];
		double				big, dum, pivinv, swap;
		int					i, icol, irow, j, k, l, ll;
		int					n = 3, m = 3;

		for (j = 0; (j < n); j++)
		  ipiv[j] = 0;
		for (i = 0; (i < n); i++) {
		  big = 0.0;
		  for (j = 0; (j < n); j++)
			if (ipiv[j] != 1)
			  for (k = 0; (k < n); k++) {
				if (ipiv[k] == 0) {
				  if (fabs(a[j][k]) >= big) {
					big = fabs(a[j][k]);
					irow = j;
					icol = k;
				  }
				}
				else if (ipiv[k] > 1) {
				  message("ERROR - Singular matrix in Gauss inversion");
				  return(ERROR);
				}
			  }
		  ++(ipiv[icol]);
		  if (irow != icol) {
			for (l = 0; (l < n); l++) {
			  swap = a[irow][l];
			  a[irow][l] = a[icol][l];
			  a[icol][l] = swap;
			}
			for (l = 0; (l < m); l++) {
			  swap = b[irow][l];
			  b[irow][l] = b[icol][l];
			  b[icol][l] = swap;
			}
		  }
		  indxr[i] = irow;
		  indxc[i] = icol;
		  if ((a[icol][icol] * a[icol][icol]) == 0.0) {
			message("ERROR - Singular matrix in Gauss inversion");
			return(ERROR);
		  }
		  pivinv = 1.0 / a[icol][icol];
		  a[icol][icol] = 1.0;
		  for (l = 0; (l < n); l++)
			a[icol][l] *= pivinv;
		  for (l = 0; (l < m); l++)
			b[icol][l] *= pivinv;
		  for (ll = 0; (ll < n); ll++)
			if (ll != icol) {
			  dum = a[ll][icol];
			  a[ll][icol] = 0.0;
			  for (l = 0; (l < n); l++)
				a[ll][l] -= a[icol][l] * dum;
			  for (l = 0; (l < m); l++)
				b[ll][l] -= b[icol][l] * dum;
			}
		}
		for (l = n - 1; (l >= 0); l--) {
		  if (indxr[l] != indxc[l])
			for (k = 0; (k < n); k++) {
			  swap = a[k][indxr[l]];
			  a[k][indxr[l]] = a[k][indxc[l]];
			  a[k][indxc[l]] = swap;
			}
		}
		return(!ERROR);
} /* End of gaussj */

/************************************************************************/
/* FUNCTION: absDeterminant												*/
/************************************************************************/
double				absDeterminant	(double				skew[][3]) {

		return(fabs(skew[0][0] * skew[1][1] * skew[2][2] + skew[1][0] * skew[2][1] * skew[0][2]
		  + skew[2][0] * skew[0][1] * skew[1][2] - skew[0][0] * skew[1][2] * skew[2][1]
		  - skew[0][1] * skew[1][0] * skew[2][2] - skew[0][2] * skew[1][1] * skew[2][0]));
} /* End of absDeterminant */

/************************************************************************/
/* FUNCTION: adornImage													*/
/************************************************************************/
void				adornImage		(float				*outPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 int				clipNx,
									 int				clipNy,
									 int				clipNz,
									 int				clipping,
									 float				backgrnd) {

		int					x, y, z;

		for (z = 0; (z < clipNz); z++) {
		  for (y = 0; (y < clipNy); y++) {
			for (x = 0; (x < clipNx); outPtr++, mskPtr++, x++)
			  if (clipping && (*mskPtr == 0.0))
				*outPtr = backgrnd;
			for (x = clipNx; (x < nx); outPtr++, mskPtr++, x++) {
			  *mskPtr = (float)FALSE;
			  if (clipping)
				*outPtr = backgrnd;
			}
		  }
		  for (y = clipNy; (y < ny); y++)
			for (x = 0; (x < nx); outPtr++, mskPtr++, x++) {
			  *mskPtr = (float)FALSE;
			  if (clipping)
				*outPtr = backgrnd;
			}
		}
		for (z = clipNz; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); outPtr++, mskPtr++, x++) {
			  *mskPtr = (float)FALSE;
			  if (clipping)
				*outPtr = backgrnd;
			}
} /* End of adornImage */

/************************************************************************/
/* FUNCTION: combineMasks												*/
/************************************************************************/
void				combineMasks	(float				*inMsk,
									 float				*outMsk,
									 enum	mBrands		maskCombine,
									 int				nx,
									 int				ny,
									 int				nz) {

		int					x, y, z;

		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inMsk++, outMsk++, x++)
			  switch (maskCombine) {
				case or:
				  *outMsk = (float)((*inMsk != 0.0F) || (*outMsk != 0.0F));
				  break;
				case nor:
				  *outMsk = (float)((*inMsk == 0.0F) && (*outMsk == 0.0F));
				  break;
				case and:
				  *outMsk = (float)((*inMsk != 0.0F) && (*outMsk != 0.0F));
				  break;
				case nand:
				  *outMsk = (float)((*inMsk == 0.0F) || (*outMsk == 0.0F));
				  break;
				case xor:
				  *outMsk = (float)(*inMsk != *outMsk);
				  break;
				case nxor:
				  *outMsk = (float)(*inMsk == *outMsk);
				  break;
			  }
} /* End of combineMasks */

/************************************************************************/
/* FUNCTION: computeMskSnr												*/
/************************************************************************/
void				computeMskSnr	(double				*mse,
									 double				*snr,
									 float				*inPtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 enum	mBrands		maskCombine,
									 int				nx,
									 int				ny,
									 int				nz) {

		double				sum;
		double				sumIn, sumNoise;
		int					combination;
		int					x, y, z;

		sum = sumIn = sumNoise = 0.0;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); inPtr++, outPtr++, inMsk++, outMsk++, x++) {
			  switch (maskCombine) {
				case or:
				  combination = (*inMsk != 0.0F) || (*outMsk != 0.0F);
				  break;
				case nor:
				  combination = (*inMsk == 0.0F) && (*outMsk == 0.0F);
				  break;
				case and:
				  combination = (*inMsk != 0.0F) && (*outMsk != 0.0F);
				  break;
				case nand:
				  combination = (*inMsk == 0.0F) || (*outMsk == 0.0F);
				  break;
				case xor:
				  combination = (*inMsk != *outMsk);
				  break;
				case nxor:
				  combination = (*inMsk == *outMsk);
				  break;
			  }
			  if (combination) {
				sumIn += (double)*inPtr * (double)*inPtr;
				sumNoise += ((double)*inPtr - (double)*outPtr) * ((double)*inPtr
				  - (double)*outPtr);
				sum += 1.0;
			  }
			}
		if ((sumNoise == 0.0) || (sumIn == 0.0))
		  *snr = 0.0;
		else
		  *snr = 10.0 * log10(sumIn / sumNoise);
		if (sum == 0.0)
		  sumNoise = 0.0;
		else
		  sumNoise /= sum;
		*mse = sumNoise;
} /* End of computeMskSnr */

/************************************************************************/
/* FUNCTION: computeRotChisq											*/
/************************************************************************/
void				computeRotChisq	(double				beta[],
									 float				*inPtr,
									 float				*mskPtr,
									 float				*outPtr,
									 float				*gradxPtr,
									 float				*gradyPtr,
									 float				*gradzPtr,
									 double				*chisq,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				xO,
									 double				yO,
									 double				zO) {

		double				gradx, grady, gradz;
		double				diff;
		double				normalization = 0.0;
		double				x, y, z;
		int					i;

		*chisq = 0.0;
		for (i = 0; (i < hessianSize); i++)
		  beta[i] = 0.0;
		for (z = -zO; (z < ((double)nz - zO)); z += 1.0)
		  for (y = -yO; (y < ((double)ny - yO)); y += 1.0)
			for (x = -xO; (x < ((double)nx - xO)); gradxPtr++, gradyPtr++, gradzPtr++,
			  inPtr++, outPtr++, mskPtr++, x += 1.0)
			  if (*mskPtr != 0.0) {
				diff = (double)*outPtr - (double)*inPtr;
				gradx = (double)*gradxPtr;
				grady = (double)*gradyPtr;
				gradz = (double)*gradzPtr;
				beta[0] -= diff * gradx;
				beta[1] -= diff * grady;
				beta[2] -= diff * gradz;
				beta[3] += diff * (-z * grady + y * gradz);
				beta[4] += diff * (z * gradx - x * gradz);
				beta[5] += diff * (-y * gradx + x * grady);
				beta[6] += diff * (x * gradx + y * grady + z * gradz);
				beta[7] += diff * (double)*inPtr;
				*chisq += diff * diff;
				normalization += 1.0;
			  }
		if (normalization != 0.0) {
		  for (i = 0; (i < hessianSize); i++)
			beta[i] /= normalization;
		  *chisq /= normalization;
		}
} /* End of computeRotChisq */

/************************************************************************/
/* FUNCTION: computeSkewChisq											*/
/************************************************************************/
void				computeSkewChisq(double				beta[],
									 float				*inPtr,
									 float				*mskPtr,
									 float				*outPtr,
									 float				*gradxPtr,
									 float				*gradyPtr,
									 float				*gradzPtr,
									 double				*chisq,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				xO,
									 double				yO,
									 double				zO) {

		double				gradx, grady, gradz;
		double				diff;
		double				normalization = 0.0;
		double				x, y, z;
		int					i;

		*chisq = 0.0;
		for (i = 0; (i < hessianSize); i++)
		  beta[i] = 0.0;
		for (z = -zO; (z < ((double)nz - zO)); z += 1.0)
		  for (y = -yO; (y < ((double)ny - yO)); y += 1.0)
			for (x = -xO; (x < ((double)nx - xO)); gradxPtr++, gradyPtr++, gradzPtr++,
			  inPtr++, outPtr++, mskPtr++, x += 1.0)
			  if (*mskPtr != 0.0) {
				diff = (double)*outPtr - (double)*inPtr;
				gradx = (double)*gradxPtr * diff;
				grady = (double)*gradyPtr * diff;
				gradz = (double)*gradzPtr * diff;
				beta[0] -= gradx;
				beta[1] -= grady;
				beta[2] -= gradz;
				beta[3] += x * gradx;
				beta[4] += y * gradx;
				beta[5] += z * gradx;
				beta[6] += x * grady;
				beta[7] += y * grady;
				beta[8] += z * grady;
				beta[9] += x * gradz;
				beta[10] += y * gradz;
				beta[11] += z * gradz;
				beta[12] += diff * (double)*inPtr;
				*chisq += diff * diff;
				normalization += 1.0;
			  }
		if (normalization != 0.0) {
		  for (i = 0; (i < hessianSize); i++)
			beta[i] /= normalization;
		  *chisq /= normalization;
		}
} /* End of computeSkewChisq */

/************************************************************************/
/* FUNCTION: copyClip													*/
/************************************************************************/
void				copyClip		(float				*inPtr,
									 int				nxIn,
									 int				nyIn,
									 int				nzIn,
									 float				*outPtr,
									 int				nxOut,
									 int				nyOut,
									 int				nzOut) {

		int					x, y, z;
		int					minNx, minNy, minNz;

		minNx = (nxIn < nxOut) ? (nxIn) : (nxOut);
		minNy = (nyIn < nyOut) ? (nyIn) : (nyOut);
		minNz = (nzIn < nzOut) ? (nzIn) : (nzOut);
		for (z = 0; (z < minNz); z++) {
		  for (y = 0; (y < minNy); y++) {
			for (x = 0; (x < minNx); outPtr++, inPtr++, x++)
			  *outPtr = *inPtr;
			while (x < nxOut) {
			  *outPtr = 0.0F;
			  outPtr++;
			  x++;
			}
			while (x < nxIn) {
			  inPtr++;
			  x++;
			}
		  }
		  while (y < nyOut) {
			for (x = 0; (x < nxOut); outPtr++, x++)
			  *outPtr = 0.0F;
			y++;
		  }
		  while (y < nyIn) {
			inPtr += nxIn;
			y++;
		  }
		}
		while (z < nzOut) {
		  for (y = 0; (y < nyOut); y++)
			for (x = 0; (x < nxOut); outPtr++, x++)
			  *outPtr = 0.0F;
		  z++;
		}
} /* End of copyClip */

/************************************************************************/
/* FUNCTION: dirMskTransform											*/
/************************************************************************/
int					dirMskTransform	(struct fitRec		*fit,
									 float				*inPtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 int				greyRendering,
									 enum	iDegree		interpolation) {

		double				*wtx, *wty, *wtz;
		double				*splinePtr, *p;
		double				*dpi, *dpj;
		float				*fpi, *fpj;
		long				*fdx, *fdy, *fdz;
		double				weightX[9], weightY[9], weightZ[9];
		long				foldX[9], foldY[9], foldZ[9];
		ptrdiff_t			d;
		double				xO, yO, zO;
		double				xz, yz, zz;
		double				xy, yy, zy;
		double				xIn, yIn, zIn;
		double				xOut, yOut, zOut;
		double				q, qi, qj;
		long				lnx = (long)nx, lny = (long)ny, lnz = (long)nz;
		long				nxy = lnx * lny;
		long				x, y, z;
		long				i, j, k;
		long				half, width;
		int					degree;

		xO = fit->origin[0];
		yO = fit->origin[1];
		zO = fit->origin[2];
		switch (interpolation) {
		  case zero:
			degree = 0;
			break;
		  case one:
			degree = 1;
			break;
		  case three:
			degree = 3;
			break;
		  default:
			message("ERROR - Unknown interpolation specification");
			return(ERROR);
		}
		half = (long)degree / 2L + 1L;
		width = 2L * half + 1L;
		switch (interpolation) {
		  case zero:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = fit->skew[0][2] * zOut - fit->dx[0] + xO;
			  yz = fit->skew[1][2] * zOut - fit->dx[1] + yO;
			  zz = fit->skew[2][2] * zOut - fit->dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + fit->skew[0][1] * yOut;
				yy = yz + fit->skew[1][1] * yOut;
				zy = zz + fit->skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + fit->skew[0][0] * xOut + 0.5;
				  yIn = yy + fit->skew[1][0] * xOut + 0.5;
				  zIn = zy + fit->skew[2][0] * xOut + 0.5;
				  x = (long)xIn;
				  if (xIn < 0.0)
					x--;
				  y = (long)yIn;
				  if (yIn < 0.0)
					y--;
				  z = (long)zIn;
				  if (zIn < 0.0)
					z--;
				  if ((0L <= x) && (x < lnx) && (0L <= y) && (y < lny) && (0L <= z)
					&& (z < lnz)) {
					d = (ptrdiff_t)(z * nxy + y * lnx + x);
					if (greyRendering)
					  *outPtr++ = (float)((double)*(inPtr + d) * exp(fit->gamma));
					else
					  *outPtr++ = *(inPtr + d);
					*outMsk++ = *(inMsk + d);
				  }
				  else {
					*outPtr++ = 0.0F;
					*outMsk++ = 0.0F;
				  }
				}
			  }
			}
			break;
		  case one:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = fit->skew[0][2] * zOut - fit->dx[0] + xO;
			  yz = fit->skew[1][2] * zOut - fit->dx[1] + yO;
			  zz = fit->skew[2][2] * zOut - fit->dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + fit->skew[0][1] * yOut;
				yy = yz + fit->skew[1][1] * yOut;
				zy = zz + fit->skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + fit->skew[0][0] * xOut + 0.5;
				  yIn = yy + fit->skew[1][0] * xOut + 0.5;
				  zIn = zy + fit->skew[2][0] * xOut + 0.5;
				  x = (long)xIn;
				  if (xIn < 0.0)
					x--;
				  y = (long)yIn;
				  if (yIn < 0.0)
					y--;
				  z = (long)zIn;
				  if (zIn < 0.0)
					z--;
				  if ((0L <= x) && (x < lnx) && (0L <= y) && (y < lny) && (0L <= z)
					&& (z < lnz)) {
					d = (ptrdiff_t)(z * nxy + y * lnx + x);
					*outMsk = *(inMsk + d);
				  }
				  else {
					d = (ptrdiff_t)(-1);
					*outMsk = 0.0F;
				  }
				  if (*outMsk++ != 0.0F) {
					xIn -= 0.5;
					yIn -= 0.5;
					zIn -= 0.5;
					wtz = weightZ;
					fdz = foldZ;
					for (k = z - 1L; (k <= (z + 1L)); k++) {
					  *wtz++ = (lnz == 1L) ? (1.0) : (BsplnWght1(k, zIn));
					  *fdz++ = (lnz == 1L) ? (0L) : (fold(k, lnz) * nxy);
					}
					wty = weightY;
					fdy = foldY;
					for (j = y - 1L; (j <= (y + 1L)); j++) {
					  *wty++ = BsplnWght1(j, yIn);
					  *fdy++ = fold(j, lny) * lnx;
					}
					wtx = weightX;
					fdx = foldX;
					for (i = x - 1L; (i <= (x + 1L)); i++) {
					  *wtx++ = BsplnWght1(i, xIn);
					  *fdx++ = fold(i, lnx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (lnz == 1L) ? (1L) : (3L); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  fpj = inPtr + (ptrdiff_t)(*fdz++);
					  qj = 0.0;
					  for (j = 3L; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						fpi = fpj + (ptrdiff_t)(*fdy++);
						qi = 0.0;
						for (i = 3L; (i-- > 0L);)
						  qi += *wtx++ * (double)*(fpi + (ptrdiff_t)(*fdx++));
						qj += qi * *wty++;
					  }
					  q += qj * *wtz++;
					}
				  }
				  else
					q = (d == (ptrdiff_t)(-1)) ? (0.0) : ((double)*(inPtr + d));
				  if (greyRendering)
					q *= exp(fit->gamma);
				  *outPtr++ = (float)q;
				}
			  }
			}
			break;
		  case three:
			splinePtr = (double *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(double));
			if (splinePtr == (double *)NULL) {
			  message("ERROR - Not enough memory for B-spline coefficients");
			  return(ERROR);
			}
			p = splinePtr;
			for (z = 0L; (z < lnz); z++)
			  for (y = 0L; (y < lny); p += (ptrdiff_t)nx, y++)
				getxF2D(inPtr, 0L, y, z, lnx, lny, p, lnx);
			if (directBsplineMirror(splinePtr, splinePtr, nx, ny, nz, degree) == ERROR) {
			  free(splinePtr);
			  message("ERROR - Unable to compute B-spline coefficients");
			  return(ERROR);
			}
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = fit->skew[0][2] * zOut - fit->dx[0] + xO;
			  yz = fit->skew[1][2] * zOut - fit->dx[1] + yO;
			  zz = fit->skew[2][2] * zOut - fit->dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + fit->skew[0][1] * yOut;
				yy = yz + fit->skew[1][1] * yOut;
				zy = zz + fit->skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + fit->skew[0][0] * xOut + 0.5;
				  yIn = yy + fit->skew[1][0] * xOut + 0.5;
				  zIn = zy + fit->skew[2][0] * xOut + 0.5;
				  x = (long)xIn;
				  if (xIn < 0.0)
					x--;
				  y = (long)yIn;
				  if (yIn < 0.0)
					y--;
				  z = (long)zIn;
				  if (zIn < 0.0)
					z--;
				  if ((0L <= x) && (x < lnx) && (0L <= y) && (y < lny) && (0L <= z)
					&& (z < lnz)) {
					d = (ptrdiff_t)(z * nxy + y * lnx + x);
					*outMsk = *(inMsk + d);
				  }
				  else {
					d = (ptrdiff_t)(-1);
					*outMsk = 0.0F;
				  }
				  if (*outMsk++ != 0.0F) {
					xIn -= 0.5;
					yIn -= 0.5;
					zIn -= 0.5;
					wtz = weightZ;
					fdz = foldZ;
					for (k = z - half; (k <= (z + half)); k++) {
					  *wtz++ = (nz == 1L) ? (1.0) : (BsplnWght(degree, k, zIn));
					  *fdz++ = (nz == 1L) ? (0L) : (fold(k, nz) * nxy);
					}
					wty = weightY;
					fdy = foldY;
					for (j = y - half; (j <= (y + half)); j++) {
					  *wty++ = BsplnWght(degree, j, yIn);
					  *fdy++ = fold(j, ny) * nx;
					}
					wtx = weightX;
					fdx = foldX;
					for (i = x - half; (i <= (x + half)); i++) {
					  *wtx++ = BsplnWght(degree, i, xIn);
					  *fdx++ = fold(i, nx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (nz == 1L) ? (1L) : (width); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  dpj = splinePtr + (ptrdiff_t)(*fdz++);
					  qj = 0.0;
					  for (j = width; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						dpi = dpj + (ptrdiff_t)(*fdy++);
						qi = 0.0;
						for (i = width; (i-- > 0L);)
						  qi += *wtx++ * (double)*(dpi + (ptrdiff_t)(*fdx++));
						qj += qi * *wty++;
					  }
					  q += qj * *wtz++;
					}
				  }
				  else
					q = (d == (ptrdiff_t)(-1)) ? (0.0) : ((double)*(inPtr + d));
				  if (greyRendering)
					q *= exp(fit->gamma);
				  *outPtr++ = (float)q;
				}
			  }
			}
			free(splinePtr);
			break;
		  default:
			message("ERROR - Unknown interpolation specification");
			return(ERROR);
		}
		return(!ERROR);
} /* End of dirMskTransform */

/************************************************************************/
/* FUNCTION: eulerRotation												*/
/************************************************************************/
void				eulerRotation	(struct	fitRec		*fit) {

		double				Rx[3][3], Ry[3][3], Rz[3][3];
		double				Rxy[3][3];
		int					i, j, k;

		for (i = 0; (i < 3); i++) {
		  for (j = 0; (j < 3); j++) {
			Rx[i][j] = 0.0;
			Ry[i][j] = 0.0;
			Rz[i][j] = 0.0;
		  }
		  Rx[i][i] = 1.0;
		  Ry[i][i] = 1.0;
		  Rz[i][i] = 1.0;
		}
		Rx[1][1] = cos(fit->phi);
		Rx[1][2] = -sin(fit->phi);
		Rx[2][1] = sin(fit->phi);
		Rx[2][2] = cos(fit->phi);
		Ry[0][0] = cos(fit->theta);
		Ry[0][2] = sin(fit->theta);
		Ry[2][0] = -sin(fit->theta);
		Ry[2][2] = cos(fit->theta);
		Rz[0][0] = cos(fit->psi);
		Rz[0][1] = -sin(fit->psi);
		Rz[1][0] = sin(fit->psi);
		Rz[1][1] = cos(fit->psi);
		for (i = 0; (i < 3); i++)
		  for (j = 0; (j < 3); j++)
			for (k = 0, Rxy[i][j] = 0.0; (k < 3); k++)
			  Rxy[i][j] += Rx[i][k] * Ry[k][j];
		for (i = 0; (i < 3); i++)
		  for (j = 0; (j < 3); j++)
			for (k = 0, fit->skew[i][j] = 0.0; (k < 3); k++)
			  fit->skew[i][j] += Rxy[i][k] * Rz[k][j];
} /* End of eulerRotation */

/************************************************************************/
/* FUNCTION: exportFit													*/
/************************************************************************/
int					exportFit		(struct fitRec		*fit) {

		char				string[MSGLENGTH];

		if (sprintf(string, " $dx1=%+17.9E, $dx2=%+17.9E, $dx3=%+17.9E",
		  fit->dx[0], fit->dx[1], fit->dx[2]) == EOF) {
		  message("ERROR - Unable to display $dx1 nor $dx2 nor $dx3");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $a11=%+17.9E, $a12=%+17.9E, $a13=%+17.9E",
		  fit->skew[0][0], fit->skew[0][1], fit->skew[0][2]) == EOF) {
		  message("ERROR - Unable to display $a11 nor $a12 nor $a13");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $a21=%+17.9E, $a22=%+17.9E, $a23=%+17.9E",
		  fit->skew[1][0], fit->skew[1][1], fit->skew[1][2]) == EOF) {
		  message("ERROR - Unable to display $a21 nor $a22 nor $a23");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $a31=%+17.9E, $a32=%+17.9E, $a33=%+17.9E",
		  fit->skew[2][0], fit->skew[2][1], fit->skew[2][2]) == EOF) {
		  message("ERROR - Unable to display $a31 nor $a32 nor $a33");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $phi=%+17.9E, $tht=%+17.9E, $psi=%+17.9E",
		  fit->phi * 180.0 / Pi, fit->theta * 180.0 / Pi, fit->psi * 180.0 / Pi) == EOF) {
		  message("ERROR - Unable to display $phi nor $tht nor $psi");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $lbd=%+17.9E, $gam=%+17.9E", fit->lambda, fit->gamma) == EOF) {
		  message("ERROR - Unable to display $lbd nor $gam");
		  return(ERROR);
		}
		message(string);
		return(!ERROR);
} /* End of exportFit */

/************************************************************************/
/* FUNCTION: exportSummary												*/
/************************************************************************/
int					exportSummary	(struct fitRec		*fit,
									 double				mse,
									 double				snr) {

		char				string[MSGLENGTH];
		double				det, skw, rot;
		double				r, s[3], n[3];
		int					i, j;

		det = exp(log(absDeterminant(fit->skew)) / 3.0);
		j = 0;
		for (s[0] = 0.0, n[j] = 0.0, i = 0; (i < 3); i++) {
		  s[0] += fit->skew[i][0] * fit->skew[i][2];
		  n[j] += fit->skew[i][j] * fit->skew[i][j];
		}
		skw = s[0] * s[0];
		for (j = 1; (j < 3); j++) {
		  for (s[j] = 0.0, n[j] = 0.0, i = 0; (i < 3); i++) {
			s[j] += fit->skew[i][j] * fit->skew[i][j - 1];
			n[j] += fit->skew[i][j] * fit->skew[i][j];
		  }
		  skw += s[j] * s[j];
		}
		j = 0;
		rot = fabs(acos(fit->skew[j][j] / sqrt(n[j])));
		for (j = 1; (j < 3); j++) {
		  r = fabs(acos(fit->skew[j][j] / sqrt(n[j])));
		  if (r > rot)
			rot = r;
		}
		rot = rot * 180.0 / Pi;
		if (sprintf(string, " $det=%+17.9E, $skw=%+17.9E, $rot=%+17.9E", det, skw, rot) == EOF) {
		  message("ERROR - Unable to display $det nor $skw nor $rot");
		  return(ERROR);
		}
		message(string);
		if (sprintf(string, " $err=%+17.9E, $snr=%+17.9E dB", mse, snr) == EOF) {
		  message("ERROR - Unable to display $eer nor $snr");
		  return(ERROR);
		}
		message(string);
		return(!ERROR);
} /* End of exportSummary */

/************************************************************************/
/* FUNCTION: fExportFit													*/
/************************************************************************/
int					fExportFit		(struct fitRec		*fit,
									 char				*fileName) {

		FILE				*fitFile;

		fitFile = fopen(fileName, "w");
		if (fitFile == (FILE *)NULL) {
		  message("ERROR - Unable to open fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$Ox1 = %+25.17E\n", fit->origin[0]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$Ox2 = %+25.17E\n", fit->origin[1]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$Ox3 = %+25.17E\n", fit->origin[2]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$dx1 = %+25.17E\n", fit->dx[0]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$dx2 = %+25.17E\n", fit->dx[1]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$dx3 = %+25.17E\n", fit->dx[2]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a11 = %+25.17E\n", fit->skew[0][0]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a12 = %+25.17E\n", fit->skew[0][1]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a13 = %+25.17E\n", fit->skew[0][2]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a21 = %+25.17E\n", fit->skew[1][0]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a22 = %+25.17E\n", fit->skew[1][1]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a23 = %+25.17E\n", fit->skew[1][2]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a31 = %+25.17E\n", fit->skew[2][0]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a32 = %+25.17E\n", fit->skew[2][1]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$a33 = %+25.17E\n", fit->skew[2][2]) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$gam = %+25.17E\n", fit->gamma) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$phi = %+25.17E\n", fit->phi * 180.0 / Pi) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$tht = %+25.17E\n", fit->theta * 180.0 / Pi) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$psi = %+25.17E\n", fit->psi * 180.0 / Pi) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fprintf(fitFile, "$lbd = %+25.17E\n", fit->lambda) == EOF) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Problem with fit file");
		  return(ERROR);
		}
		if (fclose(fitFile) == EOF) {
		  message("ERROR - Unable to close fit file");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of fExportFit */

/************************************************************************/
/* FUNCTION: importFit													*/
/************************************************************************/
int					importFit		(struct fitRec		*fit,
									 char				*fitName) {

		FILE				*fitFile;
		double				dummyDouble;
		int					index;

		fitFile = fopen(fitName, "r");
		if (fitFile == (FILE *)NULL) {
		  message("ERROR - Unable to open fit file");
		  return(ERROR);
		}
		if (fscanf(fitFile, "$Ox%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $Ox1 variable");
		  return(ERROR);
		}
		if (index != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $Ox1 variable");
		  return(ERROR);
		}
		fit->origin[0] = dummyDouble;
		if (fscanf(fitFile, "$Ox%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $Ox2 variable");
		  return(ERROR);
		}
		if (index != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $Ox2 variable");
		  return(ERROR);
		}
		fit->origin[1] = dummyDouble;
		if (fscanf(fitFile, "$Ox%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $Ox3 variable");
		  return(ERROR);
		}
		if (index != 3) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $Ox3 variable");
		  return(ERROR);
		}
		fit->origin[2] = dummyDouble;
		if (fscanf(fitFile, "$dx%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $dx1 variable");
		  return(ERROR);
		}
		if (index != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $dx1 variable");
		  return(ERROR);
		}
		fit->dx[0] = dummyDouble;
		if (fscanf(fitFile, "$dx%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $dx2 variable");
		  return(ERROR);
		}
		if (index != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $dx2 variable");
		  return(ERROR);
		}
		fit->dx[1] = dummyDouble;
		if (fscanf(fitFile, "$dx%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $dx3 variable");
		  return(ERROR);
		}
		if (index != 3) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $dx3 variable");
		  return(ERROR);
		}
		fit->dx[2] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a11 variable");
		  return(ERROR);
		}
		if (index != 11) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a11 variable");
		  return(ERROR);
		}
		fit->skew[0][0] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a12 variable");
		  return(ERROR);
		}
		if (index != 12) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a12 variable");
		  return(ERROR);
		}
		fit->skew[0][1] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a13 variable");
		  return(ERROR);
		}
		if (index != 13) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a13 variable");
		  return(ERROR);
		}
		fit->skew[0][2] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a21 variable");
		  return(ERROR);
		}
		if (index != 21) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a21 variable");
		  return(ERROR);
		}
		fit->skew[1][0] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a22 variable");
		  return(ERROR);
		}
		if (index != 22) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a22 variable");
		  return(ERROR);
		}
		fit->skew[1][1] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a23 variable");
		  return(ERROR);
		}
		if (index != 23) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a23 variable");
		  return(ERROR);
		}
		fit->skew[1][2] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a31 variable");
		  return(ERROR);
		}
		if (index != 31) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a31 variable");
		  return(ERROR);
		}
		fit->skew[2][0] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a32 variable");
		  return(ERROR);
		}
		if (index != 32) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a32 variable");
		  return(ERROR);
		}
		fit->skew[2][1] = dummyDouble;
		if (fscanf(fitFile, "$a%d = %lE\n", &index, &dummyDouble) != 2) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to access $a33 variable");
		  return(ERROR);
		}
		if (index != 33) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $a33 variable");
		  return(ERROR);
		}
		fit->skew[2][2] = dummyDouble;
		if (fscanf(fitFile, "$gam = %lE\n", &dummyDouble) != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $gam variable");
		  return(ERROR);
		}
		fit->gamma = dummyDouble;
		if (fscanf(fitFile, "$phi = %lE\n", &dummyDouble) != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $phi variable");
		  return(ERROR);
		}
		fit->phi = dummyDouble * Pi / 180.0F;
		if (fscanf(fitFile, "$tht = %lE\n", &dummyDouble) != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $tht variable");
		  return(ERROR);
		}
		fit->theta = dummyDouble * Pi / 180.0F;
		if (fscanf(fitFile, "$psi = %lE\n", &dummyDouble) != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $psi variable");
		  return(ERROR);
		}
		fit->psi = dummyDouble * Pi / 180.0F;
		if (fscanf(fitFile, "$lbd = %lE\n", &dummyDouble) != 1) {
		  if (fclose(fitFile) == EOF)
			message("ERROR - Unable to close fit file");
		  message("ERROR - Unable to recover $lbd variable");
		  return(ERROR);
		}
		fit->lambda = dummyDouble;
		if (fclose(fitFile) == EOF) {
		  message("ERROR - Unable to close fit file");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of importFit */

/************************************************************************/
/* FUNCTION: initialEstimate											*/
/************************************************************************/
void				initialEstimate	(struct fitRec		*fit,
									 int				nx,
									 int				ny,
									 int				nz) {

		fit->origin[0] = (double)((nx - 1) / 2);
		fit->origin[1] = (double)((ny - 1) / 2);
		fit->origin[2] = (double)((nz - 1) / 2);
		fit->dx[0] = 0.0;
		fit->dx[1] = 0.0;
		fit->dx[2] = 0.0;
		fit->skew[0][0] = 1.0;
		fit->skew[0][1] = 0.0;
		fit->skew[0][2] = 0.0;
		fit->skew[1][0] = 0.0;
		fit->skew[1][1] = 1.0;
		fit->skew[1][2] = 0.0;
		fit->skew[2][0] = 0.0;
		fit->skew[2][1] = 0.0;
		fit->skew[2][2] = 1.0;
		fit->gamma = 0.0;
		fit->phi = 0.0;
		fit->theta = 0.0;
		fit->psi = 0.0;
		fit->lambda = 0.0;
} /* End of initialEstimate */

/************************************************************************/
/* FUNCTION: invertFit													*/
/************************************************************************/
int					invertFit		(struct fitRec		*dirFit,
									 struct fitRec		*invFit) {

		double				dirSkew[3][3];
		int					i, j;

		for (j = 0; (j < 3); j++) {
		  invFit->origin[j] = dirFit->origin[j];
		  for (i = 0; (i < 3); i++) {
			dirSkew[i][j] = dirFit->skew[i][j];
			invFit->skew[i][j] = 0.0;
		  }
		  invFit->skew[j][j] = 1.0;
		}
		if (gaussj(dirSkew, invFit->skew) == ERROR) {
		  message("ERROR - Non reversible transformation");
		  return(ERROR);
		}
		for (j = 0; (j < 3); j++) {
		  invFit->dx[j] = 0.0;
		  for (i = 0; (i < 3); i++)
			invFit->dx[j] -= invFit->skew[j][i] * dirFit->dx[i];
		}
		invFit->gamma = -dirFit->gamma;
		invFit->phi = 0.0;
		invFit->theta = 0.0;
		invFit->psi = 0.0;
		invFit->lambda = -dirFit->lambda;
		return(!ERROR);
} /* End of invertFit */

/************************************************************************/
/* FUNCTION: invTransform												*/
/************************************************************************/
int					invTransform	(struct fitRec		*fit,
									 float				*inPtr,
									 double				*splinePtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 enum	iDegree		interpolation) {

		double				*wtx, *wty, *wtz;
		double				*dpi, *dpj;
		float				*fpi, *fpj;
		long				*fdx, *fdy, *fdz;
		double				weightX[5], weightY[5], weightZ[5];
		long				foldX[5], foldY[5], foldZ[5];
		struct fitRec		invFit;
		ptrdiff_t			d;
		double				xO, yO, zO;
		double				xz, yz, zz;
		double				xy, yy, zy;
		double				xIn, yIn, zIn;
		double				xOut, yOut, zOut;
		double				q, qi, qj;
		long				lnx = (long)nx, lny = (long)ny, lnz = (long)nz;
		long				nxy = lnx * lny;
		long				x, y, z;
		long				i, j, k;

		if (invertFit(fit, &invFit) == ERROR) {
		  message("ERROR - Unable to compute the backward transformation");
		  return(ERROR);
		}
		xO = invFit.origin[0];
		yO = invFit.origin[1];
		zO = invFit.origin[2];
		switch (interpolation) {
		  case one:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invFit.skew[0][2] * zOut - invFit.dx[0] + xO;
			  yz = invFit.skew[1][2] * zOut - invFit.dx[1] + yO;
			  zz = invFit.skew[2][2] * zOut - invFit.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invFit.skew[0][1] * yOut;
				yy = yz + invFit.skew[1][1] * yOut;
				zy = zz + invFit.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + invFit.skew[0][0] * xOut + 0.5;
				  yIn = yy + invFit.skew[1][0] * xOut + 0.5;
				  zIn = zy + invFit.skew[2][0] * xOut + 0.5;
				  x = (long)xIn;
				  if (xIn < 0.0)
					x--;
				  y = (long)yIn;
				  if (yIn < 0.0)
					y--;
				  z = (long)zIn;
				  if (zIn < 0.0)
					z--;
				  if ((0L <= x) && (x < lnx) && (0L <= y) && (y < lny) && (0L <= z)
					&& (z < lnz)) {
					d = (ptrdiff_t)(z * nxy + y * lnx + x);
					*outMsk = *(inMsk + d);
				  }
				  else {
					d = (ptrdiff_t)(-1);
					*outMsk = 0.0F;
				  }
				  if (*outMsk++ != 0.0F) {
					xIn -= 0.5;
					yIn -= 0.5;
					zIn -= 0.5;
					wtz = weightZ;
					fdz = foldZ;
					for (k = z - 1L; (k <= (z + 1L)); k++) {
					  *wtz++ = (lnz == 1L) ? (1.0) : (BsplnWght1(k, zIn));
					  *fdz++ = (lnz == 1L) ? (0L) : (fold(k, lnz) * nxy);
					}
					wty = weightY;
					fdy = foldY;
					for (j = y - 1L; (j <= (y + 1L)); j++) {
					  *wty++ = BsplnWght1(j, yIn);
					  *fdy++ = fold(j, lny) * lnx;
					}
					wtx = weightX;
					fdx = foldX;
					for (i = x - 1L; (i <= (x + 1L)); i++) {
					  *wtx++ = BsplnWght1(i, xIn);
					  *fdx++ = fold(i, lnx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (lnz == 1L) ? (1L) : (3L); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  fpj = inPtr + (ptrdiff_t)(*fdz++);
					  qj = 0.0;
					  for (j = 3L; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						fpi = fpj + (ptrdiff_t)(*fdy++);
						qi = 0.0;
						for (i = 3L; (i-- > 0L);)
						  qi += *wtx++ * (double)*(fpi + (ptrdiff_t)(*fdx++));
						qj += qi * *wty++;
					  }
					  q += qj * *wtz++;
					}
				  }
				  else
					q = (d == (ptrdiff_t)(-1)) ? (0.0) : ((double)*(inPtr + d)
					  * exp(invFit.gamma));
				  *outPtr++ = (float)q;
				}
			  }
			}
			break;
		  case three:
			for (zOut = -zO; (zOut < ((double)nz - zO)); zOut += 1.0) {
			  xz = invFit.skew[0][2] * zOut - invFit.dx[0] + xO;
			  yz = invFit.skew[1][2] * zOut - invFit.dx[1] + yO;
			  zz = invFit.skew[2][2] * zOut - invFit.dx[2] + zO;
			  for (yOut = -yO; (yOut < ((double)ny - yO)); yOut += 1.0) {
				xy = xz + invFit.skew[0][1] * yOut;
				yy = yz + invFit.skew[1][1] * yOut;
				zy = zz + invFit.skew[2][1] * yOut;
				for (xOut = -xO; (xOut < ((double)nx - xO)); xOut += 1.0) {
				  xIn = xy + invFit.skew[0][0] * xOut + 0.5;
				  yIn = yy + invFit.skew[1][0] * xOut + 0.5;
				  zIn = zy + invFit.skew[2][0] * xOut + 0.5;
				  x = (long)xIn;
				  if (xIn < 0.0)
					x--;
				  y = (long)yIn;
				  if (yIn < 0.0)
					y--;
				  z = (long)zIn;
				  if (zIn < 0.0)
					z--;
				  if ((0L <= x) && (x < lnx) && (0L <= y) && (y < lny) && (0L <= z)
					&& (z < lnz)) {
					d = (ptrdiff_t)(z * nxy + y * lnx + x);
					*outMsk = *(inMsk + d);
				  }
				  else {
					d = (ptrdiff_t)(-1);
					*outMsk = 0.0F;
				  }
				  if (*outMsk++ != 0.0F) {
					xIn -= 0.5;
					yIn -= 0.5;
					zIn -= 0.5;
					wtz = weightZ;
					fdz = foldZ;
					for (k = z - 2L; (k <= (z + 2L)); k++) {
					  *wtz++ = (lnz == 1L) ? (1.0) : (BsplnWght3(k, zIn));
					  *fdz++ = (lnz == 1L) ? (0L) : (fold(k, lnz) * nxy);
					}
					wty = weightY;
					fdy = foldY;
					for (j = y - 2L; (j <= (y + 2L)); j++) {
					  *wty++ = BsplnWght3(j, yIn);
					  *fdy++ = fold(j, lny) * lnx;
					}
					wtx = weightX;
					fdx = foldX;
					for (i = x - 2L; (i <= (x + 2L)); i++) {
					  *wtx++ = BsplnWght3(i, xIn);
					  *fdx++ = fold(i, lnx);
					}
					q = 0.0;
					wtz = weightZ;
					fdz = foldZ;
					for (k = (lnz == 1L) ? (1L) : (5L); (k-- > 0L);) {
					  wty = weightY;
					  fdy = foldY;
					  dpj = splinePtr + (ptrdiff_t)(*fdz++);
					  qj = 0.0;
					  for (j = 5L; (j-- > 0L);) {
						wtx = weightX;
						fdx = foldX;
						dpi = dpj + (ptrdiff_t)(*fdy++);
						qi = 0.0;
						for (i = 5L; (i-- > 0L);)
						  qi += *wtx++ * *(dpi + (ptrdiff_t)(*fdx++));
						qj += qi * *wty++;
					  }
					  q += qj * *wtz++;
					}
				  }
				  else
					q = (d == (ptrdiff_t)(-1)) ? (0.0) : ((double)*(inPtr + d)
					  * exp(invFit.gamma));
				  *outPtr++ = (float)q;
				}
			  }
			}
			break;
		  default:
			message("ERROR - Unknown interpolation specification");
			return(ERROR);
		}
		return(!ERROR);
} /* End of invTransform */

/************************************************************************/
/* FUNCTION: maskFromData												*/
/************************************************************************/
int					maskFromData	(float				*inPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				sx,
									 double				sy,
									 double				sz) {

		float				*p;
		short				*tmp, *q;
		struct qParam		quantData;
		double				r, rounded;
		int					x, y, z;

		if (anisoGauss3d(inPtr, mskPtr, (long)nx, (long)ny, (long)nz, sx, sy, sz) == ERROR) {
		  message("ERROR - Not enough memory for low-pass filtering");
		  return(ERROR);
		}
		tmp = (short *)malloc((size_t)nx * (size_t)ny * (size_t)nz * sizeof(short));
		if (tmp == (short *)NULL) {
		  message("ERROR - Not enough memory for temporary mask");
		  return(ERROR);
		}
		quantData.inPtr.s = tmp;
		quantData.outPtr.s = tmp;
		quantData.nx = nx;
		quantData.ny = ny;
		quantData.nz = nz;
		quantData.floatInput = FALSE;
		quantData.floatOutput = FALSE;
		quantData.mode = LloydMax;
		quantData.labels = from;
		quantData.bins = 2U;
		quantData.first = (short)0;
		quantData.step = (short)1;
		quantData.epsilon = (double)FLT_EPSILON;
		p = mskPtr;
		q = tmp;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); p++, q++, x++) {
			  r = (double)*p + 0.5;
			  rounded = (double)((long)r);
			  if ((r < 0.0) && (rounded != r))
				rounded -= 1.0;
			  if (rounded < (double)SHRT_MIN)
				*q = (short)SHRT_MIN;
			  else if ((double)SHRT_MAX < rounded)
				*q = (short)SHRT_MAX;
			  else
				*q = (short)rounded;
			}
		if (histQuant(&quantData) == ERROR) {
		  free(tmp);
		  message("ERROR - Unable to execute vector quantization of data");
		  return(ERROR);
		}
		p = mskPtr;
		q = tmp;
		for (z = 0; (z < nz); z++)
		  for (y = 0; (y < ny); y++)
			for (x = 0; (x < nx); p++, q++, x++)
			  *p = (float)*q;
		free(tmp);
		return(!ERROR);
} /* End of maskFromData */
