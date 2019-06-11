#include		<stddef.h>
#include		<stdlib.h>
#include		<string.h>
#include		<math.h>
#include		"phil.h"

#include		"BsplnTrf.h"

#include		"convolve.h"
#include		"getPut.h"

/************************************************************************/
static	int			directBspline	(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree,
									 enum	bBrands		boundary);
static	int			inverseBspline	(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree,
									 enum	bBrands		boundary);

/************************************************************************/
static	int			directBspline	(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree,
									 enum	bBrands		boundary) {

		double				*data;
		double				realPoles[3];
		double				gain;
		long				np;
		long				x, y, z;

		switch (degree) {
		  case 0:
		  case 1:
			outPtr = (double *)memcpy(outPtr, inPtr, (size_t)(nx * ny * nz * (long)sizeof(double)));
			return(!ERROR);
		  case 2:
			np = 1L;
			realPoles[0] = sqrt(8.0) - 3.0;
			gain = 24.0 - sqrt(512.0);
			break;
		  case 3:
			np = 1L;
			realPoles[0] = sqrt(3.0) - 2.0;
			gain = 12.0 - sqrt(108.0);
			break;
		  case 4:
			np = 2L;
			realPoles[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			realPoles[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			gain = (1.0 - realPoles[0]) * (1.0 - realPoles[1]);
			gain *= gain;
			break;
		  case 5:
			np = 2L;
			realPoles[0] = 0.5 * (sqrt(270.0 - sqrt(70980.0)) + sqrt(105.0) - 13.0);
			realPoles[1] = 0.5 * (sqrt(270.0 + sqrt(70980.0)) - sqrt(105.0) - 13.0);
			gain = (1.0 - realPoles[0]) * (1.0 - realPoles[1]);
			gain *= gain;
			break;
		  case 6:
			np = 3L;
			realPoles[0] = -0.488294589303044755130118038883789062112279161239377608394;
			realPoles[1] = -0.081679271076237512597937765737059080653379610398148178525368;
			realPoles[2] = -0.00141415180832581775108724397655859252786416905534669851652709;
			gain = 2.598975999348577818390170516255374207847876853191217652822;
			break;
		  case 7:
			np = 3L;
			realPoles[0] = -0.5352804307964381655424037816816460718339231523426924148812;
			realPoles[1] = -0.122554615192326690515272264359357343605486549427295558490763;
			realPoles[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
			gain = 3.0248282036441843886795463832305782146916878615537002580987;
			break;
		  default:
			message("ERROR - Unknown degree");
			return(ERROR);
		}
		data = outPtr;
		for (z = 0L; (z < nz); z++)
		  for (y = 0L; (y < ny); inPtr += (ptrdiff_t)nx, data += (ptrdiff_t)nx, y++)
			switch (boundary) {
			  case null:
				iirConvolveFinite(inPtr, data, nx, gain, realPoles, np);
				break;
			  case periodic:
				iirConvolvePeriodic(inPtr, data, nx, gain, realPoles, np);
				break;
			  case mirror:
				iirConvolveMirror(inPtr, data, nx, gain, realPoles, np);
				break;
			  default:
				message("ERROR - Unknown boundary");
				return(ERROR);
			}
		if (ny > 1L) {
		  data = (double *)malloc((size_t)ny * sizeof(double));
		  if (data == (double *)NULL) {
			message("ERROR - Unable to allocate y data");
			return(ERROR);
		  }
		  for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  getyD2D(outPtr, x, 0L, z, nx, ny, data, ny);
			  switch (boundary) {
				case null:
				  iirConvolveFinite(data, data, ny, gain, realPoles, np);
				  break;
				case periodic:
				  iirConvolvePeriodic(data, data, ny, gain, realPoles, np);
				  break;
				case mirror:
				  iirConvolveMirror(data, data, ny, gain, realPoles, np);
				  break;
				default:
				  free(data);
				  message("ERROR - Unknown boundary");
				  return(ERROR);
			  }
			  putyD2D(outPtr, x, 0L, z, nx, ny, data, ny);
			}
		  free(data);
		}
		if (nz > 1L) {
		  data = (double *)malloc((size_t)nz * sizeof(double));
		  if (data == (double *)NULL) {
			message("ERROR - Unable to allocate z data");
			return(ERROR);
		  }
		  for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  getzD2D(outPtr, x, y, 0L, nx, ny, data, nz);
			  switch (boundary) {
				case null:
				  iirConvolveFinite(data, data, nz, gain, realPoles, np);
				  break;
				case periodic:
				  iirConvolvePeriodic(data, data, nz, gain, realPoles, np);
				  break;
				case mirror:
				  iirConvolveMirror(data, data, nz, gain, realPoles, np);
				  break;
				default:
				  free(data);
				  message("ERROR - Unknown boundary");
				  return(ERROR);
			  }
			  putzD2D(outPtr, x, y, 0L, nx, ny, data, nz);
			}
		  free(data);
		}
		return(!ERROR);
} /* End of directBspline */

/************************************************************************/
static	int			inverseBspline	(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree,
									 enum	bBrands		boundary) {

		double				*inData, *outData;
		double				kernel[7];
		long				nk;
		long				x, y, z;

		switch (degree) {
		  case 0:
		  case 1:
			outPtr = (double *)memcpy(outPtr, inPtr, (size_t)(nx * ny * nz * (long)sizeof(double)));
			return(!ERROR);
		  case 2:
			nk = 3L;
			kernel[0] = 1.0 / 8.0;
			kernel[1] = 6.0 / 8.0;
			kernel[2] = 1.0 / 8.0;
			break;
		  case 3:
			nk = 3L;
			kernel[0] = 1.0 / 6.0;
			kernel[1] = 4.0 / 6.0;
			kernel[2] = 1.0 / 6.0;
			break;
		  case 4:
			nk = 5L;
			kernel[0] = 1.0 / 384.0;
			kernel[1] = 76.0 / 384.0;
			kernel[2] = 230.0 / 384.0;
			kernel[3] = 76.0 / 384.0;
			kernel[4] = 1.0 / 384.0;
			break;
		  case 5:
			nk = 5L;
			kernel[0] = 1.0 / 120.0;
			kernel[1] = 26.0 / 120.0;
			kernel[2] = 66.0 / 120.0;
			kernel[3] = 26.0 / 120.0;
			kernel[4] = 1.0 / 120.0;
			break;
		  case 6:
			nk = 7L;
			kernel[0] = 1.0 / 46080.0;
			kernel[1] = 722.0 / 46080.0;
			kernel[2] = 10543.0 / 46080.0;
			kernel[3] = 23548.0 / 46080.0;
			kernel[4] = 10543.0 / 46080.0;
			kernel[5] = 722.0 / 46080.0;
			kernel[6] = 1.0 / 46080.0;
			break;
		  case 7:
			nk = 7L;
			kernel[0] = 1.0 / 5040.0;
			kernel[1] = 120.0 / 5040.0;
			kernel[2] = 1191.0 / 5040.0;
			kernel[3] = 2416.0 / 5040.0;
			kernel[4] = 1191.0 / 5040.0;
			kernel[5] = 120.0 / 5040.0;
			kernel[6] = 1.0 / 5040.0;
			break;
		  default:
			message("ERROR - Unknown degree");
			return(ERROR);
		}
		outData = outPtr;
		for (z = 0L; (z < nz); z++)
		  for (y = 0L; (y < ny); inPtr += (ptrdiff_t)nx, outData += (ptrdiff_t)nx, y++)
			switch (boundary) {
			  case null:
				firConvolveFinite(inPtr, outData, nx, 0L, kernel, nk);
				break;
			  case periodic:
				firConvolvePeriodic(inPtr, outData, nx, 0L, kernel, nk);
				break;
			  case mirror:
				firConvolveMirror(inPtr, outData, nx, 0L, kernel, nk);
				break;
			  default:
				message("ERROR - Unknown boundary");
				return(ERROR);
			}
		if (ny > 1L) {
		  inData = (double *)malloc((size_t)ny * sizeof(double));
		  if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate y inData");
			return(ERROR);
		  }
		  outData = (double *)malloc((size_t)ny * sizeof(double));
		  if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate y outData");
			return(ERROR);
		  }
		  for (z = 0L; (z < nz); z++)
			for (x = 0L; (x < nx); x++) {
			  getyD2D(outPtr, x, 0L, z, nx, ny, inData, ny);
			  switch (boundary) {
				case null:
				  firConvolveFinite(inData, outData, ny, 0L, kernel, nk);
				  break;
				case periodic:
				  firConvolvePeriodic(inData, outData, ny, 0L, kernel, nk);
				  break;
				case mirror:
				  firConvolveMirror(inData, outData, ny, 0L, kernel, nk);
				  break;
				default:
				  free(outData);
				  free(inData);
				  message("ERROR - Unknown boundary");
				  return(ERROR);
			  }
			  putyD2D(outPtr, x, 0L, z, nx, ny, outData, ny);
			}
		  free(outData);
		  free(inData);
		}
		if (nz > 1L) {
		  inData = (double *)malloc((size_t)nz * sizeof(double));
		  if (inData == (double *)NULL) {
			message("ERROR - Unable to allocate z inData");
			return(ERROR);
		  }
		  outData = (double *)malloc((size_t)nz * sizeof(double));
		  if (outData == (double *)NULL) {
			free(inData);
			message("ERROR - Unable to allocate z outData");
			return(ERROR);
		  }
		  for (y = 0L; (y < ny); y++)
			for (x = 0L; (x < nx); x++) {
			  getzD2D(outPtr, x, y, 0L, nx, ny, inData, nz);
			  switch (boundary) {
				case null:
				  firConvolveFinite(inData, outData, nz, 0L, kernel, nk);
				  break;
				case periodic:
				  firConvolvePeriodic(inData, outData, nz, 0L, kernel, nk);
				  break;
				case mirror:
				  firConvolveMirror(inData, outData, nz, 0L, kernel, nk);
				  break;
				default:
				  free(outData);
				  free(inData);
				  message("ERROR - Unknown boundary");
				  return(ERROR);
			  }
			  putzD2D(outPtr, x, y, 0L, nx, ny, outData, nz);
			}
		  free(outData);
		  free(inData);
		}
		return(!ERROR);
} /* End of inverseBspline */

/************************************************************************/
/* FUNCTION: directBsplineFinite										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may share data storage with the input			*/
/************************************************************************/
int					directBsplineFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (directBspline(inPtr, outPtr, nx, ny, nz, degree, null) == ERROR) {
		  message("ERROR - Unable to perform direct B-spline transform (finite boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of directBsplineFinite */

/************************************************************************/
/* FUNCTION: directBsplineMirror										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may share data storage with the input			*/
/************************************************************************/
int					directBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (directBspline(inPtr, outPtr, nx, ny, nz, degree, mirror) == ERROR) {
		  message("ERROR - Unable to perform direct B-spline transform (mirror boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of directBsplineMirror */

/************************************************************************/
/* FUNCTION: directBsplinePeriodic										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may share data storage with the input			*/
/************************************************************************/
int					directBsplinePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (directBspline(inPtr, outPtr, nx, ny, nz, degree, periodic) == ERROR) {
		  message("ERROR - Unable to perform direct B-spline transform (periodic boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of directBsplinePeriodic */

/************************************************************************/
/* FUNCTION: inverseBsplineFinite										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may NOT share data storage with the input		*/
/************************************************************************/
int					inverseBsplineFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (inverseBspline(inPtr, outPtr, nx, ny, nz, degree, null) == ERROR) {
		  message("ERROR - Unable to perform inverse B-spline transform (finite boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of inverseBsplineFinite */

/************************************************************************/
/* FUNCTION: inverseBsplineMirror										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may NOT share data storage with the input		*/
/************************************************************************/
int					inverseBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (inverseBspline(inPtr, outPtr, nx, ny, nz, degree, mirror) == ERROR) {
		  message("ERROR - Unable to perform inverse B-spline transform (mirror boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of inverseBsplineMirror */

/************************************************************************/
/* FUNCTION: inverseBsplinePeriodic										*/
/*----------------------------------------------------------------------*/
/* Comment:	The output may NOT share data storage with the input		*/
/************************************************************************/
int					inverseBsplinePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree) {

		if (inverseBspline(inPtr, outPtr, nx, ny, nz, degree, periodic) == ERROR) {
		  message("ERROR - Unable to perform inverse B-spline transform (periodic boundary)");
		  return(ERROR);
		}
		return(!ERROR);
} /* End of inverseBsplinePeriodic */
