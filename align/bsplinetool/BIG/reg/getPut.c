/************************************************************************/
/*																		*/
/* Task:		Extraction and insertion into volumes of line arrays	*/
/*				Full set of conversions (double)-(float)-(short)		*/
/*																		*/
/************************************************************************/

#include		<stddef.h>
#include		<stdlib.h>
#include		<string.h>

#define			MEMORY_CHECK_NO
#undef			MEMORY_CHECK_FULL
#undef			MEMORY_MEDX

#ifdef			MEMORY_CHECK_NO
#endif
#ifdef			MEMORY_CHECK_FULL
#include		"access.h"
#endif
#ifdef			MEMORY_MEDX
#include		"medx.h"
#endif

#include		"getPut.h"

/************************************************************************/
/* FUNCTION: allocLinD													*/
/************************************************************************/
double*				allocLinD		(long				n) {

#ifdef	MEMORY_CHECK_NO
		return((double *)malloc((size_t)n * sizeof(double)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*lin;

		lin = (void *)NULL;
		if (allocateLine(&lin, NIHdouble, n) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (double) line data");
		  return((double *)NULL);
		}
		return((double *)lin);
#endif
#ifdef	MEMORY_MEDX
		return(AllocateLineDouble((int)n));
#endif
} /* End of allocLinD */

/************************************************************************/
/* FUNCTION: allocLinF													*/
/************************************************************************/
float*				allocLinF		(long				n) {

#ifdef	MEMORY_CHECK_NO
		return((float *)malloc((size_t)n * sizeof(float)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*lin;

		lin = (void *)NULL;
		if (allocateLine(&lin, NIHfloat, n) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (float) line data");
		  return((float *)NULL);
		}
		return((float *)lin);
#endif
#ifdef	MEMORY_MEDX
		return(AllocateLineFloat((int)n));
#endif
} /* End of allocLinF */

/************************************************************************/
/* FUNCTION: allocLinS													*/
/************************************************************************/
short*				allocLinS		(long				n) {

#ifdef	MEMORY_CHECK_NO
		return((short *)malloc((size_t)n * sizeof(short)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*lin;

		lin = (void *)NULL;
		if (allocateLine(&lin, NIHshort, n) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (short) line data");
		  return((short *)NULL);
		}
		return((short *)lin);
#endif
#ifdef	MEMORY_MEDX
		return(AllocateLineShort((int)n));
#endif
} /* End of allocLinS */

/************************************************************************/
/* FUNCTION: allocVolD													*/
/************************************************************************/
double*				allocVolD		(long				nx,
									 long				ny,
									 long				nz) {

#ifdef	MEMORY_CHECK_NO
		return((double *)malloc((size_t)(nx * ny * nz) * sizeof(double)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*vol;

		vol = (void *)NULL;
		if (allocateVolume(&vol, NIHdouble, nx, ny, nz) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (double) volume data");
		  return((double *)NULL);
		}
		return((double *)vol);
#endif
#ifdef	MEMORY_MEDX
		return((double *)CreateVolume(64, 0, (int)nx, (int)ny, (int)nz));
#endif
} /* End of allocVolD */

/************************************************************************/
/* FUNCTION: allocVolF													*/
/************************************************************************/
float*				allocVolF		(long				nx,
									 long				ny,
									 long				nz) {

#ifdef	MEMORY_CHECK_NO
		return((float *)malloc((size_t)(nx * ny * nz) * sizeof(float)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*vol;

		vol = (void *)NULL;
		if (allocateVolume(&vol, NIHfloat, nx, ny, nz) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (float) volume data");
		  return((float *)NULL);
		}
		return((float *)vol);
#endif
#ifdef	MEMORY_MEDX
		return((float *)CreateVolume(32, 0, (int)nx, (int)ny, (int)nz));
#endif
} /* End of allocVolF */

/************************************************************************/
/* FUNCTION: allocVolS													*/
/************************************************************************/
short*				allocVolS		(long				nx,
									 long				ny,
									 long				nz) {

#ifdef	MEMORY_CHECK_NO
		return((short *)malloc((size_t)(nx * ny * nz) * sizeof(short)));
#endif
#ifdef	MEMORY_CHECK_FULL
		void				*vol;

		vol = (void *)NULL;
		if (allocateVolume(&vol, NIHshort, nx, ny, nz) == ERROR) {
		  errorReport("ERROR - Not enough memory to allocate a (short) volume data");
		  return((short *)NULL);
		}
		return((short *)vol);
#endif
#ifdef	MEMORY_MEDX
		return((short *)CreateVolume(16, 0, (int)nx, (int)ny, (int)nz));
#endif
} /* End of allocVolS */

/************************************************************************/
/* FUNCTION: freeLinD													*/
/************************************************************************/
void				freeLinD		(double				*line) {

#ifdef	MEMORY_CHECK_NO
		free(line);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)line, NIHdouble, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (double) line data");
		else
		  freeLine((struct NIHline *)line);
#endif
#ifdef	MEMORY_MEDX
		FreeLine(line);
#endif
} /* End of freeLinD */

/************************************************************************/
/* FUNCTION: freeLinF													*/
/************************************************************************/
void				freeLinF		(float				*line) {

#ifdef	MEMORY_CHECK_NO
		free(line);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)line, NIHfloat, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (float) line data");
		else
		  freeLine((struct NIHline *)line);
#endif
#ifdef	MEMORY_MEDX
		FreeLine(line);
#endif
} /* End of freeLinF */

/************************************************************************/
/* FUNCTION: freeLinS													*/
/************************************************************************/
void				freeLinS		(short				*line) {

#ifdef	MEMORY_CHECK_NO
		free(line);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)line, NIHshort, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (short) line data");
		else
		  freeLine((struct NIHline *)line);
#endif
#ifdef	MEMORY_MEDX
		FreeLine(line);
#endif
} /* End of freeLinS */

/************************************************************************/
/* FUNCTION: freeVolD													*/
/************************************************************************/
void				freeVolD		(double				*volume) {

#ifdef	MEMORY_CHECK_NO
		free(volume);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volume, NIHdouble, (long *)NULL, (long *)NULL,
		  (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (double) volume data");
		else
		  freeVolume((struct NIHvolume *)volume);
#endif
#ifdef	MEMORY_MEDX
		FreeVolume((Canvas *)volume);
#endif
} /* End of freeVolD */

/************************************************************************/
/* FUNCTION: freeVolF													*/
/************************************************************************/
void				freeVolF		(float				*volume) {

#ifdef	MEMORY_CHECK_NO
		free(volume);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volume, NIHfloat, (long *)NULL, (long *)NULL,
		  (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (float) volume data");
		else
		  freeVolume((struct NIHvolume *)volume);
#endif
#ifdef	MEMORY_MEDX
		FreeVolume((Canvas *)volume);
#endif
} /* End of freeVolF */

/************************************************************************/
/* FUNCTION: freeVolS													*/
/************************************************************************/
void				freeVolS		(short				*volume) {

#ifdef	MEMORY_CHECK_NO
		free(volume);
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volume, NIHshort, (long *)NULL, (long *)NULL,
		  (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to free a (short) volume data");
		else
		  freeVolume((struct NIHvolume *)volume);
#endif
#ifdef	MEMORY_MEDX
		FreeVolume((Canvas *)volume);
#endif
} /* End of freeVolS */

/************************************************************************/
/* FUNCTION: getPxlD2D													*/
/************************************************************************/
double				getPxlD2D		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return(*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (double) volume data");
		  return(0.0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0);
		}
		return(pxl.data.d);
#endif
#ifdef	MEMORY_MEDX
		return(GetDoubleVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlD2D */

/************************************************************************/
/* FUNCTION: getPxlD2F													*/
/************************************************************************/
float				getPxlD2F		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((float)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (double) volume data");
		  return(0.0F);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0F);
		}
		return(pxl.data.f);
#endif
#ifdef	MEMORY_MEDX
		return(GetFloatVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlD2F */

/************************************************************************/
/* FUNCTION: getPxlD2S													*/
/************************************************************************/
float				getPxlD2S		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((short)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (double) volume data");
		  return((short)0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return((short)0);
		}
		return(pxl.data.s);
#endif
#ifdef	MEMORY_MEDX
		return(GetShortVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlD2S */

/************************************************************************/
/* FUNCTION: getPxlF2D													*/
/************************************************************************/
double				getPxlF2D		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((double)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (float) volume data");
		  return(0.0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0);
		}
		return(pxl.data.d);
#endif
#ifdef	MEMORY_MEDX
		return(GetDoubleVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlF2D */

/************************************************************************/
/* FUNCTION: getPxlF2F													*/
/************************************************************************/
float				getPxlF2F		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return(*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (float) volume data");
		  return(0.0F);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0F);
		}
		return(pxl.data.f);
#endif
#ifdef	MEMORY_MEDX
		return(GetFloatVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlF2F */

/************************************************************************/
/* FUNCTION: getPxlF2S													*/
/************************************************************************/
short				getPxlF2S		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((double)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (float) volume data");
		  return((short)0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return((short)0);
		}
		return(pxl.data.s);
#endif
#ifdef	MEMORY_MEDX
		return(GetShortVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlF2S */

/************************************************************************/
/* FUNCTION: getPxlS2D													*/
/************************************************************************/
double				getPxlS2D		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((double)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (short) volume data");
		  return(0.0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0);
		}
		return(pxl.data.d);
#endif
#ifdef	MEMORY_MEDX
		return(GetDoubleVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlS2D */

/************************************************************************/
/* FUNCTION: getPxlS2F													*/
/************************************************************************/
float				getPxlS2F		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return((float)*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (short) volume data");
		  return(0.0F);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return(0.0F);
		}
		return(pxl.data.f);
#endif
#ifdef	MEMORY_MEDX
		return(GetFloatVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlS2F */

/************************************************************************/
/* FUNCTION: getPxlS2S													*/
/************************************************************************/
short				getPxlS2S		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny) {

#ifdef	MEMORY_CHECK_NO
		return(*(volumeSrc + (ptrdiff_t)(x + nx * (y + ny * z))));
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR) {
		  errorReport("ERROR - Unable to access a (short) volume data");
		  return((short)0);
		}
		if (getPixel((struct NIHvolume *)volumeSrc, &pxl, x, y, z) == ERROR) {
		  errorReport("ERROR - Unable to get a pixel data");
		  return((short)0);
		}
		return(pxl.data.s);
#endif
#ifdef	MEMORY_MEDX
		return(GetShortVoxel((Canvas *)volumeSrc, (int)x, (int)y, (int)z));
#endif
} /* End of getPxlS2S */

/************************************************************************/
/* FUNCTION: getxD2D													*/
/************************************************************************/
void				getxD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		rowDst = (double *)memcpy(rowDst, volumeSrc, (size_t)n * sizeof(double));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowDouble((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxD2D */

/************************************************************************/
/* FUNCTION: getxD2F													*/
/************************************************************************/
void				getxD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (float)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowFloat((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxD2F */

/************************************************************************/
/* FUNCTION: getxD2S													*/
/************************************************************************/
void				getxD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (short)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowShort((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxD2S */

/************************************************************************/
/* FUNCTION: getxF2D													*/
/************************************************************************/
void				getxF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (double)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowDouble((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxF2D */

/************************************************************************/
/* FUNCTION: getxF2F													*/
/************************************************************************/
void				getxF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		rowDst = (float *)memcpy(rowDst, volumeSrc, (size_t)n * sizeof(float));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowFloat((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxF2F */

/************************************************************************/
/* FUNCTION: getxF2S													*/
/************************************************************************/
void				getxF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (short)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowShort((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxF2S */

/************************************************************************/
/* FUNCTION: getxS2D													*/
/************************************************************************/
void				getxS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (double)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowDouble((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxS2D */

/************************************************************************/
/* FUNCTION: getxS2F													*/
/************************************************************************/
void				getxS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowDst + (ptrdiff_t)n;
		while (rowDst < end)
		  *rowDst++ = (float)*volumeSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowFloat((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxS2F */

/************************************************************************/
/* FUNCTION: getxS2S													*/
/************************************************************************/
void				getxS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		rowDst = (short *)memcpy(rowDst, volumeSrc, (size_t)n * sizeof(short));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)rowDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getRow((struct NIHvolume *)volumeSrc, (struct NIHline *)rowDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetRowShort((Canvas *)volumeSrc, rowDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getxS2S */

/************************************************************************/
/* FUNCTION: getyD2D													*/
/************************************************************************/
void				getyD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnDouble((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyD2D */

/************************************************************************/
/* FUNCTION: getyD2F													*/
/************************************************************************/
void				getyD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (float)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnFloat((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyD2F */

/************************************************************************/
/* FUNCTION: getyD2S													*/
/************************************************************************/
void				getyD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (short)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnShort((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyD2S */

/************************************************************************/
/* FUNCTION: getyF2D													*/
/************************************************************************/
void				getyF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (double)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnDouble((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyF2D */

/************************************************************************/
/* FUNCTION: getyF2F													*/
/************************************************************************/
void				getyF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnFloat((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyF2F */

/************************************************************************/
/* FUNCTION: getyF2S													*/
/************************************************************************/
void				getyF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (short)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnShort((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyF2S */

/************************************************************************/
/* FUNCTION: getyS2D													*/
/************************************************************************/
void				getyS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (double)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnDouble((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyS2D */

/************************************************************************/
/* FUNCTION: getyS2F													*/
/************************************************************************/
void				getyS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = (float)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnFloat((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyS2F */

/************************************************************************/
/* FUNCTION: getyS2S													*/
/************************************************************************/
void				getyS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeSrc += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnDst + (ptrdiff_t)n; (columnDst < end); volumeSrc += (ptrdiff_t)nx)
		  *columnDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)columnDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)columnDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to get a line data");
#endif
#ifdef	MEMORY_MEDX
		GetColumnShort((Canvas *)volumeSrc, columnDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getyS2S */

/************************************************************************/
/* FUNCTION: getzD2D													*/
/************************************************************************/
void				getzD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarDouble((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzD2D */

/************************************************************************/
/* FUNCTION: getzD2F													*/
/************************************************************************/
void				getzD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (float)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarFloat((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzD2F */

/************************************************************************/
/* FUNCTION: getzD2S													*/
/************************************************************************/
void				getzD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (short)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarShort((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzD2S */

/************************************************************************/
/* FUNCTION: getzF2D													*/
/************************************************************************/
void				getzF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (double)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarDouble((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzF2D */

/************************************************************************/
/* FUNCTION: getzF2F													*/
/************************************************************************/
void				getzF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarFloat((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzF2F */

/************************************************************************/
/* FUNCTION: getzF2S													*/
/************************************************************************/
void				getzF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (short)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarShort((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzF2S */

/************************************************************************/
/* FUNCTION: getzS2D													*/
/************************************************************************/
void				getzS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (double)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarDouble((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzS2D */

/************************************************************************/
/* FUNCTION: getzS2F													*/
/************************************************************************/
void				getzS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = (float)*volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarFloat((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzS2F */

/************************************************************************/
/* FUNCTION: getzS2S													*/
/************************************************************************/
void				getzS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeSrc += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarDst + (ptrdiff_t)n; (pillarDst < end); volumeSrc += (ptrdiff_t)nxny)
		  *pillarDst++ = *volumeSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkVol((struct NIHvolume *)volumeSrc, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (checkLin((struct NIHline *)pillarDst, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (getColumn((struct NIHvolume *)volumeSrc, (struct NIHline *)pillarDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to access a line data");
#endif
#ifdef	MEMORY_MEDX
		GetPillarShort((Canvas *)volumeSrc, pillarDst, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of getzS2S */

/************************************************************************/
/* FUNCTION: putPxlD2D													*/
/************************************************************************/
void				putPxlD2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.d = pixelSrc;
		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutDoubleVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlD2D */

/************************************************************************/
/* FUNCTION: putPxlD2F													*/
/************************************************************************/
void				putPxlD2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (float)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.d = pixelSrc;
		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutDoubleVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlD2F */

/************************************************************************/
/* FUNCTION: putPxlD2S													*/
/************************************************************************/
void				putPxlD2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (short)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.d = pixelSrc;
		pxl.dataType = NIHdouble;
		if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutDoubleVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlD2S */

/************************************************************************/
/* FUNCTION: putPxlF2D													*/
/************************************************************************/
void				putPxlF2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (double)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.f = pixelSrc;
		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutFloatVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlF2D */

/************************************************************************/
/* FUNCTION: putPxlF2F													*/
/************************************************************************/
void				putPxlF2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.f = pixelSrc;
		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutFloatVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlF2F */

/************************************************************************/
/* FUNCTION: putPxlF2S													*/
/************************************************************************/
void				putPxlF2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (short)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.f = pixelSrc;
		pxl.dataType = NIHfloat;
		if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutFloatVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlF2S */

/************************************************************************/
/* FUNCTION: putPxlS2D													*/
/************************************************************************/
void				putPxlS2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (double)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.s = pixelSrc;
		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutShortVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlS2D */

/************************************************************************/
/* FUNCTION: putPxlS2F													*/
/************************************************************************/
void				putPxlS2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = (float)pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.s = pixelSrc;
		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutShortVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlS2F */

/************************************************************************/
/* FUNCTION: putPxlS2S													*/
/************************************************************************/
void				putPxlS2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc) {

#ifdef	MEMORY_CHECK_NO
		*(volumeDst + (ptrdiff_t)(x + nx * (y + ny * z))) = pixelSrc;
#endif
#ifdef	MEMORY_CHECK_FULL
		struct	NIHpixel	pxl;

		pxl.data.s = pixelSrc;
		pxl.dataType = NIHshort;
		if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putPixel(&pxl, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a pixel data");
#endif
#ifdef	MEMORY_MEDX
		PutShortVoxel((Canvas *)volumeDst, (int)x, (int)y, (int)z, pixelSrc);
#endif
} /* End of putPxlS2S */

/************************************************************************/
/* FUNCTION: putxD2D													*/
/************************************************************************/
void				putxD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		volumeDst = (double *)memcpy(volumeDst, rowSrc, (size_t)n * sizeof(double));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowDouble((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxD2D */

/************************************************************************/
/* FUNCTION: putxD2F													*/
/************************************************************************/
void				putxD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (float)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowDouble((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxD2F */

/************************************************************************/
/* FUNCTION: putxD2S													*/
/************************************************************************/
void				putxD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (short)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowDouble((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxD2S */

/************************************************************************/
/* FUNCTION: putxF2D													*/
/************************************************************************/
void				putxF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (double)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowFloat((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxF2D */

/************************************************************************/
/* FUNCTION: putxF2F													*/
/************************************************************************/
void				putxF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		volumeDst = (float *)memcpy(volumeDst, rowSrc, (size_t)n * sizeof(float));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowFloat((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxF2F */

/************************************************************************/
/* FUNCTION: putxF2S													*/
/************************************************************************/
void				putxF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (short)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowFloat((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxF2S */

/************************************************************************/
/* FUNCTION: putxS2D													*/
/************************************************************************/
void				putxS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (double)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowShort((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxS2D */

/************************************************************************/
/* FUNCTION: putxS2F													*/
/************************************************************************/
void				putxS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		end = rowSrc + (ptrdiff_t)n;
		while (rowSrc < end)
		  *volumeDst++ = (float)*rowSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowShort((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxS2F */

/************************************************************************/
/* FUNCTION: putxS2S													*/
/************************************************************************/
void				putxS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		volumeDst = (short *)memcpy(volumeDst, rowSrc, (size_t)n * sizeof(short));
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)rowSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)rowSrc, (struct NIHvolume *)volumeDst, x, y, z) == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutRowShort((Canvas *)volumeDst, rowSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putxS2S */

/************************************************************************/
/* FUNCTION: putyD2D													*/
/************************************************************************/
void				putyD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = *columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnDouble((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyD2D */

/************************************************************************/
/* FUNCTION: putyD2F													*/
/************************************************************************/
void				putyD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (float)*columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnDouble((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyD2F */

/************************************************************************/
/* FUNCTION: putyD2S													*/
/************************************************************************/
void				putyD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (short)*columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnDouble((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyD2S */

/************************************************************************/
/* FUNCTION: putyF2D													*/
/************************************************************************/
void				putyF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (double)*columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnFloat((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyF2D */

/************************************************************************/
/* FUNCTION: putyF2F													*/
/************************************************************************/
void				putyF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = *columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnFloat((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyF2F */

/************************************************************************/
/* FUNCTION: putyF2S													*/
/************************************************************************/
void				putyF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = *columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnFloat((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyF2S */

/************************************************************************/
/* FUNCTION: putyS2D													*/
/************************************************************************/
void				putyS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (double)*columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnShort((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyS2D */

/************************************************************************/
/* FUNCTION: putyS2F													*/
/************************************************************************/
void				putyS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = (float)*columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnShort((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyS2F */

/************************************************************************/
/* FUNCTION: putyS2S													*/
/************************************************************************/
void				putyS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;

		volumeDst += (ptrdiff_t)(x + nx * (y + ny * z));
		for (end = columnSrc + (ptrdiff_t)n; (columnSrc < end); volumeDst += (ptrdiff_t)nx)
		  *volumeDst = *columnSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)columnSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)columnSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutColumnShort((Canvas *)volumeDst, columnSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putyS2S */

/************************************************************************/
/* FUNCTION: putzD2D													*/
/************************************************************************/
void				putzD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = *pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarDouble((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzD2D */

/************************************************************************/
/* FUNCTION: putzD2F													*/
/************************************************************************/
void				putzD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (float)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarDouble((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzD2F */

/************************************************************************/
/* FUNCTION: putzD2S													*/
/************************************************************************/
void				putzD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		double				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (short)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHdouble, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (double) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarDouble((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzD2S */

/************************************************************************/
/* FUNCTION: putzF2D													*/
/************************************************************************/
void				putzF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (double)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarFloat((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzF2D */

/************************************************************************/
/* FUNCTION: putzF2F													*/
/************************************************************************/
void				putzF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = *pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarFloat((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzF2F */

/************************************************************************/
/* FUNCTION: putzF2S													*/
/************************************************************************/
void				putzF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		float				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (short)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHfloat, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (float) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarFloat((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzF2S */

/************************************************************************/
/* FUNCTION: putzS2D													*/
/************************************************************************/
void				putzS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (double)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHdouble, &nx, &ny, (long *)NULL)
		  == ERROR)
		  errorReport("ERROR - Unable to access a (double) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarShort((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzS2D */

/************************************************************************/
/* FUNCTION: putzS2F													*/
/************************************************************************/
void				putzS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = (float)*pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHfloat, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (float) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarShort((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzS2F */

/************************************************************************/
/* FUNCTION: putzS2S													*/
/************************************************************************/
void				putzS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n) {

#ifdef	MEMORY_CHECK_NO
		short				*end;
		long				nxny = nx * ny;

		volumeDst += (ptrdiff_t)(x + nx * y + nxny * z);
		for (end = pillarSrc + (ptrdiff_t)n; (pillarSrc < end); volumeDst += (ptrdiff_t)nxny)
		  *volumeDst = *pillarSrc++;
#endif
#ifdef	MEMORY_CHECK_FULL
		if (checkLin((struct NIHline *)pillarSrc, NIHshort, &n) == ERROR)
		  errorReport("ERROR - Unable to access a (short) line data");
		else if (checkVol((struct NIHvolume *)volumeDst, NIHshort, &nx, &ny, (long *)NULL) == ERROR)
		  errorReport("ERROR - Unable to access a (short) volume data");
		else if (putRow((struct NIHline *)pillarSrc, (struct NIHvolume *)volumeDst, x, y, z)
		  == ERROR)
		  errorReport("ERROR - Unable to put a line data");
#endif
#ifdef	MEMORY_MEDX
		PutPillarShort((Canvas *)volumeDst, pillarSrc, (int)x, (int)y, (int)z, (int)n);
#endif
} /* End of putzS2S */
