/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	double*		allocLinD		(long				n);
extern	float*		allocLinF		(long				n);
extern	short*		allocLinS		(long				n);
extern	double*		allocVolD		(long				nx,
									 long				ny,
									 long				nz);
extern	float*		allocVolF		(long				nx,
									 long				ny,
									 long				nz);
extern	short*		allocVolS		(long				nx,
									 long				ny,
									 long				nz);
extern	void		freeLinD		(double				*line);
extern	void		freeLinF		(float				*line);
extern	void		freeLinS		(short				*line);
extern	void		freeVolD		(double				*volume);
extern	void		freeVolF		(float				*volume);
extern	void		freeVolS		(short				*volume);
extern	double		getPxlD2D		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	float		getPxlD2F		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	float		getPxlD2S		(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	double		getPxlF2D		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	float		getPxlF2F		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	short		getPxlF2S		(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	double		getPxlS2D		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	float		getPxlS2F		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	short		getPxlS2S		(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny);
extern	void		getxD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n);
extern	void		getxD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n);
extern	void		getxD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n);
extern	void		getxF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n);
extern	void		getxF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n);
extern	void		getxF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n);
extern	void		getxS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowDst,
									 long				n);
extern	void		getxS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowDst,
									 long				n);
extern	void		getxS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowDst,
									 long				n);
extern	void		getyD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n);
extern	void		getyD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n);
extern	void		getyD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n);
extern	void		getyF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n);
extern	void		getyF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n);
extern	void		getyF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n);
extern	void		getyS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnDst,
									 long				n);
extern	void		getyS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnDst,
									 long				n);
extern	void		getyS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnDst,
									 long				n);
extern	void		getzD2D			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n);
extern	void		getzD2F			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n);
extern	void		getzD2S			(double				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n);
extern	void		getzF2D			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n);
extern	void		getzF2F			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n);
extern	void		getzF2S			(float				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n);
extern	void		getzS2D			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarDst,
									 long				n);
extern	void		getzS2F			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarDst,
									 long				n);
extern	void		getzS2S			(short				*volumeSrc,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarDst,
									 long				n);
extern	void		putPxlD2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc);
extern	void		putPxlD2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc);
extern	void		putPxlD2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				pixelSrc);
extern	void		putPxlF2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc);
extern	void		putPxlF2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc);
extern	void		putPxlF2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				pixelSrc);
extern	void		putPxlS2D		(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc);
extern	void		putPxlS2F		(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc);
extern	void		putPxlS2S		(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				pixelSrc);
extern	void		putxD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n);
extern	void		putxD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n);
extern	void		putxD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*rowSrc,
									 long				n);
extern	void		putxF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n);
extern	void		putxF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n);
extern	void		putxF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*rowSrc,
									 long				n);
extern	void		putxS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n);
extern	void		putxS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n);
extern	void		putxS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*rowSrc,
									 long				n);
extern	void		putyD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n);
extern	void		putyD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n);
extern	void		putyD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*columnSrc,
									 long				n);
extern	void		putyF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n);
extern	void		putyF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n);
extern	void		putyF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*columnSrc,
									 long				n);
extern	void		putyS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n);
extern	void		putyS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n);
extern	void		putyS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*columnSrc,
									 long				n);
extern	void		putzD2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n);
extern	void		putzD2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n);
extern	void		putzD2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillarSrc,
									 long				n);
extern	void		putzF2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n);
extern	void		putzF2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n);
extern	void		putzF2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 float				*pillarSrc,
									 long				n);
extern	void		putzS2D			(double				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n);
extern	void		putzS2F			(float				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n);
extern	void		putzS2S			(short				*volumeDst,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 short				*pillarSrc,
									 long				n);

/************************************************************************/
/* MACROS																*/
/************************************************************************/

#undef				GETPXLD2D
#ifdef	MEMORY_CHECK_NO
#define				GETPXLD2D(		 V, X, Y, Z, NX, NY) ( \
		*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) \
) /* End of macro GETPXLD2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLD2D(		 V, X, Y, Z, NX, NY) ( \
		getPxlD2D((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLD2D */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLD2D(		 V, X, Y, Z, NX, NY) ( \
		(GetDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLD2D */
#endif

#undef				GETPXLD2F
#ifdef	MEMORY_CHECK_NO
#define				GETPXLD2F(		 V, X, Y, Z, NX, NY) ( \
		(float)*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) \
) /* End of macro GETPXLD2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLD2F(		 V, X, Y, Z, NX, NY) ( \
		getPxlD2F((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLD2F */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLD2F(		 V, X, Y, Z, NX, NY) ( \
		(GetFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLD2F */
#endif

#undef				GETPXLD2S
#ifdef	MEMORY_CHECK_NO
#define				GETPXLD2S(		 V, X, Y, Z, NX, NY) ( \
		(short)*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) \
) /* End of macro GETPXLD2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLD2S(		 V, X, Y, Z, NX, NY) ( \
		getPxlD2S((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLD2S */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLD2S(		 V, X, Y, Z, NX, NY) ( \
		(GetShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLD2S */
#endif

#undef				GETPXLF2D
#ifdef	MEMORY_CHECK_NO
#define				GETPXLF2D(		 V, X, Y, Z, NX, NY) ( \
		(double)(*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z))))) \
) /* End of macro GETPXLF2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLF2D(		 V, X, Y, Z, NX, NY) ( \
		getPxlF2D((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLF2D */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLF2D(		 V, X, Y, Z, NX, NY) ( \
		(GetDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLF2D */
#endif

#undef				GETPXLF2F
#ifdef	MEMORY_CHECK_NO
#define				GETPXLF2F(		 V, X, Y, Z, NX, NY) ( \
		*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) \
) /* End of macro GETPXLF2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLF2F(		 V, X, Y, Z, NX, NY) ( \
		getPxlF2F((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLF2F */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLF2F(		 V, X, Y, Z, NX, NY) ( \
		(GetFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLF2F */
#endif

#undef				GETPXLF2S
#ifdef	MEMORY_CHECK_NO
#define				GETPXLF2S(		 V, X, Y, Z, NX, NY) ( \
		(short)(*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z))))) \
) /* End of macro GETPXLF2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLF2S(		 V, X, Y, Z, NX, NY) ( \
		getPxlF2S((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLF2S */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLF2S(		 V, X, Y, Z, NX, NY) ( \
		(GetShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLF2S */
#endif

#undef				GETPXLS2D
#ifdef	MEMORY_CHECK_NO
#define				GETPXLS2D(		 V, X, Y, Z, NX, NY) ( \
		(double)(*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z))))) \
) /* End of macro GETPXLS2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLS2D(		 V, X, Y, Z, NX, NY) ( \
		getPxlS2D((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLS2D */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLS2D(		 V, X, Y, Z, NX, NY) ( \
		(GetDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLS2D */
#endif

#undef				GETPXLS2F
#ifdef	MEMORY_CHECK_NO
#define				GETPXLS2F(		 V, X, Y, Z, NX, NY) ( \
		(float)(*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z))))) \
) /* End of macro GETPXLS2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLS2F(		 V, X, Y, Z, NX, NY) ( \
		getPxlS2F((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLS2F */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLS2F(		 V, X, Y, Z, NX, NY) ( \
		(GetFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLS2F */
#endif

#undef				GETPXLS2S
#ifdef	MEMORY_CHECK_NO
#define				GETPXLS2S(		 V, X, Y, Z, NX, NY) ( \
		*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) \
) /* End of macro GETPXLS2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				GETPXLS2S(		 V, X, Y, Z, NX, NY) ( \
		getPxlS2S((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY)) \
) /* End of macro GETPXLS2S */
#endif
#ifdef	MEMORY_MEDX
#define				GETPXLS2S(		 V, X, Y, Z, NX, NY) ( \
		(GetShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z))) \
) /* End of macro GETPXLS2S */
#endif

#undef				PUTPXLD2D
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLD2D(		 V, X, Y, Z, NX, NY, P) ( \
		*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (double)(P) \
) /* End of macro PUTPXLD2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLD2D(		 V, X, Y, Z, NX, NY) ( \
		putPxlD2D((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (double)P) \
) /* End of macro PUTPXLD2D */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLD2D(		 V, X, Y, Z, NX, NY) ( \
		(PutDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (double)P)) \
) /* End of macro PUTPXLD2D */
#endif

#undef				PUTPXLD2F
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLD2F(		 V, X, Y, Z, NX, NY, P) ( \
		*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (float)((double)(P)) \
) /* End of macro PUTPXLD2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLD2F(		 V, X, Y, Z, NX, NY) ( \
		putPxlD2F((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (double)P) \
) /* End of macro PUTPXLD2F */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLD2F(		 V, X, Y, Z, NX, NY) ( \
		(PutDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (double)P)) \
) /* End of macro PUTPXLD2F */
#endif

#undef				PUTPXLD2S
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLD2S(		 V, X, Y, Z, NX, NY, P) ( \
		*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (short)((double)(P)) \
) /* End of macro PUTPXLD2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLD2S(		 V, X, Y, Z, NX, NY) ( \
		putPxlD2S((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (double)P) \
) /* End of macro PUTPXLD2S */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLD2S(		 V, X, Y, Z, NX, NY) ( \
		(PutDoubleVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (double)P)) \
) /* End of macro PUTPXLD2S */
#endif

#undef				PUTPXLF2D
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLF2D(		 V, X, Y, Z, NX, NY, P) ( \
		*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (double)((float)(P)) \
) /* End of macro PUTPXLF2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLF2D(		 V, X, Y, Z, NX, NY) ( \
		putPxlF2D((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (float)P) \
) /* End of macro PUTPXLF2D */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLF2D(		 V, X, Y, Z, NX, NY) ( \
		(PutFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (float)P)) \
) /* End of macro PUTPXLF2D */
#endif

#undef				PUTPXLF2F
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLF2F(		 V, X, Y, Z, NX, NY, P) ( \
		*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (float)(P) \
) /* End of macro PUTPXLF2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLF2F(		 V, X, Y, Z, NX, NY) ( \
		putPxlF2F((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (float)P) \
) /* End of macro PUTPXLF2F */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLF2F(		 V, X, Y, Z, NX, NY) ( \
		(PutFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (float)P)) \
) /* End of macro PUTPXLF2F */
#endif

#undef				PUTPXLF2S
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLF2S(		 V, X, Y, Z, NX, NY, P) ( \
		*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (short)((float)(P)) \
) /* End of macro PUTPXLF2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLF2S(		 V, X, Y, Z, NX, NY) ( \
		putPxlF2S((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (float)P) \
) /* End of macro PUTPXLF2S */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLF2S(		 V, X, Y, Z, NX, NY) ( \
		(PutFloatVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (float)P)) \
) /* End of macro PUTPXLF2S */
#endif

#undef				PUTPXLS2D
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLS2D(		 V, X, Y, Z, NX, NY, P) ( \
		*((double *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (double)((short)(P)) \
) /* End of macro PUTPXLS2D */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLS2D(		 V, X, Y, Z, NX, NY) ( \
		putPxlS2D((double *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (short)P) \
) /* End of macro PUTPXLS2D */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLS2D(		 V, X, Y, Z, NX, NY) ( \
		(PutShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (short)P)) \
) /* End of macro PUTPXLS2D */
#endif

#undef				PUTPXLS2F
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLS2F(		 V, X, Y, Z, NX, NY, P) ( \
		*((float *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (float)((short)(P)) \
) /* End of macro PUTPXLS2F */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLS2F(		 V, X, Y, Z, NX, NY) ( \
		putPxlS2F((float *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (short)P) \
) /* End of macro PUTPXLS2F */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLS2F(		 V, X, Y, Z, NX, NY) ( \
		(PutShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (short)P)) \
) /* End of macro PUTPXLS2F */
#endif

#undef				PUTPXLS2S
#ifdef	MEMORY_CHECK_NO
#define				PUTPXLS2S(		 V, X, Y, Z, NX, NY, P) ( \
		*((short *)(V) + (ptrdiff_t)((long)(X) + (long)(NX) * ((long)(Y) + (long)(NY) \
		  * (long)(Z)))) = (short)(P) \
) /* End of macro PUTPXLS2S */
#endif
#ifdef	MEMORY_CHECK_FULL
#define				PUTPXLS2S(		 V, X, Y, Z, NX, NY) ( \
		putPxlS2S((short *)(V), (long)(X), (long)(Y), (long)(Z), (long)(NX), (long)(NY), \
		  (short)P) \
) /* End of macro PUTPXLS2S */
#endif
#ifdef	MEMORY_MEDX
#define				PUTPXLS2S(		 V, X, Y, Z, NX, NY) ( \
		(PutShortVoxel((Canvas *)(V), (int)(X), (int)(Y), (int)(Z), (short)P)) \
) /* End of macro PUTPXLS2S */
#endif

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* access.h																*/
/* configs.h															*/
/* medx.h																*/
/************************************************************************/
