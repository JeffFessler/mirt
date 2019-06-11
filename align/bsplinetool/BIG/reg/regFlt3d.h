/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	int			regFloat3d		(struct	rParam		*regDataPtr);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* reg0.c																*/
/* reg1.c																*/
/* reg2.c																*/
/************************************************************************/

/* reg0.c */
extern	double		absDeterminant	(double				skew[][3]);

/* reg0.c */
extern	void		adornImage		(float				*outPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 int				clipNx,
									 int				clipNy,
									 int				clipNz,
									 int				clipping,
									 float				backgrnd);

/* reg0.c */
extern	void		computeMskSnr	(double				*mse,
									 double				*snr,
									 float				*inPtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 enum	mBrands		maskCombine,
									 int				nx,
									 int				ny,
									 int				nz);

/* reg1.c */
extern	void		convertOrigin	(struct	fitRec		*fit,
									 double				origin[]);

/* reg0.c */
extern	int			dirMskTransform	(struct fitRec		*fit,
									 float				*inPtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 int				greyRendering,
									 enum	iDegree		interpolation);

/* reg1.c */
extern	void		downscaleFit	(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze);

/* reg0.c */
extern	int			exportFit		(struct	fitRec		*fit);

/* reg0.c */
extern	int			exportSummary	(struct	fitRec		*fit,
									 double				mse,
									 double				snr);

/* reg0.c */
extern	int			fExportFit		(struct	fitRec		*fit,
									 char				*fileName);

/* reg1.c */
extern	void		freePyramids	(float				*pyrHdl1[],
									 float				*pyrHdl2[],
									 float				*mskHdl1[],
									 float				*mskHdl2[],
									 int				levels);

/* reg0.c */
extern	int			importFit		(struct	fitRec		*fit,
									 char				*fitName);

/* reg0.c */
extern	void		initialEstimate	(struct fitRec		*fit,
									 int				nx,
									 int				ny,
									 int				nz);

/* reg0.c */
extern	int			invertFit		(struct	fitRec		*dirFit,
									 struct	fitRec		*invFit);

/* reg0.c */
extern	int			maskFromData	(float				*inPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				sx,
									 double				sy,
									 double				sz);

/* reg2.c */
extern	int			optimize		(struct fitRec		*fit,
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
									 int				nz);

/* reg1.c */
extern	int			pyramid			(float				*inPtr,
									 int				levels,
									 float				*pyrHdl[],
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				nx,
									 int				ny,
									 int				nz,
									 int				zSqueeze,
									 enum	iDegree		interpolation);

/* reg1.c */
extern	void		upscaleFit		(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze);
