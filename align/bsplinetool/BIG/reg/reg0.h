/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	double		absDeterminant	(double				skew[][3]);
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
extern	void		combineMasks	(float				*inMsk,
									 float				*outMsk,
									 enum	mBrands		maskCombine,
									 int				nx,
									 int				ny,
									 int				nz);
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
extern	void		computeRotChisq	(double				beta[],
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
									 double				zO);
extern	void		computeSkewChisq(double				beta[],
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
									 double				zO);
extern	void		copyClip		(float				*inPtr,
									 int				nxIn,
									 int				nyIn,
									 int				nzIn,
									 float				*outPtr,
									 int				nxOut,
									 int				nyOut,
									 int				nzOut);
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
extern	void		eulerRotation	(struct	fitRec		*fit);
extern	int			exportFit		(struct	fitRec		*fit);
extern	int			exportSummary	(struct	fitRec		*fit,
									 double				mse,
									 double				snr);
extern	int			fExportFit		(struct	fitRec		*fit,
									 char				*fileName);
extern	int			importFit		(struct	fitRec		*fit,
									 char				*fitName);
extern	void		initialEstimate	(struct fitRec		*fit,
									 int				nx,
									 int				ny,
									 int				nz);
extern	int			invertFit		(struct	fitRec		*dirFit,
									 struct	fitRec		*invFit);
extern	int			invTransform	(struct fitRec		*fit,
									 float				*inPtr,
									 double				*splinePtr,
									 float				*outPtr,
									 float				*inMsk,
									 float				*outMsk,
									 int				nx,
									 int				ny,
									 int				nz,
									 enum	iDegree		interpolation);
extern	int			maskFromData	(float				*inPtr,
									 float				*mskPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 double				sx,
									 double				sy,
									 double				sz);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* BsplnTrf.c															*/
/* BsplnWgt.c															*/
/* getPut.c																*/
/* intp1.c																*/
/* quant.c																*/
/* reg1.c																*/
/************************************************************************/

/* intp1.c */
extern	long		assign_double_command_var
									(char				*name,
									 double				dvalue);

/* BsplnWgt.c */
extern	double		BsplnWght		(int				degree,
									 long				i,
									 double				x);

/* BsplnWgt.c */
extern	double		BsplnWght1		(long				i,
									 double				x);

/* BsplnWgt.c */
extern	double		BsplnWght3		(long				i,
									 double				x);

/* reg1.c */
extern	void		convertOrigin	(struct	fitRec		*fit,
									 double				origin[]);

/* BsplnTrf.c */
extern	int			directBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);

/* BsplnWgt.c */
extern	long		fold			(long				i,
									 long				n);

/* getPut.c */
extern	void		getxF2D			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		getyF2D			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		getzF2D			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* quant.c */
extern	int			histQuant		(struct	qParam		*quantDataPtr);

/* convolve.c */
extern	void		iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);

/* getPut.c */
extern	void		putxD2F			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		putyD2F			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		putzD2F			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* reg1.c */
extern	void		upscaleFit		(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze);
