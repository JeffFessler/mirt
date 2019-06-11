/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	void		computeSkewHessian
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
extern	int			spatialGradient	(float				*gradHdl[],
									 float				*inPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 enum	iDegree		interpolation);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* BsplnTrf.c															*/
/* convolve.c															*/
/* getPut.c																*/
/* reg0.c																*/
/* reg3.c																*/
/************************************************************************/

/* reg0.c */
extern	double		absDeterminant	(double				skew[][3]);

/* reg0.c */
extern	void		combineMasks	(float				*inMsk,
									 float				*outMsk,
									 enum	mBrands		maskCombine,
									 int				nx,
									 int				ny,
									 int				nz);

/* reg0.c */
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

/* reg0.c */
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

/* BsplnTrf.c */
extern	int			directBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);

/* reg0.c */
extern	void		eulerRotation	(struct	fitRec		*fit);

/* convolve.c */
extern	void		firConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);

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
extern	void		getyD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*column,
									 long				n);

/* getPut.c */
extern	void		getzD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillar,
									 long				n);

/* reg0.c */
extern	int			invertFit		(struct	fitRec		*dirFit,
									 struct	fitRec		*invFit);

/* reg0.c */
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

/* reg3.c */
extern	int			marquardt		(struct	rOper		*directives,
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
									 int				nz);

/* getPut.c */
extern	void		putxD2F			(float				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*column,
									 long				n);

/* getPut.c */
extern	void		putyD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*column,
									 long				n);

/* getPut.c */
extern	void		putzD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*pillar,
									 long				n);
