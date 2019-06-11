/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

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

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* reg0.c																*/
/* svdcmp.c																*/
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

/* reg0.c */
extern	void		eulerRotation	(struct	fitRec		*fit);

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

/* svdcmp.c */
extern	int			svbksb			(double				*u,
									 double				w[],
									 double				*v,
									 long				m,
									 long				n,
									 double				b[],
									 double				x[]);

/* svdcmp.c */
extern	int			svdcmp			(double				*a,
									 long				m,
									 long				n,
									 double				w[],
									 double				*v);
