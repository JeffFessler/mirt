/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	void		convertOrigin	(struct	fitRec		*fit,
									 double				origin[]);
extern	void		downscaleFit	(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				l,
									 int				zSqueeze);
extern	void		freePyramids	(float				*pyrHdl1[],
									 float				*pyrHdl2[],
									 float				*mskHdl1[],
									 float				*mskHdl2[],
									 int				levels);
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
extern	void		upscaleFit		(struct	fitRec		*fit,
									 long				nxRange[],
									 long				nyRange[],
									 long				nzRange[],
									 int				fromLevel,
									 int				zSqueeze);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* pyrFilt.c															*/
/* pyrGetSz.c															*/
/* reg0.c																*/
/************************************************************************/

/* reg0.c */
extern	void		copyClip		(float				*inPtr,
									 int				nxIn,
									 int				nyIn,
									 int				nzIn,
									 float				*outPtr,
									 int				nxOut,
									 int				nyOut,
									 int				nzOut);

/* pyrFilt.c */
extern	void		pyr_filters		(float				g[],
									 short				*ng,
									 float				h[],
									 short				*nh,
									 short				idegree);

/* pyrGetSz.c */
extern	int			pyr_getsize		(int				nx,
									 int				ny,
									 int				nz,
									 int				nlevels,
									 long				nxa[],
									 long				nya[],
									 long				nza[]);
