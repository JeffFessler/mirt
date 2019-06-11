/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/

enum	bBrands {
  null,
  mirror,
  periodic,
  adhoc
};

/*--- Functions --------------------------------------------------------*/

extern	int			directBsplineFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
extern	int			directBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
extern	int			directBsplinePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
extern	int			inverseBsplineFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
extern	int			inverseBsplineMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);
extern	int			inverseBsplinePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				nx,
									 long				ny,
									 long				nz,
									 int				degree);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* convolve.c															*/
/* getPut.c																*/
/************************************************************************/

/* convolve.c */
extern	void		firConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);

/* convolve.c */
extern	void		firConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);

/* convolve.c */
extern	void		firConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);

/* getPut.c */
extern	void		getyD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		getzD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* convolve.c */
extern	void		iirConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);

/* convolve.c */
extern	void		iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);

/* convolve.c */
extern	void		iirConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);

/* getPut.c */
extern	void		putyD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);

/* getPut.c */
extern	void		putzD2D			(double				*volume,
									 long				x,
									 long				y,
									 long				z,
									 long				nx,
									 long				ny,
									 double				*row,
									 long				n);
