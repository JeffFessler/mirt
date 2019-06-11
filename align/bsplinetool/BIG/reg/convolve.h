/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/
/* None */

/*--- Functions --------------------------------------------------------*/

extern	void		firConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);
extern	void		firConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);
extern	void		firConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 long				shift,
									 double				*kernel,
									 long				nk);
extern	void		iirConvolveFinite
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);
extern	void		iirConvolveMirror
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);
extern	void		iirConvolvePeriodic
									(double				*inPtr,
									 double				*outPtr,
									 long				n,
									 double				gain,
									 double				*realPoles,
									 long				np);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* None																	*/
/************************************************************************/
