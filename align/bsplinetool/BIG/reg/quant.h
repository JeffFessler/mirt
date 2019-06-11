/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/
/* None */

/*--- Types ------------------------------------------------------------*/

enum	qBrands {
  linear,
  equalize,
  LloydMax,
  none,
  from
};

struct	qParam {
  union {
	short			*s;
	float			*f;
  }					inPtr;
  union {
	short			*s;
	float			*f;
  }					outPtr;
  float				*tstMsk, *refMsk;
  char				*fileName;
  enum	qBrands		mode, labels;
  double			entropy;
  double			mse, snr;
  double			epsilon;
  float				threshold;
  float				foregrnd, backgrnd;
  float				oMin, oMax;
  int				floatInput, floatOutput;
  int				nx, ny, nz;
  int				parzenDegree;
  short				iMin, iMax;
  short				first, step;
  unsigned			bins;
};

struct	histList {
  struct	histList
					*next;
  double			occurences;
  float				value;
};

/*--- Functions --------------------------------------------------------*/

extern	int			histEntropy		(struct	qParam		*entDataPtr);
extern	int			histogram		(struct qParam		*histDataPtr);
extern	int			histQuant		(struct	qParam		*quantDataPtr);
extern	int			parzenEntropy	(struct qParam		*entDataPtr);
extern	void		quantizeFromBins(float				*inPtr,
									 float				*binPtr,
									 short				*outPtr,
									 int				nx,
									 int				ny,
									 int				nz,
									 unsigned			bins);
extern	int			sliceQuant		(struct	qParam		*sliceDataPtr);
extern	int			thrQuant		(struct	qParam		*thrDataPtr);

/************************************************************************/
/* Externals															*/
/*----------------------------------------------------------------------*/
/* BsplnWgt.c															*/
/************************************************************************/

/* BsplnWgt.c */
extern	double		BsplnWght		(int				degree,
									 long				i,
									 double				x);
