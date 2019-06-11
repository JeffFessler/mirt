/************************************************************************/
/* Public implementation follows										*/
/************************************************************************/

/*--- Defines ----------------------------------------------------------*/

#undef				ONlevel
#define				ONlevel			((int)127)

#undef				OFFlevel
#define				OFFlevel		((int)(-128))

#undef				maxPyrLevels
#define				maxPyrLevels	((int)12)

#undef				maxTapsNumber
#define				maxTapsNumber	((int)20)

#undef				hessianSize
#define				hessianSize		((int)13)		/* max(struct fitRec) */

/*--- Types ------------------------------------------------------------*/

enum	cBrands {
  gravity,
  Marquardt
};

enum	kBrands {
  blank,
  provided,
  computed
};

enum	mBrands {
  or,
  nor,
  and,
  nand,
  xor,
  nxor
};

enum	iDegree {
  zero,
  one,
  three
};

struct	fitRec {
  double			dx[3];
  double			skew[3][3];
  double			gamma;
  double			phi, theta, psi;
  double			lambda;
  double			origin[3];
};

struct	rOper {
  enum	mBrands		maskCombine;
  enum	kBrands		referenceMask, testMask;
  enum	iDegree		interpolation;
  enum	cBrands		convergence;
  int				clipping;
  int				importFit, exportFit;
  int				xTrans, yTrans, zTrans;
  int				xRot, yRot, zRot;
  int				isoScaling;
  int				xSkew, ySkew, zSkew;
  int				matchGrey;
  int				greyRendering;
  int				zapMean;
  int				zSqueeze;
};

struct	rParam {
  struct	rOper	directives;
  float				*inPtr1, *inPtr2;
  float				*inMsk1, *inMsk2;
  float				*outPtr, *mskPtr;
  char				*inFit, *outFit;
  double			sx, sy, sz;
  double			firstLambda, lambdaScale;
  double			minGain;
  double			epsilon;
  float				backgrnd;
  int				nx, ny, nz;
  int				levels, lastLevel;
};

/*--- Functions --------------------------------------------------------*/
/* None */
