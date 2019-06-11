/*
mvmult.c
Matlab mex for multiplying matrix times vector
where the matrix is the system matrix for MRI
defined by parameters
we,kx,ky,t,s
Brad Sutton and Jeff Fessler - 2/15/00

Matlab Usage: mvmult(kx,ky,t,we,x,xval,yval,trans,fov,n)
Where A is system matrix given by kx,ky,t,we,fov,n
and the function returns A*x for trans = 0.
The parameters are defined as follows:
	kx is a column vector of all kx values for each sample
	ky is the corresponding ky values at each sample in a column vector
	t is a time column vector
	we is an xnx1 inhomogeneity field map, with entries corresponding to
	the locations specified in the xval and yval vectors
	x is the vector to be multiplied, generally complex
	xval xvalues of length np
	yval yvalues of length np
	trans is a 1 if A'*x is desired
	fov is the field of view in cm
	n is the matrix size, ie. 64, 128
*/



# include "mex.h"
# include "math.h"
# include "matrix.h"
# include <stdlib.h>

# define PI 3.14159265
# define GAM 2*PI*4257.5
/* that is in Hz/G for protons */

// Sinc Function
double sinc(double x)
{
	if (fabs(x) < 0.0001)
		return 1.0;
	else
		return sin(PI*x)/(PI*x);
}


void sysmult(double yr[],double yi[],double kx[],double ky[],
		double t[],double we[],double xr[],double xi[],
		double xval[],double yval[], int tm,int np,
		double trans[],
		double fov[], double NN[])
{
 double dx, dy, sumr, sumi, prefix, tpi, argument, cosarg, sinarg, rex, imx,
    scx, scy, fctr;
 int i, j, p, q;

mexPrintf("tm = %i, np = %i \n",tm, np);

 fctr = 1.0;  /*(fov[0]*fov[0])/(NN[0]*NN[0]); */

  tpi = 2*PI;
  if (trans[0] == 0)
    {
      for (i=0;i<tm;i++)
	{
	sumr = 0;
	sumi = 0;
	scx = sinc(kx[i]/NN[0]);
	scy = sinc(ky[i]/NN[0]);
	for (j=0;j<np;j++)
	{
		argument = (tpi*(kx[i]*xval[j]+ky[i]*yval[j]))+(we[j]*t[i]);
		cosarg = cos(argument);
		sinarg = sin(argument);
		rex = xr[j];
		imx = xi[j];
		sumr += ((cosarg*rex)+(sinarg*imx));
		sumi += ((-sinarg*rex)+(cosarg*imx));
	}
	yr[i] = fctr*scx*scy*sumr;   /* *prefix; */
	yi[i] = fctr*scx*scy*sumi;   /* *prefix; */
	}
  } else {    /* We want the conjugate transpose of A */
       for (j=0;j<np;j++)
	{
		sumr = 0;
		sumi = 0;
		for (i=0;i<tm;i++)
		{
		scx = sinc(kx[i]/NN[0]);
		scy = sinc(ky[i]/NN[0]);
		argument = (tpi*(kx[i]*xval[j]+ky[i]*yval[j]))+(we[j]*t[i]);
		cosarg = cos(argument);
		sinarg = sin(argument);
		rex = xr[i];
		imx = xi[i];
		sumr += scx*scy*((cosarg*rex)-(sinarg*imx)); // *prefix
		sumi += scx*scy*((sinarg*rex)+(cosarg*imx)); // *prefix
		}
	 yr[j] = fctr*sumr;
	 yi[j] = fctr*sumi;
	 }
    }
}



/*
 *Gateway routine - Interface with Matlab
 */

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  double *kx, *ky, *t, *we, *xr, *xi, *xval, *yval, *trans;
  double *yr, *yi, *fov, *NN;
  int kxm, kxn, kym, kyn, tn, tm, xm, xn, np, xvn, yp, yvn, wem, wen,
	nfov, mfov, nN, mN;


  if (nlhs > 1){
    mexErrMsgTxt("One output argument required.");
  } else {if (nrhs != 10) {
    mexErrMsgTxt("Ten input arguments required");
  }}
  /*  Get the sizes of all matrices input, to determine if proper */
  /*	    input is given					  */

  kxm = mxGetM(prhs[0]);
  kxn = mxGetN(prhs[0]);
  kym = mxGetM(prhs[1]);
  kyn = mxGetN(prhs[1]);
  tm = mxGetM(prhs[2]);
  tn = mxGetN(prhs[2]);
  wem = mxGetM(prhs[3]);
  wen = mxGetN(prhs[3]);
  np = mxGetM(prhs[5]);
  xvn = mxGetN(prhs[5]);
  yp = mxGetM(prhs[6]);
  yvn = mxGetN(prhs[6]);
  xm = mxGetM(prhs[4]);
  xn = mxGetN(prhs[4]);
  mfov = mxGetM(prhs[8]);
  nfov = mxGetN(prhs[8]);
  nN = mxGetN(prhs[9]);
  mN = mxGetM(prhs[9]);

  /* Also check size of we */

  if ((mfov != 1)||(nfov !=1)){
    mexErrMsgTxt("FOV must be a scalar");
  }

  if ((mN != 1)||(nN !=1)){
    mexErrMsgTxt("Matrix size must be a scalar");
  }

  /*  If the input matrices are not of correct size, print error message */
  if ((kxm != kym)||(kxm != tm)||(kxn != 1)||(kyn != 1)||(tn != 1)){
    mexErrMsgTxt("Vectors kx,ky, and t must be n_timex1");
  }

  /* Check sizes of xval and yval and we to make sure they agree */
  if ((xvn != 1)||(yvn !=1)||(np != yp)){
    mexErrMsgTxt("Xval and Yval vectors must be npx1");
  }

	// Size of output matrix depends on whether or not the transpose
	// is desired.
 trans = mxGetPr(prhs[7]);

  // Check size of input matrix
  if ((xn != 1)){
    mexErrMsgTxt("Input vector must be a column vector");
  }

  if ((xm != np)&(trans[0] == 0)){
    mexErrMsgTxt("x Must be np x 1 when trans = 0");
  }

  if ((xm != tm)&(trans[0] == 1)){
    mexErrMsgTxt("x must be tm*snum when trans = 1");
  }

  if ((wem != np)||(wen != 1)){
    mexErrMsgTxt("Field Map is not appropriate size, must be npx1");
  }


  /* Must create output matrix for values of yr and yi  */
 if (trans[0] == 1){
    plhs[0] = mxCreateDoubleMatrix(np,1,mxCOMPLEX);
 }else{
    plhs[0] = mxCreateDoubleMatrix(tm,1,mxCOMPLEX);
 }

  // Load in all variables, must load real and imaginary parts
  // separately

  kx = mxGetPr(prhs[0]);
  ky = mxGetPr(prhs[1]);
  t = mxGetPr(prhs[2]);
  we = mxGetPr(prhs[3]);
  if (!(xr = mxGetPr(prhs[4]))){
  mexErrMsgTxt("Real Part of input vector not detected");
  }
  if (!(xi = mxGetPi(prhs[4]))){
  xi=mxCalloc(xm, sizeof(double));
  }

  xval = mxGetPr(prhs[5]);
  yval = mxGetPr(prhs[6]);
  trans = mxGetPr(prhs[7]);
  fov = mxGetPr(prhs[8]);
  NN = mxGetPr(prhs[9]);

  yr = mxGetPr(plhs[0]);
  yi = mxGetPi(plhs[0]);

  sysmult(yr, yi, kx, ky, t, we, xr, xi, xval, yval, tm, np, trans, fov, NN);
}
