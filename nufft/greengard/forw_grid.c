#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	FMINTAU	prhs[0]
#define	MSP	prhs[1]
#define M1      prhs[2]
#define M2      prhs[3]
#define MR      prhs[4]
#define E2X     prhs[5]
#define E2Y     prhs[6]
#define E3      prhs[7]
#define OMLEN   prhs[8]

/* Output Arguments */

#define	FJ	plhs[0]

#define PI 3.14159265

static void dogridding(double *fmintaur,
		       double *fmintaui,
		       double *fjr,
		       double *fji,
		       int msp,
		       double *m1,
		       double *m2,
		       int mr,
		       int omlen,
		       double *e2x,
		       double *e2y,
		       double *e3)
{
  
  register int j,l1,l2;
  register int rowind,colind;
  register double kern1,kern2;
  register double fjr_tmp,fji_tmp;
  register double *fmintaur_p,*fmintaui_p;
  register double *m1_p,*m2_p;
  register double *e2x_p,*e2y_p;
  register double *fjr_p,*fji_p;
  register double *e3_p1,*e3_p2;
  register double inve2x,e2x_posl1,e2x_negl1;
  register double inve2y,e2y_posl2,e2y_negl2;

  m1_p = m1;
  m2_p = m2;
  e2x_p = e2x;
  e2y_p = e2y;
  fjr_p = fjr;
  fji_p = fji;

  for(j=0;j<omlen;j++)
    {
      fjr_tmp = 0.0;fji_tmp = 0.0;
      e3_p1 = e3;

      /* case l2 = 0 */
      colind = *m2_p+mr;
      colind = colind % mr;
      kern2 = *e3_p1++;
      fmintaur_p = &fmintaur[colind*mr];
      fmintaui_p = &fmintaui[colind*mr];

      e3_p2 = e3;
      
      /* case l1 = 0 */
      rowind = *m1_p+mr;
      rowind = rowind % mr;
      kern1 = kern2 * (*e3_p2++);
      fjr_tmp += *(fmintaur_p + rowind) * kern1;
      fji_tmp += *(fmintaui_p + rowind) * kern1;
      
      /* case 0<abs(l1)<msp */
      inve2x  = 1.0 / *e2x_p;
      e2x_posl1 = *e2x_p;
      e2x_negl1 = inve2x;
      
      for(l1=1;l1<msp;l1++)
	{
	  /* positive l1 */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_posl1*(*e3_p2);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;
	  e2x_posl1 *= *e2x_p;
	  
	  /* negative l1 */
	  rowind = (*m1_p-l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_negl1*(*e3_p2++);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;
	  e2x_negl1 *= inve2x;
	}
      
      /* case l1 = msp */
      /* positive l1 only */
      rowind = (*m1_p+l1)+mr;
      rowind = rowind % mr;
      kern1 = kern2*e2x_posl1*(*e3_p2);
      fjr_tmp += *(fmintaur_p + rowind) * kern1;
      fji_tmp += *(fmintaui_p + rowind) * kern1;    
      
      /* case 0<abs(l2)<msp */
      inve2y = 1.0 / *e2y_p;
      e2y_posl2 = *e2y_p;
      e2y_negl2 = inve2y;

      for(l2=1;l2<msp;l2++)
	{
	  /* positive l2 */
	  colind = (*m2_p + l2)+mr; /* make sure result of mod is not negative */
	  colind = colind % mr;
	  kern2 = e2y_posl2*(*e3_p1);
	  fmintaur_p = &fmintaur[colind*mr];
	  fmintaui_p = &fmintaui[colind*mr];
	  e3_p2 = e3;

	  /* case l1 = 0 */
	  rowind = *m1_p+mr;
	  rowind = rowind % mr;
	  kern1 = kern2 * (*e3_p2++);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;

	  /* case 0<abs(l1)<msp */
	  inve2x  = 1.0 / *e2x_p;
	  e2x_posl1 = *e2x_p;
	  e2x_negl1 = inve2x;

	  for(l1=1;l1<msp;l1++)
	    {
	      /* positive l1 */
	      rowind = (*m1_p+l1)+mr;
	      rowind = rowind % mr;
	      kern1 = kern2*e2x_posl1*(*e3_p2);
	      fjr_tmp += *(fmintaur_p + rowind) * kern1;
	      fji_tmp += *(fmintaui_p + rowind) * kern1;
	      e2x_posl1 *= *e2x_p;

	      /* negative l1 */
	      rowind = (*m1_p-l1)+mr;
	      rowind = rowind % mr;
	      kern1 = kern2*e2x_negl1*(*e3_p2++);
	      fjr_tmp += *(fmintaur_p + rowind) * kern1;
	      fji_tmp += *(fmintaui_p + rowind) * kern1;
	      e2x_negl1 *= inve2x;
	    }

	  /* case l1 = msp */
	  /* positive l1 only */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_posl1*(*e3_p2);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;    
	  
	  e2y_posl2 *= *e2y_p;

	  /* negative l2 */
	  colind = (*m2_p - l2)+mr; /* make sure result of mod is not negative */
	  colind = colind % mr;
	  kern2 = e2y_negl2*(*e3_p1++);
	  fmintaur_p = &fmintaur[colind*mr];
	  fmintaui_p = &fmintaui[colind*mr];
	  e3_p2 = e3;

	  /* case l1 = 0 */
	  rowind = *m1_p+mr;
	  rowind = rowind % mr;
	  kern1 = kern2 * (*e3_p2++);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;

	  /* case 0<abs(l1)<msp */
	  inve2x  = 1.0 / *e2x_p;
	  e2x_posl1 = *e2x_p;
	  e2x_negl1 = inve2x;

	  for(l1=1;l1<msp;l1++)
	    {
	      /* positive l1 */
	      rowind = (*m1_p+l1)+mr;
	      rowind = rowind % mr;
	      kern1 = kern2*e2x_posl1*(*e3_p2);
	      fjr_tmp += *(fmintaur_p + rowind) * kern1;
	      fji_tmp += *(fmintaui_p + rowind) * kern1;
	      e2x_posl1 *= *e2x_p;

	      /* negative l1 */
	      rowind = (*m1_p-l1)+mr;
	      rowind = rowind % mr;
	      kern1 = kern2*e2x_negl1*(*e3_p2++);
	      fjr_tmp += *(fmintaur_p + rowind) * kern1;
	      fji_tmp += *(fmintaui_p + rowind) * kern1;
	      e2x_negl1 *= inve2x;
	    }

	  /* case l1 = msp */
	  /* positive l1 only */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_posl1*(*e3_p2);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;    

	  e2y_negl2 *= inve2y;

	}

      /* case l2 = msp */
      /* positive l2 only */
      /* positive l2 */
      colind = (*m2_p + l2)+mr; /* make sure result of mod is not negative */
      colind = colind % mr;
      kern2 = e2y_posl2*(*e3_p1);
      fmintaur_p = &fmintaur[colind*mr];
      fmintaui_p = &fmintaui[colind*mr];
      e3_p2 = e3;
      
      /* case l1 = 0 */
      rowind = *m1_p+mr;
      rowind = rowind % mr;
      kern1 = kern2 * (*e3_p2++);
      fjr_tmp += *(fmintaur_p + rowind) * kern1;
      fji_tmp += *(fmintaui_p + rowind) * kern1;
      
      /* case 0<abs(l1)<msp */
      inve2x  = 1.0 / *e2x_p;
      e2x_posl1 = *e2x_p;
      e2x_negl1 = inve2x;
      
      for(l1=1;l1<msp;l1++)
	{
	  /* positive l1 */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_posl1*(*e3_p2);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;
	  e2x_posl1 *= *e2x_p;
	  
	  /* negative l1 */
	  rowind = (*m1_p-l1)+mr;
	  rowind = rowind % mr;
	  kern1 = kern2*e2x_negl1*(*e3_p2++);
	  fjr_tmp += *(fmintaur_p + rowind) * kern1;
	  fji_tmp += *(fmintaui_p + rowind) * kern1;
	  e2x_negl1 *= inve2x;
	}
      
      /* case l1 = msp */
      /* positive l1 only */
      rowind = (*m1_p+l1)+mr;
      rowind = rowind % mr;
      kern1 = kern2*e2x_posl1*(*e3_p2);
      fjr_tmp += *(fmintaur_p + rowind) * kern1;
      fji_tmp += *(fmintaui_p + rowind) * kern1;  

      *fjr++ = fjr_tmp;
      *fji++ = fji_tmp;
      m1_p++;
      m2_p++;
      e2x_p++;
      e2y_p++;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  double *fmintaur,*fmintaui,*fjr,*fji,*e2x,*e2y,*e3; 
  double *m1,*m2;

  /* Check for proper number of arguments */
  
  if (nrhs != 9)
    mexErrMsgTxt("Nine input arguments required."); 
  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments."); 
  if (!mxIsComplex(FMINTAU)) 
    mexErrMsgTxt("fmintau must be complex.");
  
  /* Assign pointers to the various parameters */ 
  fmintaur = mxGetPr(FMINTAU);
  fmintaui = mxGetPi(FMINTAU);
  e2x = mxGetPr(E2X);
  e2y = mxGetPr(E2Y);
  e3 = mxGetPr(E3);
  m1 = mxGetPr(M1);
  m2 = mxGetPr(M2);
  
  int msp = mxGetScalar(MSP);
  int mr = mxGetScalar(MR);
  int omlen = mxGetScalar(OMLEN);
  
  /* Create a matrix for the return argument */ 
  FJ = mxCreateDoubleMatrix(omlen, 1, mxCOMPLEX); 
  fjr = mxGetPr(FJ);
  fji = mxGetPi(FJ);
    
  /* Do the actual computations in a subroutine */
  dogridding(fmintaur,fmintaui,fjr,fji,msp,m1,m2,mr,omlen,e2x,e2y,e3);
  
  return;
  
}
