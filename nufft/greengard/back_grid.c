#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	V0	prhs[0]
#define	MSP	prhs[1]
#define M1      prhs[2]
#define M2      prhs[3]
#define MR      prhs[4]
#define E2X     prhs[5]
#define E2Y     prhs[6]
#define E3      prhs[7]
#define OMLEN   prhs[8]

/* Output Arguments */

#define	FK	plhs[0]

#define PI 3.14159265

static void dogridding(double *v0r,
		       double *v0i,
		       double *Fkr,
		       double *Fki,
		       int msp,
		       double *m1,
		       double *m2,
		       int mr,
		       int omlen,
		       double *e2x,
		       double *e2y,
		       double *e3)
{
  
  register int j,l2,l1;
  register int rowind,colind;
  register double vyr,vyi,kern1,kern2;

  register double *m1_p,*m2_p;
  register double *e2x_p,*e2y_p;
  register double *e3_p1,*e3_p2;
  register double *Fkr_p,*Fki_p;
  register double *v0r_p,*v0i_p;
  register double inve2x,e2x_posl1,e2x_negl1;
  register double inve2y,e2y_posl2,e2y_negl2;

  m1_p = m1;
  m2_p = m2;
  e2x_p = e2x;
  e2y_p = e2y;
  v0r_p = v0r;
  v0i_p = v0i;

  for(j=0;j<omlen;j++)
    {
      e3_p1 = e3;

      /* case l2 = 0 */
      colind = *m2_p+mr;
      colind = colind % mr;
      kern2 = *e3_p1++;
      vyr = *v0r_p * kern2;
      vyi = *v0i_p * kern2;
      Fkr_p = &Fkr[colind*mr];
      Fki_p = &Fki[colind*mr];
	  
      e3_p2 = e3;
      
      /* case l1 = 0 */
      rowind = *m1_p+mr;
      rowind = rowind % mr;
      kern1 = *e3_p2++;
      *(Fkr_p + rowind) += vyr * kern1;
      *(Fki_p + rowind) += vyi * kern1;
      
      /* case 0<abs(l1)<msp */
      inve2x = 1.0 / *e2x_p;
      e2x_posl1 = *e2x_p;
      e2x_negl1 = inve2x;
      
      for(l1=1;l1<msp;l1++)
	{
	  /* positive l1 */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_posl1*(*e3_p2);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	  e2x_posl1 *= *e2x_p;
	  
	  /* negative l1 */
	  rowind = (*m1_p-l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_negl1*(*e3_p2++);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	  e2x_negl1 *= inve2x;
	}
      
      /* case l1 = msp */
      /* positive l1 only */
      rowind = (*m1_p+l1)+mr;
      rowind = rowind % mr;
      kern1 = e2x_posl1*(*e3_p2);
      *(Fkr_p + rowind) += vyr * kern1;
      *(Fki_p + rowind) += vyi * kern1;
      
      /* case 0<abs(l2)<msp */
      inve2y = 1.0 / *e2y_p;
      e2y_posl2 = *e2y_p;
      e2y_negl2 = inve2y;

      for(l2=1;l2<msp;l2++)
	{
	  /* positive l2 */
	  colind = (*m2_p+l2)+mr;
	  colind = colind % mr;
	  kern2 = e2y_posl2*(*e3_p1);
	  vyr = *v0r_p * kern2;
	  vyi = *v0i_p * kern2;
	  Fkr_p = &Fkr[colind*mr];
	  Fki_p = &Fki[colind*mr];
	  
	  e3_p2 = e3;
	  
	  /* case l1 = 0 */
	  rowind = *m1_p+mr;
	  rowind = rowind % mr;
	  kern1 = *e3_p2++;
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	      
	  /* case 0<abs(l1)<msp */
	  inve2x = 1.0 / *e2x_p;
	  e2x_posl1 = *e2x_p;
	  e2x_negl1 = inve2x;
      
	  for(l1=1;l1<msp;l1++)
	    {
	      /* positive l1 */
	      rowind = (*m1_p+l1)+mr;
	      rowind = rowind % mr;
	      kern1 = e2x_posl1*(*e3_p2);
	      *(Fkr_p + rowind) += vyr * kern1;
	      *(Fki_p + rowind) += vyi * kern1;
	      e2x_posl1 *= *e2x_p;

	      /* negative l1 */
	      rowind = (*m1_p-l1)+mr;
	      rowind = rowind % mr;
	      kern1 = e2x_negl1*(*e3_p2++);
	      *(Fkr_p + rowind) += vyr * kern1;
	      *(Fki_p + rowind) += vyi * kern1;
	      e2x_negl1 *= inve2x;
	    }

	  /* case l1 = msp */
	  /* positive l1 only */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_posl1*(*e3_p2);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;

	  e2y_posl2 *= *e2y_p;

	  /* negative l2 */
	  colind = (*m2_p-l2)+mr;
	  colind = colind % mr;
	  kern2 = e2y_negl2*(*e3_p1++);
	  vyr = *v0r_p * kern2;
	  vyi = *v0i_p * kern2;
	  Fkr_p = &Fkr[colind*mr];
	  Fki_p = &Fki[colind*mr];
	  
	  e3_p2 = e3;
	  
	  /* case l1 = 0 */
	  rowind = *m1_p+mr;
	  rowind = rowind % mr;
	  kern1 = *e3_p2++;
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	      
	  /* case 0<abs(l1)<msp */
	  inve2x = 1.0 / *e2x_p;
	  e2x_posl1 = *e2x_p;
	  e2x_negl1 = inve2x;
      
	  for(l1=1;l1<msp;l1++)
	    {
	      /* positive l1 */
	      rowind = (*m1_p+l1)+mr;
	      rowind = rowind % mr;
	      kern1 = e2x_posl1*(*e3_p2);
	      *(Fkr_p + rowind) += vyr * kern1;
	      *(Fki_p + rowind) += vyi * kern1;
	      e2x_posl1 *= *e2x_p;

	      /* negative l1 */
	      rowind = (*m1_p-l1)+mr;
	      rowind = rowind % mr;
	      kern1 = e2x_negl1*(*e3_p2++);
	      *(Fkr_p + rowind) += vyr * kern1;
	      *(Fki_p + rowind) += vyi * kern1;
	      e2x_negl1 *= inve2x;
	    }

	  /* case l1 = msp */
	  /* positive l1 only */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_posl1*(*e3_p2);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;

	  e2y_negl2 *= inve2y;
	}

      /* case l2 = 0 */
      /* positive l2 only */
      colind = (*m2_p+l2)+mr;
      colind = colind % mr;
      kern2 = e2y_posl2*(*e3_p1);
      vyr = *v0r_p * kern2;
      vyi = *v0i_p * kern2;
      Fkr_p = &Fkr[colind*mr];
      Fki_p = &Fki[colind*mr];
      
      e3_p2 = e3;
      
      /* case l1 = 0 */
      rowind = *m1_p+mr;
      rowind = rowind % mr;
      kern1 = *e3_p2++;
      *(Fkr_p + rowind) += vyr * kern1;
      *(Fki_p + rowind) += vyi * kern1;
      
      /* case 0<abs(l1)<msp */
      inve2x = 1.0 / *e2x_p;
      e2x_posl1 = *e2x_p;
      e2x_negl1 = inve2x;
      
      for(l1=1;l1<msp;l1++)
	{
	  /* positive l1 */
	  rowind = (*m1_p+l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_posl1*(*e3_p2);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	  e2x_posl1 *= *e2x_p;
	  
	  /* negative l1 */
	  rowind = (*m1_p-l1)+mr;
	  rowind = rowind % mr;
	  kern1 = e2x_negl1*(*e3_p2++);
	  *(Fkr_p + rowind) += vyr * kern1;
	  *(Fki_p + rowind) += vyi * kern1;
	  e2x_negl1 *= inve2x;
	}
      
      /* case l1 = msp */
      /* positive l1 only */
      rowind = (*m1_p+l1)+mr;
      rowind = rowind % mr;
      kern1 = e2x_posl1*(*e3_p2);
      *(Fkr_p + rowind) += vyr * kern1;
      *(Fki_p + rowind) += vyi * kern1;
      
      m1_p++;
      m2_p++;
      e2x_p++;
      e2y_p++;
      v0r_p++;
      v0i_p++;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  double *v0r,*v0i,*Fkr,*Fki,*e2x,*e2y,*e3; 
  double *m1,*m2;

  /* Check for proper number of arguments */
  
  if (nrhs != 9)
    mexErrMsgTxt("Nine input arguments required."); 
  if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments."); 
  if (!mxIsComplex(V0)) 
    mexErrMsgTxt("V0 must be complex.");
  
  /* Assign pointers to the various parameters */ 
  v0r = mxGetPr(V0);
  v0i = mxGetPi(V0);
  e2x = mxGetPr(E2X);
  e2y = mxGetPr(E2Y);
  e3 = mxGetPr(E3);
  m1 = mxGetPr(M1);
  m2 = mxGetPr(M2);
  
  int msp = mxGetScalar(MSP);
  int mr = mxGetScalar(MR);
  int omlen = mxGetScalar(OMLEN);
  
  /* Create a matrix for the return argument */ 
  FK = mxCreateDoubleMatrix(mr, mr, mxCOMPLEX); 
  Fkr = mxGetPr(FK);
  Fki = mxGetPi(FK);
    
  /* Do the actual computations in a subroutine */
  dogridding(v0r,v0i,Fkr,Fki,msp,m1,m2,mr,omlen,e2x,e2y,e3);
  
  return;
  
}
