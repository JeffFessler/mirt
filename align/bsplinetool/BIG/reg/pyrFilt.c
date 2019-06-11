#include <stdlib.h>

#include <stdio.h>
#define DEBUG 0
/*************************************************************
 * Title : pyrFilt.c
 *
 * PURPOSE :  Initializes down- and up-sampling filter arrays for
 * least squares splines of degree 0 to 3.  (little l_2 norm)
 * g : reduce filter
 * h : expand filter
 *
 * NOTE : filter arrays should be defined as 
 * 			float g[20],h[20];
 *
 *  Michael Unser / BEIP	MAY-92
 *
 *************************************************************/
void pyr_filters(float g[],short *ng,float h[],short *nh,short idegree);
void pyr_filters(float g[],short *ng,float h[],short *nh,short idegree)
/* float g[],h[];	filter arrays
   short *ng,*nh;	number of taps
	short idegree;	degree of the spline */	
{
	switch (idegree) {
	case 0 :
		*ng=1; *nh=1;
		break;

	case 1 :
		g[0]=0.707107; g[1]=0.292893; g[2]=-0.12132; g[3]=-0.0502525;
		g[4]=0.0208153; g[5]=0.00862197; g[6]=-0.00357134;
		g[7]=-0.0014793; g[8]=0.000612745;
		*ng=9;
		h[0]=1.; h[1]=0.5;
		*nh=2;
		break;
	case 2 :
		g[0]=0.617317; g[1]=0.310754; g[2]=-0.0949641; g[3]=-0.0858654;
		g[4]=0.0529153; g[5]=0.0362437; g[6]=-0.0240408;
		g[7]=-0.0160987; g[8]=0.0107498; g[9]=0.00718418;
		g[10]=-0.00480004; g[11]=-0.00320734; g[12]=0.00214306;
		g[13]=0.00143195; g[14]=-0.0009568; g[15]=-0.000639312;
		*ng=16;
		h[0]=1.; h[1]=0.585786; h[2]=0; h[3]=-0.100505; h[4]=0;
		h[5]=0.0172439; h[6]=0; h[7]=-0.00295859; h[8]=0;
		h[9]=0.000507614;
		*nh=10;
		break;
	case 3 :
		g[0]=0.596797; g[1]=0.313287; g[2]=-0.0827691; g[3]=-0.0921993;
		g[4]=0.0540288; g[5]=0.0436996; g[6]=-0.0302508;
		g[7]=-0.0225552; g[8]=0.0162251; g[9]=0.0118738;
		g[10]=-0.00861788; g[11]=-0.00627964; g[12]=0.00456713;
		g[13]=0.00332464; g[14]=-0.00241916; g[15]=-0.00176059;
		g[16]=0.00128128; g[17]=0.000932349; g[18]=-0.000678643;
		g[19]=-0.000493682;
		*ng=20;
		h[0]=1.; h[1]=0.600481; h[2]=0; h[3]=-0.127405; h[4]=0;
		h[5]=0.034138; h[6]=0; h[7]=-0.00914725; h[8]=0;
		h[9]=0.002451; h[10]=0; h[11]=-0.000656743;
		*nh=12;
		break;
	default :
		*ng=1; *nh=1;
		break;
	}
}
