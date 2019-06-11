/* Bloch simulation: capable with parallel simu*/
/* Originally download from Brian Hargreaves's website */
/* Hao Sun made the following change: (Jul 22, 2012)
 *(1)match field map with position
 *(2)correct the minus sign error in the b1imag in the blochsimfz sub function
 *(3)add parallel simulation capability
 */
#include "mex.h"
#include <stdio.h>
#include <math.h>

#define GAMMA   26751.3
#define TWOPI	6.283185307179586

#define DEBUG



void multmatvec(double *mat, double *vec, double *matvec)

/* Multiply 3x3 matrix by 3x1 vector. */

{
    *matvec++ = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
    *matvec++ = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
    *matvec++ = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
}



void addvecs(double *vec1, double *vec2, double *vecsum)

/* Add two 3x1 Vectors */

{
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
}




void adjmat(double *mat, double *adj)

/* ======== Adjoint of a 3x3 matrix ========= */

{
    *adj++ = (mat[4]*mat[8]-mat[7]*mat[5]);
    *adj++ =-(mat[1]*mat[8]-mat[7]*mat[2]);
    *adj++ = (mat[1]*mat[5]-mat[4]*mat[2]);
    *adj++ =-(mat[3]*mat[8]-mat[6]*mat[5]);
    *adj++ = (mat[0]*mat[8]-mat[6]*mat[2]);
    *adj++ =-(mat[0]*mat[5]-mat[3]*mat[2]);
    *adj++ = (mat[3]*mat[7]-mat[6]*mat[4]);
    *adj++ =-(mat[0]*mat[7]-mat[6]*mat[1]);
    *adj++ = (mat[0]*mat[4]-mat[3]*mat[1]);
}


void zeromat(double *mat)

/* ====== Set a 3x3 matrix to all zeros	======= */

{
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
}


void eyemat(double *mat)

/* ======== Return 3x3 Identity Matrix  ========= */

{
    zeromat(mat);
    mat[0]=1;
    mat[4]=1;
    mat[8]=1;
    
}

double detmat(double *mat)

/* ======== Determinant of a 3x3 matrix ======== */

{
    double det;
    
    det = mat[0]*mat[4]*mat[8];
    det+= mat[3]*mat[7]*mat[2];
    det+= mat[6]*mat[1]*mat[5];
    det-= mat[0]*mat[7]*mat[5];
    det-= mat[3]*mat[1]*mat[8];
    det-= mat[6]*mat[4]*mat[2];
    
    return det;
}


void scalemat(double *mat, double scalar)

/* ======== multiply a matrix by a scalar ========= */

{
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
}


void invmat(double *mat, double *imat)

/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
    double det;
    int count;
    
    det = detmat(mat);	/* Determinant */
    adjmat(mat, imat);	/* Adjoint */
    
    for (count=0; count<9; count++)
        *imat = *imat++ / det;
}


void addmats(double *mat1, double *mat2, double *matsum)

/* ====== Add two 3x3 matrices.	====== */

{
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
}


double multmats(double *mat1, double *mat2, double *matproduct)

/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */

{
    *matproduct++ = mat1[0]*mat2[0] + mat1[3]*mat2[1] + mat1[6]*mat2[2];
    *matproduct++ = mat1[1]*mat2[0] + mat1[4]*mat2[1] + mat1[7]*mat2[2];
    *matproduct++ = mat1[2]*mat2[0] + mat1[5]*mat2[1] + mat1[8]*mat2[2];
    *matproduct++ = mat1[0]*mat2[3] + mat1[3]*mat2[4] + mat1[6]*mat2[5];
    *matproduct++ = mat1[1]*mat2[3] + mat1[4]*mat2[4] + mat1[7]*mat2[5];
    *matproduct++ = mat1[2]*mat2[3] + mat1[5]*mat2[4] + mat1[8]*mat2[5];
    *matproduct++ = mat1[0]*mat2[6] + mat1[3]*mat2[7] + mat1[6]*mat2[8];
    *matproduct++ = mat1[1]*mat2[6] + mat1[4]*mat2[7] + mat1[7]*mat2[8];
    *matproduct++ = mat1[2]*mat2[6] + mat1[5]*mat2[7] + mat1[8]*mat2[8];
}


double calcrotmat(double nx, double ny, double nz, double *rmat)

/* Find the rotation matrix that rotates |n| radians about
 * the vector given by nx,ny,nz				*/
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;
    
    /* printf("calcrotmat(%6.3f,%6.3f,%6.3f) -> phi = %6.3f\n",nx,ny,nz,phi); */
    
    phi = sqrt(nx*nx+ny*ny+nz*nz);
    
    if (phi == 0.0)
    {
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
    }
    
    
    
    else
    {
        /* First define Cayley-Klein parameters 	*/
        hp = phi/2;
        cp = cos(hp);
        sp = sin(hp)/phi;	/* /phi because n is unit length in defs. */
        ar = cp;
        ai = -nz*sp;
        br = ny*sp;
        bi = -nx*sp;
        
        /* Make auxiliary variables to speed this up	*/
        
        arar = ar*ar;
        aiai = ai*ai;
        arai2 = 2*ar*ai;
        brbr = br*br;
        bibi = bi*bi;
        brbi2 = 2*br*bi;
        arbi2 = 2*ar*bi;
        aibr2 = 2*ai*br;
        arbr2 = 2*ar*br;
        aibi2 = 2*ai*bi;
        
        
        /* Make rotation matrix.	*/
        
        *rmat++ = arar-aiai-brbr+bibi;
        *rmat++ = -arai2-brbi2;
        *rmat++ = -arbr2+aibi2;
        *rmat++ =  arai2-brbi2;
        *rmat++ = arar-aiai+brbr-bibi;
        *rmat++ = -aibr2-arbi2;
        *rmat++ =  arbr2+aibi2;
        *rmat++ =  arbi2-aibr2;
        *rmat++ = arar+aiai-brbr-bibi;
    }
}



void zerovec(double *vec)

/*	Set a 3x1 vector to all zeros	*/

{
    *vec++=0;
    *vec++=0;
    *vec++=0;
}


int times2intervals( double *endtimes, double *intervals, long n)
/* ------------------------------------------------------------
 * Function takes the given endtimes of intervals, and
 * returns the interval lengths in an array, assuming that
 * the first interval starts at 0.
 *
 * If the intervals are all greater than 0, then this
 * returns 1, otherwise it returns 0.
 * ------------------------------------------------------------ */

{
    int allpos;
    int count;
    double lasttime;
    
    allpos=1;
    lasttime = 0.0;
    
    for (count = 0; count < n; count++)
    {
        intervals[count] = endtimes[count]-lasttime;
        lasttime = endtimes[count];
        if (intervals[count] <= 0)
            allpos =0;
    }
    
    return (allpos);
}






int blochsim(double *b1real, double *b1imag,
        double *xgrad, double *ygrad, double *zgrad, double *tsteps,
        int ntime, double *e1, double *e2, double df,
        double dx, double dy, double dz,
        double *mx, double *my, double *mz, int mode)
        
        /* Go through time for one df and one dx,dy,dz.		*/
        
{
    int count;
    int tcount;
    double gammadx;
    double gammady;
    double gammadz;
    double rotmat[9];
    double amat[9], bvec[3];	/* A and B propagation matrix and vector */
    double arot[9], brot[3];	/* A and B after rotation step. */
    double decmat[9];		/* Decay matrix for each time step. */
    double decvec[3];		/* Recovery vector for each time step. */
    double rotx,roty,rotz;		/* Rotation axis coordinates. */
    double mstart[3];
    double mfinish[3];
    double imat[9], mvec[3];
    double mcurr0[3];		/* Current magnetization before rotation. */
    double mcurr1[3];		/* Current magnetization before decay. */
    
    eyemat(amat); 		/* A is the identity matrix.	*/
    eyemat(imat); 		/* I is the identity matrix.	*/
    
    zerovec(bvec);
    zerovec(decvec);
    zeromat(decmat);
    
    gammadx = dx*GAMMA;	/* Convert to Hz/cm */
    gammady = dy*GAMMA;	/* Convert to Hz/cm */
    gammadz = dz*GAMMA;	/* Convert to Hz/cm */
    
    
    mcurr0[0] = *mx;		/* Set starting x magnetization */
    mcurr0[1] = *my;		/* Set starting y magnetization */
    mcurr0[2] = *mz;		/* Set starting z magnetization */
    
    
    for (tcount = 0; tcount < ntime; tcount++)
    {
        /*	Rotation 	*/
        /*Hao:  the rotation in this code is right handed. So all the rotation matrix is negated */
        rotz = -(*xgrad++ * gammadx + *ygrad++ * gammady + *zgrad++ * gammadz +
                df*TWOPI ) * *tsteps;
        rotx = (- *b1real++ * GAMMA * *tsteps);
        roty = (- *b1imag++ * GAMMA * *tsteps++); /*Hao: changed to minus*/
        calcrotmat(rotx, roty, rotz, rotmat);
        
        if (mode == 1 || mode == 101)
        {
            multmats(rotmat,amat,arot);
            multmatvec(rotmat,bvec,brot);
        }
        else
            multmatvec(rotmat,mcurr0,mcurr1);
        
        
        /* 	Decay	*/
        
        decvec[2]= 1- *e1;
        decmat[0]= *e2;
        decmat[4]= *e2++;
        decmat[8]= *e1++;
        
        if (mode == 1 || mode == 101)
        {
            multmats(decmat,arot,amat);
            multmatvec(decmat,brot,bvec);
            addvecs(bvec,decvec,bvec);
        }
        else
        {
            multmatvec(decmat,mcurr1,mcurr0);
            addvecs(mcurr0,decvec,mcurr0);
        }
        
        /*
         * printf("rotmat = [%6.3f  %6.3f  %6.3f ] \n",rotmat[0],rotmat[3],
         * rotmat[6]);
         * printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[1],rotmat[4],
         * rotmat[7]);
         * printf("         [%6.3f  %6.3f  %6.3f ] \n",rotmat[2],rotmat[5],
         * rotmat[8]);
         * printf("A = [%6.3f  %6.3f  %6.3f ] \n",amat[0],amat[3],amat[6]);
         * printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[1],amat[4],amat[7]);
         * printf("    [%6.3f  %6.3f  %6.3f ] \n",amat[2],amat[5],amat[8]);
         * printf(" B = <%6.3f,%6.3f,%6.3f> \n",bvec[0],bvec[1],bvec[2]);
         * printf("<mx,my,mz> = <%6.3f,%6.3f,%6.3f> \n",
         * amat[6] + bvec[0], amat[7] + bvec[1], amat[8] + bvec[2]);
         *
         * printf("\n");
         */
        
        if (mode == 10 || mode == 110)		/* Sample output at times.  */
            /* Only do this if transient! */
        {
            *mx = mcurr0[0];
            *my = mcurr0[1];
            *mz = mcurr0[2];
            
            mx++;
            my++;
            mz++;
        }
    }
    
    
    
    /* If only recording the endpoint, either store the last
     * point, or calculate the steady-state endpoint. */
    
    if (mode == 0 || mode == 100)		/* Indicates start at given m, or m0. */
    {
        *mx = mcurr0[0];
        *my = mcurr0[1];
        *mz = mcurr0[2];
    }
    
    else if (mode == 1 || mode == 101)	/* Indicates to find steady-state magnetization */
    {
        scalemat(amat,-1.0);		/* Negate A matrix 	*/
        addmats(amat,imat,amat);	/* Now amat = (I-A)		*/
        invmat(amat,imat);		/* Reuse imat as inv(I-A) 	*/
        multmatvec(imat,bvec,mvec);	/* Now M = inv(I-A)*B		*/
        *mx = mvec[0];
        *my = mvec[1];
        *mz = mvec[2];
    }
    
    
}



int blochsimfz(double *b1real, double *b1imag, double *xgrad, double *ygrad, double *zgrad,
        double *tsteps,
        int ntime, double t1, double t2, double *dfreq, int nfreq,
        double *dxpos, double *dypos, double *dzpos, int npos,
        double *mx, double *my, double *mz, int mode, int ncoils, double *dsensr, double *dsensi,
        int printflag, int uniSense)
        
        
        
{
    int count;
    int poscount;
    int fcount;
    int totpoints;
    int totcount = 0;
    int coilcount;
    int timecount;
    
    int ntout;
    
    double *e1;
    double *e2;
    double *e1ptr;
    double *e2ptr;
    double *tstepsptr;
    double *dxptr, *dyptr, *dzptr;
    double *b1comr, *b1comi, *b1comrTmp, *b1comiTmp;
    double *b1realTmp, *b1imagTmp;
    double *dsensrTmp, *dsensiTmp;
    
    if (mode & 2)
        ntout = ntime;
    else
        ntout = 1;
    
    /* First calculate the E1 and E2 values at each time step. */
    
    e1 = (double *) malloc(ntime * sizeof(double));
    e2 = (double *) malloc(ntime * sizeof(double));
    
    e1ptr = e1;
    e2ptr = e2;
    tstepsptr = tsteps;
    /*printf("mode = %d.\n", mode);*/
    /*printf("unit sense = %d.\n", uniSense);*/
    for (count=0; count < ntime; count++)
    {
        *e1ptr++ = exp(- *tstepsptr / t1);
        *e2ptr++ = exp(- *tstepsptr++ / t2);
    }
    
    totpoints = npos;
    
    dxptr = dxpos;
    dyptr = dypos;
    dzptr = dzpos;
    
    b1comr = (double *) malloc(ntime * sizeof(double));
    b1comi = (double *) malloc(ntime * sizeof(double));
    
    for (poscount=0; poscount < npos; poscount++)
        
    {
        
        /* used to weighted combine b1 waves of different coil*/
        
        b1comrTmp = b1comr;
        b1comiTmp = b1comi;
        for(timecount = 0; timecount < ntime; timecount++)
        {
            *b1comrTmp = 0;
            *b1comiTmp = 0;
            b1realTmp = b1real+timecount;
            b1imagTmp = b1imag+timecount;
            dsensrTmp = dsensr;
            dsensiTmp = dsensi;
            for(coilcount = 0; coilcount < ncoils; coilcount++)
            {
                *b1comrTmp = *b1comrTmp + *b1realTmp * (*dsensrTmp) - *b1imagTmp * (*dsensiTmp);
                *b1comiTmp = *b1comiTmp + *b1imagTmp * (*dsensrTmp) + *b1realTmp * (*dsensiTmp);
                b1realTmp += ntime;
                b1imagTmp += ntime;
                if(uniSense != 1)
                {
                    dsensrTmp++;
                    dsensiTmp++;
                }
            }
            
            /* printf("the b1com is %6.3f, %6.3f \n",*b1comrTmp,*b1comiTmp); */
            b1comrTmp++;
            b1comiTmp++;
        }
        
        if(uniSense != 1)
        {
            dsensr += ncoils;
            dsensi += ncoils;
        }
        
        if (mode == 11 || mode == 111)	/* Steady state AND record all time points. */
        /* original is mode == 3, seems not working correctly, change to mode == 11 */
        {	/* First go through and find steady state, then
         * repeat as if transient starting at steady st.*/
            
            blochsim(b1comr, b1comi, xgrad, ygrad, zgrad, tsteps, ntime,
                    e1, e2, *dfreq, *dxptr, *dyptr,
                    *dzptr, mx, my, mz, 1);
            
            blochsim(b1comr, b1comi, xgrad, ygrad, zgrad, tsteps, ntime,
                    e1, e2, *dfreq++, *dxptr++, *dyptr++,
                    *dzptr++, mx, my, mz, 10);
            /* Hao: change the *dfreq and *dfreq++ order */
        }
        else
        {
            blochsim(b1comr, b1comi, xgrad, ygrad, zgrad, tsteps, ntime,
                    e1, e2, *dfreq++, *dxptr++, *dyptr++,
                    *dzptr++, mx, my, mz, mode);
        }
        
        mx += ntout;
        my += ntout;
        mz += ntout;
        
        
        totcount++;
        if ((totpoints > 40000) && ( ((10*totcount)/totpoints)> (10*(totcount-1)/totpoints) ) && printflag == 1)
            printf("%d%% Complete.\n",(100*totcount/totpoints));
    }
    
    free(e1);
    free(e2);
    free(b1comr);
    free(b1comi);
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

/* bloch(b1,gradxyz,dt,t1,t2,df,dx,dy,dz,mode,mx,my,mz) */
{
    double *b1r;	/* Real-part of B1 field.	*/
    double *b1i;	/* Imag-part of B1 field.	*/
    double *gx;	/* X-axis gradient. 		*/
    double *gy;	/* Y-axis gradient. 		*/
    double *gz;	/* Z-axis gradient. 		*/
    double *tp;	/* Time steps (s)		*/
    double *ti;	/* Time intervals (s) 		*/
    double t1;	/* T1 time constant (s)	 	*/
    double t2;	/* T2 time constant (s)		*/
    double *df;	/* Off-resonance Frequencies (Hz)	*/
    double *dx;	/* X Positions (cm)			*/
    double *dy;	/* Y Positions (cm)			*/
    double *dz;	/* Z Positions (cm)			*/
    int md;		/* Mode - 0=from M0, 1=steady-state	*/
    double tstep;	/* Time step, if single parameter */
    double *mxin;	/* Input points given by matlab*/
    double *myin;
    double *mzin;
    
    double *mxout;	/* used to pass the data of mxin to mx at the beginning, in the later computation routine. mx my mz will be updated*/
    double *myout;
    double *mzout;
    
    double *mx;	/* Output Arrays */
    double *my;
    double *mz;
    
    int gyaflag=0;	/* 1 if gy was allocated. */
    int gzaflag=0;	/* 1 if gy was allocated. */
    int dyaflag=0;	/* 1 if dy was allocated. */
    int dzaflag=0;	/* 1 if dy was allocated. */
    
    int ntime;	/* Number of time points. 	 */
    int ntout;	/* Number of time poitns at output. */
    int outsize[3];	/* Output matrix sizes		*/
    
    int ngrad;	/* Number of gradient dimensions */
    int nf;
    int npos;	/* Number of positions.  Calculated from nposN and nposM, depends on them. */
    int nposM;	/* Height of passed position matrix. */
    int nposN;	/* Width of passed position matrix. */
    int nfnpos;	/* Number of frequencies * number of positions. */
    int ntnpos;     /*  Number of output times * number of positions. */
    int ntnfnpos;	/* Number of output times *frequencies*number of positions. */
    int count;
    
    int ncoils;      /*number of coils */
    int nsens;       /*number of sensitivity point */
    double *dsensr, *dsensi; /* pointer to sensitivity maps */
    int nb1;      /* number of b1 points */
    int uniSense; /* set to 1 if sensitivies are all 1 */
    
    int printflag; /*set to 1 if print out information*/
    
    
    /* ===== Mode, defaults to 0 (print info, simulate single endpoint, transient). ==== */
    
    if (nrhs > 7)
        md = (int)(*mxGetPr(prhs[7]));
    else
        md = 0;
    
    if ((md & 4)==0)
    {
        printflag = 1;
    }
    else
        printflag = 0;
    
    if (printflag == 1)
    {
        if ((md & 1)==0)
            printf("Simulation from Initial Condition.\n");
        else
            printf("Simulation of Steady-State.\n");
        
        
        if ((md & 2)==0)
            printf("Simulation to Endpoint. \n");
        else
            printf("Simulation over Time.\n");
    }
    
    
    if(printflag == 1){
        printf("---------------------------------------\n");
        printf("3D-position + frequency Bloch Simulator\n");
        printf("---------------------------------------\n\n");
    }
    
    
    
    /* ====================== RF (B1) =========================
     * :  If complex, split up.  If real, allocate an imaginary part. ==== */
    nb1 = mxGetM(prhs[0]) * mxGetN(prhs[0]);	/* Number of Time, RF, and Grad points */
    
    
    if (mxIsComplex(prhs[0]))
    {
        b1r = mxGetPr(prhs[0]);
        b1i = mxGetPi(prhs[0]);
    }
    else
    {
        b1r = mxGetPr(prhs[0]);
        b1i = (double *)malloc(nb1 * sizeof(double));
        for (count=0; count < nb1; count++)
            b1i[count]=0.0;
    }
    if(printflag == 1){
        printf("%d B1 points.\n",nb1);
    }
    
    
    
    /* ======================= Gradients ========================= */
    
    ngrad = mxGetM(prhs[1]) * mxGetN(prhs[1]);	/* Number of Time,Grad points */
    ntime = ngrad/3;   /* the input have to contain gx,gy,gz */
    gx = mxGetPr(prhs[1]);				/* X-gradient is first N points. */
    gy = gx + ntime;
    gz = gy + ntime;
    
    ncoils = nb1/ntime;
    /* printf(" %d, %d, %d \n",nb1,ntime,ngrad); */
    if(printflag == 1){
        printf("%d coils. \n",ncoils);
    }
    
    if(printflag == 1){
        if (gx == NULL)
            printf("ERROR:  gx is not allocated. \n");
        if (gy == NULL)
            printf("ERROR:  gy is not allocated. \n");
        if (gz == NULL)
            printf("ERROR:  gz is not allocated. \n");
    }
    
    
    /* === Time points ===== */
    
    /*	THREE Cases:
     * 1) Single value given -> this is the interval length for all.
     * 2) List of intervals given.
     * 3) Monotonically INCREASING list of end times given.
     *
     * For all cases, the goal is for tp to have the intervals.
     */
    
    ti = NULL;
    tp = mxGetPr(prhs[2]);
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)	/* === Case 1 === */
    {
        tp = (double *)malloc(ntime * sizeof(double));
        tstep = *(mxGetPr(prhs[2]));
        for (count =0; count < ntime; count++)
            tp[count]=tstep;
    }
    else if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != ntime)
    {
        if(printflag == 1){
            printf("Time-point length differs from B1 length.\n");
        }
    }
    else
    {
        tp = mxGetPr(prhs[2]);
        ti = (double *)malloc(ntime * sizeof(double));
        if (( times2intervals( tp, ti, ntime )))
        {
            if(printflag == 1){
                printf("Times are monotonically increasing. \n");
            }
            tp = ti;
        }
    }
    
    
    /* === Relaxation Times ===== */
    
    t1 = *mxGetPr(prhs[3]);
    t2 = *mxGetPr(prhs[4]);
    
    
    /* === Frequency Points ===== */
    
    df = mxGetPr(prhs[5]);
    nf = mxGetM(prhs[5]) * mxGetN(prhs[5]);
    
    if(printflag == 1){
        printf("%d Frequency points.\n",nf);
    }
    
    
    /* === Position Points ===== */
    
    nposM = mxGetM(prhs[6]);
    nposN = mxGetN(prhs[6]);
    
    if(printflag == 1){
        printf("Position vector is: %d x %d. \n",nposM,nposN);
    }
    
    if (nposN==3)			/* Assume 3 position dimensions given */
    {
        npos = nposM;
        if(printflag == 1){
            printf("Assuming %d 3-Dimensional Positions\n",npos);
        }
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = dy + npos;
    }
    
    else if (nposN==2)		/* Assume only 2 position dimensions given */
    {
        npos = nposM;
        if(printflag == 1){
            printf("Assuming %d 2-Dimensional Positions\n",npos);
        }
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = (double *)malloc(npos * sizeof(double));
        dzaflag=1;
        for (count=0; count < npos; count++)
            dz[count]=0.0;
    }
    
    else				/* Either 1xN, Nx1 or something random.  In all these
     * cases we assume that 1 position is given, because it
     * is too much work to try to figure out anything else! */
    {
        npos = nposM * nposN;
        if(printflag == 1){
            printf("Assuming %d 1-Dimensional Positions\n",npos);
        }
        dx = mxGetPr(prhs[6]);
        dy = (double *)malloc(npos * sizeof(double));
        dz = (double *)malloc(npos * sizeof(double));
        dyaflag=1;
        dzaflag=1;
        for (count=0; count < npos; count++)
        {
            dy[count]=0.0;
            dz[count]=0.0;
        }
        if(printflag == 1){
            if ((nposM !=1) && (nposN!=1))
            {
                printf("Position vector should be 1xN, Nx1, Nx2 or Nx3. \n");
                printf(" -> Assuming 1 position dimension is given. \n");
            }
        }
    }
    
    if (dx == NULL)
        printf("ERROR:  dx is not allocated. \n");
    if (dy == NULL)
        printf("ERROR:  dy is not allocated. \n");
    if (dz == NULL)
        printf("ERROR:  dz is not allocated. \n");
    
    nfnpos = nf*npos;	/* Just used to speed things up below. 	*/
    if (nf != npos)
        printf("Warning: dimension mismatch between freq. and location. \n");
    
    /*============ load sens ===============*/
    if (nrhs > 8)
    {
        nsens = mxGetM(prhs[8]) * mxGetN(prhs[8]);
        if (nsens > 1)
        {
            uniSense = 0; 
            if(printflag == 1){
                printf("Sens vector is: %d x %d. \n",mxGetM(prhs[8]),mxGetN(prhs[8]));
            }
            if((nsens/npos) != ncoils)
                printf("Warning: sens dimenstion seems not correct. \n");
            
            if (mxIsComplex(prhs[8]))
            {
                dsensr = mxGetPr(prhs[8]);
                dsensi = mxGetPi(prhs[8]);
                /* printf("SENS IS COMPLEX. \n");   */
            }
            else
            {
                dsensr = mxGetPr(prhs[8]);
                dsensi = (double *)malloc(nsens * sizeof(double));
                for (count = 0; count < nsens; count++)
                    dsensi[count] = 0.0;
            }
        }
        else
        {
            uniSense = 1;
            nsens = npos * ncoils;
            dsensr = (double *)malloc(sizeof(double));
            dsensi = (double *)malloc(sizeof(double));
            *dsensr = 1.0;
            *dsensi = 0.0;
        }
        
    }
    else
    {
        printf("Warning: sens matrix is not provided. Assume uniform sense. \n");
        uniSense = 1;
        nsens = npos * ncoils;
        dsensr = (double *)malloc(sizeof(double));
        dsensi = (double *)malloc(sizeof(double));
        *dsensr = 1.0; 
        *dsensi = 0.0; 
    }
    
    
    
    
    
    /* ===== Allocate Output Magnetization vectors arrays.	*/
    if (md & 2)
        ntout = ntime;		/* Include time points.	*/
    else
        ntout = 1;
    
    if(printflag == 1){
        printf("Mode = %d, %d Output Time Points \n",md,ntout);
    }
    ntnpos = ntout*npos;
    ntnfnpos = ntout*nfnpos;
    
    plhs[0] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* Mx, output. */
    plhs[1] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* My, output. */
    plhs[2] = mxCreateDoubleMatrix(ntnpos,1,mxREAL);	/* Mz, output. */
    
    mx = mxGetPr(plhs[0]);
    my = mxGetPr(plhs[1]);
    mz = mxGetPr(plhs[2]);
    
    mxout = mx;
    myout = my;
    mzout = mz;
    
    /* ===== If Initial Magnetization is given... */
    
    if ( (nrhs > 11) &&
            (mxGetM(prhs[9]) * mxGetN(prhs[9]) == npos) &&
            (mxGetM(prhs[10]) * mxGetN(prhs[10]) == npos) &&
            (mxGetM(prhs[11]) * mxGetN(prhs[11]) == npos)  )
        
        /* Set output magnetization to that passed.
         * If multiple time points, then just the
         * first is set.				*/
        
        
    {
        if(printflag == 1){
            printf("Using Specified Initial Magnetization.\n");
        }
        
        mxin = mxGetPr(prhs[9]);
        myin = mxGetPr(prhs[10]);
        mzin = mxGetPr(prhs[11]);
        for (count =0; count < npos; count++)
        {
            *mxout = *mxin++;  /*this opeartion first assign hte data pointed by mxin to data pointed by mxout, then increase the pointer(address value) of mxin*/
            *myout = *myin++;
            *mzout = *mzin++;
            mxout += ntout;
            myout += ntout;
            mzout += ntout;
        }
    }
    else
    {
        if(printflag == 1){
            if (nrhs > 11) 	 /* Magnetization given, but wrong size! */
            {
                printf("Initial magnetization passed, but not Npositions x Nfreq. \n");
            }
            printf(" --> Using [0; 0; 1] for initial magnetization. \n");
        }
        for (count =0; count < npos; count++)
        {
            *mxout = 0;	/* Set magnetization to Equilibrium */
            *myout = 0;
            *mzout = 1;
            mxout += ntout;
            myout += ntout;
            mzout += ntout;
        }
    }
    
    
    /* ======= Do The Simulation! ====== */
    
    if(printflag == 1){
        printf("Calling blochsimfz() function in Mex file.\n");
    }
    
    blochsimfz(b1r,b1i,gx,gy,gz,tp,ntime,t1,t2,df,nf,dx,dy,dz,npos,mx,my,mz,md,ncoils,dsensr,dsensi,printflag,uniSense);
    
    
    /* ======= Reshape Output Matrices ====== */
    
    if ((ntout > 1) && (nf > 1) && (npos > 1))
    {
        outsize[0]=ntout;
        /*outsize[1]=npos;*/
        outsize[1]=nf;
        mxSetDimensions(plhs[0],outsize,2);  /* Set to 3D array. */
        mxSetDimensions(plhs[1],outsize,2);  /* Set to 3D array. */
        mxSetDimensions(plhs[2],outsize,2);  /* Set to 3D array. */
    }
    else			/* Basically "squeeze" the matrix. */
    {
        if (ntout > 1)
        {
            outsize[0]=ntout;
            /*outsize[1]=npos*nf;*/
            outsize[1]=npos;
        }
        else  /*this is the case for SPINS simulation*/
        {
            outsize[0]=npos;
            /*outsize[1]=nf;*/
        }
        mxSetDimensions(plhs[0],outsize,1);  /* Set to 2D array. */
        mxSetDimensions(plhs[1],outsize,1);  /* Set to 2D array. */
        mxSetDimensions(plhs[2],outsize,1);  /* Set to 2D array. */
    }
    
    
    /* ====== Free up allocated memory, if necessary. ===== */
    
    if (!mxIsComplex(prhs[0]))
        free(b1i);	/* We had to allocate this before. */
    
    if (nrhs > 8)
    {
        if (nsens > 1)
        {
            if (!mxIsComplex(prhs[8]))
                free(dsensi);	/* We had to allocate this before. */
        }
        else
        {
            /*
            free(dsensr);
            free(dsensi);
             */
        }
    }
    else
    {
        /*
        free(dsensr);
        free(dsensi);
         */
    }
    
    
    
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)
        free(tp);	/* We had to allocate this. */
    
    if (ti != NULL)
        free(ti);
    
    if (dyaflag==1)
        free(dy);
    if (dzaflag==1)
        free(dz);
    if (gyaflag==1)
        free(gy);
    if (gzaflag==1)
        free(gz);
    
}






