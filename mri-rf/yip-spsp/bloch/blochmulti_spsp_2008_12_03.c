/* Bloch simulator
This code takes in time sequences of gradients and B1 field, and 
simulates the excitation pattern over a 3 dimensional volume. 

Originally from DC Noll (5/18/03)
Adopted by CY Yip
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXNPNTS 100000
#define MAXFOVXY 128
#define MAXFOVZ 64
#define MAXNSHOTS 16
#define PI 3.14159265358979323846
#define GAM 26751.         /*Gyromagnetic ratio?*/ 

int i,j,k,f,z,tindex;
short phrf[MAXNPNTS];      /*Phase of the RF? I guess this guy goes from -1 (corr to -pi) to 1 (corr to pi) */
short magrf[MAXNPNTS];     /*Magnitude of the RF?*/
float b1xar[MAXNPNTS];     /*Real part of B1?*/
float b1yar[MAXNPNTS];     /*Imag part of B1?*/
short gx[MAXNPNTS];        /*Gradients*/
short gy[MAXNPNTS];
short gz[MAXNPNTS];
/*float freqmap[256][256];*/
float m[3], tempm[3];
short mfinal[3*MAXFOVZ*MAXFOVXY*MAXFOVXY*MAXNSHOTS];
float dm[3],dm1[3],dm2[3],dm3[3],dm4[3];
float arf;                /*Flip angle thing??*/
FILE *ofile, *fid;
char fname[50],outfname[50];
float tfinal,fa;    /*gscale is multiplied to the grads*/
int fovf,fovz,sumrf,shotn,nshots,pulsen,npulses;
float filescale;
int npnts,nzpnts;
float deltat,tfinal,dgdtmax,pwgza,pwgz,simfovz;
double pnttime, gscale, trfpw; /*changed to double type in MGH, 12/3/08, cy*/
char ifn[50], ofn[50], sfn[50];/*, mfn[50]; */
int nz, nf; /*ny;*/
float fstep,zstep;
/*int trfpw;*/
int simnpnts;
/*int simfovxy,simfovx,simfovy;*/
int simfovf;


bloch(int shind/*So I guess it is shot index?*/) {
     
     /* integration by rotational forms */
     float cphi, sphi, cpsi, spsi, ct, st;
     float Bmag, Btrans;
     float Mx0, My0, Mz0;

     int pindex,gxyindex;
     float bx, by, bz;
     float b1x, b1y, xgrad, ygrad, zgrad;
     float phi;
    
     phi=2.0*PI*((float)(shind)/nshots)*trfpw; /*What's that?*/

     for (tindex = 0; tindex< simnpnts; tindex++) {
     /*Within this shot, there are npnts where rf is defined.
       tindex goes thru those npnts*/
     pindex=tindex+shind*npnts;     
     /*the index in the original time line, not within the shot*/
     
        gxyindex=tindex;

/*     xgrad=gscale/filescale*(cos(phi)*gx[gxyindex]-sin(phi)*gy[gxyindex]);
     ygrad=gscale/filescale*(cos(phi)*gy[gxyindex]+sin(phi)*gx[gxyindex]); */
     
     /*Why multiply the gradients by gscale?*/
     xgrad=gscale*gx[pindex]/filescale;   
     ygrad=gscale*gy[pindex]/filescale;   
     zgrad=gscale*gz[pindex]/filescale;   
     
     /* Define B-fields */
     bx=(b1xar[pindex])*deltat*GAM;  /*dM = M x (B*dt*gamma)*/
     by=(b1yar[pindex])*deltat*GAM;
     /* bz=((xgrad*xstep+ygrad*ystep+zgrad*zstep)*GAM + freqmap[y][x]*2*PI)*deltat; /*freqmap is the field map. Why [y][x]?*/
     bz= ((zgrad*zstep)*GAM + 2*PI*fstep)*deltat;

     /* rotations */
     Bmag = sqrt(bx*bx+by*by+bz*bz);  /*Magnitude of B field*/
     Btrans = sqrt(bx*bx+by*by);      /*Mag in the trans plane*/ 
     ct = (Bmag > 0 )? bz/Bmag : 1;  /*t is theta, wrt +ve z axis*/
     st = sqrt(1.0 - ct*ct);         /*sine of t*/
     cphi = (Btrans > 0 )? bx/Btrans : 1; /*phi: ang. wrt to +x*/
     /* sphi = sqrt(1.0 - cphi*cphi); */  /*Why not?*/
     sphi = (Btrans > 0 )? by/Btrans : 0; 
     
	 /*M is normalized. after we we rotate the rotational axis to be aligned with B,
	 how much (angle, psi) does M move? dM/|M| = M/|M| x (gamma*B*dt). |1|*|gam*B*dt|sin(pi/2)=dM
	 |dM|=|gam*B*dt|~=the arc length of the unit circle.*/
	 cpsi = cos(Bmag);     
     spsi = sin(Bmag);

     Mx0 = m[0];
     My0 = m[1];
     Mz0 = m[2];
     m[0] =  cphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) - sphi*(-spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
     m[1] =   sphi*(ct*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0))+st*(ct*Mz0+st*(sphi*My0+cphi*Mx0))) + cphi*(-spsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + cpsi*(cphi*My0-sphi*Mx0));
     m[2] =  ct*(ct*Mz0+st*(sphi*My0+cphi*Mx0)) - st*(cpsi*(ct*(sphi*My0+cphi*Mx0)-st*Mz0) + spsi*(cphi*My0-sphi*Mx0));
     }
  	 /*But no relaxation considered!*/

     return 0;
     
 }


print_parmeters() {

  printf("files in=%s, out=%s\n",ifn,ofn);/*,mfn);*/
   printf("fa=%f\n",fa);
   printf("simfovf=%d, simfovz=%f, fovf=%d, fovz=%d\n",simfovf,simfovz,fovf,fovz);
   printf("tfinal=%f\n",tfinal);
   printf("pnttime=%1.8f\n",pnttime);
   printf("npnts=%d\n",npnts);
   printf("arf=%f\n",arf);
   printf("nshots=%d\n",nshots);
   printf("npulses=%d\n",npulses);
   printf("simnpnts=%d\n",simnpnts);
   printf("therefore sim time = %1.8f\n",simnpnts*pnttime);   

   return 0;

}
 

flip_angle() {

   float totrf, temprf;

   totrf=0.;
   for (i=0;i<npnts;i++) {
     temprf=totrf;
     totrf=temprf+magrf[i]/filescale;	 /*An integral in fact*/
   }

   /*arf is like a factor to make B1's mag right for archieving flip angle fa.
    My insight is that at the spatial location where there is no delta Bo,
    integrating B1 * gamma gives you the flip angle...*/
   arf=(fa/360)/(GAM/2/PI*pnttime*totrf);

   return 0;

}

read_pulse(int pn) {

 /*  sprintf(fname,"%s.%d.%d.%d.%d.%d.gx",ifn,npnts,fovxy,fovz,nshots,1); */
   sprintf(fname, "%s.%d.%d.%d.%d.%d.gx",ifn,npnts,fovf,fovz,nshots,pn);
   fid=fopen(fname,"r");
/*   fread(gx,2,npnts,fid); */
   fread(gx+(npnts*nshots/npulses)*(pn-1),2,npnts*nshots/npulses,fid);
   fclose(fid);
 /*  sprintf(fname,"%s.%d.%d.%d.%d.%d.gy",ifn,npnts,fovxy,fovz,nshots,1); */
   sprintf(fname,"%s.%d.%d.%d.%d.%d.gy",ifn,npnts,fovf,fovz,nshots,pn);
   fid=fopen(fname,"r");
/*   fread(gy,2,npnts,fid); */
   fread(gy+(npnts*nshots/npulses)*(pn-1),2,npnts*nshots/npulses,fid);
   fclose(fid);
   sprintf(fname,"%s.%d.%d.%d.%d.%d.gz",ifn,npnts,fovf,fovz,nshots,pn);
   fid=fopen(fname,"r");
   fread(gz+(npnts*nshots/npulses)*(pn-1),2,npnts*nshots/npulses,fid);
   fclose(fid);
   sprintf(fname,"%s.%d.%d.%d.%d.%d.mag",ifn,npnts,fovf,fovz,nshots,pn);
   fid=fopen(fname,"r");
   fread(magrf+(npnts*nshots/npulses)*(pn-1),2,npnts*nshots/npulses,fid);
   fclose(fid);
   sprintf(fname,"%s.%d.%d.%d.%d.%d.ph",ifn,npnts,fovf,fovz,nshots,pn);
   fid=fopen(fname,"r");
   fread(phrf+(npnts*nshots/npulses)*(pn-1),2,npnts*nshots/npulses,fid);
   fclose(fid);

   return 0;

}

main(int argc,char **argv) {
         
   /* Read the command line */
   if(argc == 1) {
     fprintf(stderr,"Usage: blochmulti_spsp_03Dec08 infile outfile flipangle(degrees) nf nz simfovz(cm) npnts simfovf(Hz) dfovf(Hz) dfovz(cm) nshots npulses trfpw simnpnts gmax pnttime\n");
     exit (-1);
   }
   argv++;
   strcpy(ifn,*argv++);		/*infile*/
   strcpy(ofn,*argv++);     /*outfile*/
   /*strcpy(mfn,*argv++);     mapfile*/
   sscanf(*argv++,"%f",&fa);  /*filp angle*/
   nf=atoi(*argv++);     /*number of points in the f direction*/
   /*ny=atoi(*argv++);     number of points in the y direction*/
   nz=atoi(*argv++);     /*number of points in the z direction*/
   sscanf(*argv++,"%f",&simfovz); /* bfovz*/
   npnts=atoi(*argv++);    /*number of points in a shot?*/
   simfovf=atoi(*argv++);    /*bfovf*/
   fovf=atoi(*argv++);    /*dfovf*/
   fovz=atoi(*argv++);     /*dfovz*/
   nshots=atoi(*argv++);   /*number of shots*/
   npulses=atoi(*argv++);
   trfpw=atof(*argv++);         /*?*/
   simnpnts=atoi(*argv++);   
   gscale=atof(*argv++);
   pnttime=atof(*argv++); 
 
   /* Open files and read in pulses */
   for (pulsen=0;pulsen<npulses;pulsen++)
     read_pulse(pulsen+1);

   /* fid=fopen(mfn,"r");*/
      /*   for (y=0;y<ny;y++) */
   /*     fread(freqmap[y],4,nx,fid);*/
   /*  fclose(fid);*/
 
   /* Fix parameters */
   /*gscale=4.0;*/ /*gauss/cm*/
   /*pnttime=.000004; */   /*delta t?*/            
   tfinal=pnttime*(npnts+1);    /*tfinal = the ending time?*/
   deltat=pnttime;    /*=pnttime?*/
   filescale=32767.;
   /*simfovx=1*simfovxy;*/       /*Spatial*/ 
   /*simfovy=1*simfovxy;*/
   

   /* Calculate the flip angle */
   flip_angle();
   /* Print parmeters*/ 
   print_parmeters();
     
  
   
   
   /* Begin shots loop */
   for (shotn=0;shotn<nshots;shotn++) /*This is the shot loop*/{

	   for (k = 0; k< npnts; k++) /*This is time loop in a shot*/{
		  i=k+shotn*npnts;  /*This counts time?*/
		  b1xar[i]=magrf[i]/filescale*cos(phrf[i]/filescale*PI); /*arf is the factor cal.
		  by flip_angle()*/
		  b1yar[i]=magrf[i]/filescale*sin(phrf[i]/filescale*PI);
	   }


    /* Begin xyz loops */
    for (z=0;z<nz;z++) {
      /* for (y=0;y<ny;y++) {*/
        for (f=0;f<nf;f++) {
 
         /* Initial conditions */
		 /* Notice that for each shot you initialize. Then, well, do you sum the effect
		 of all the shots? */
	 m[0]=0;
         m[1]=0; 
         m[2]=1;
     

	 /*Where in the grid?
         x_i = FOVx/nx * (x-nx/2)
	 So it sweeps from -ve to +ve.*/

         fstep=((float)(simfovf))*(f-nf/2)/nf; /*Hz*/
	 /* ystep=((float)(simfovy))*(y-ny/2)/ny;*/ /*cm*/
         zstep=(float)(simfovz)*(z-nz/2)/nz; /*cm*/
        
         /* Print out some counters */
         /*if(x==0) printf("shotn=%d y=%d (%f cm) z=%d (%f cm)\n",shotn,y,ystep,z,zstep);*/
         
         /*Apply the Bloch Equation*/ 
         bloch(shotn);
      
         /* Add to final file 
	    Somehow mfinal is defined to be 1-dimensional,
	    why not multi-D?*/
         for (k=0;k<3;k++) {
           int i;
           /*i=k+x*3+y*nx*3+z*nx*ny*3+shotn*nz*nx*ny*3;*/
           i = k+f*3+z*nf*3+shotn*nz*nf*3;   /*CY:is that right?*/
           mfinal[i]=(short)(m[k]*filescale);
         }

        } /* End f loop */
	/* }  End y loop */
    } /* End z loop */
 
   } /* End shot loop */
   /* Write to output file */
   /*sprintf(outfname,"%s.%d.%d.%d.%d.dat",ofn,nx,ny,nz,nshots);*/
   sprintf(outfname,"%s.%d.%d.%d.dat",ofn,nf,nz,nshots);
   ofile=fopen(outfname,"w");
   /*fwrite(mfinal,3*nx*ny*nz*nshots,sizeof(*mfinal),ofile);*/
   fwrite(mfinal,3*nf*nz*nshots,sizeof(*mfinal),ofile);
   fclose(ofile);

} /* End main */
