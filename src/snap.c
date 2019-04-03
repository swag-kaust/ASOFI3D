/*------------------------------------------------------------------------
 *   Write 3D snapshot for current timestep  to disk                                   
 *
 *  ----------------------------------------------------------------------*/

#include "data_structures.h"
#include "fd.h"
#include "globvar.h"


void snap(FILE *fp, int nt, int nsnap, int format, int type, 
        Velocity *v, Tensor3d *s,
        float ***u, float ***pi,
        int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2, 
        int ny2, int nz2) {

	/* 
	different data formats of output available:
	format=1  :  SU (IEEE)
	format=2  :  ASCII
	format=3  :  BINARY (IEEE)
	
	different types:
	type=1 : values in vx, vy, and vz
	type=2 : -(sxx+syy+szz) (pressure field)
	type=3 : divergence of vx, vy and vz (energy of compressional waves)
	         and curl of vx, vy and vz (energy of shear waves)
	type=4 : both particle velocities (type=1) and energy (type=3)
	*/


	
	char xfile[STRING_SIZE], yfile[STRING_SIZE], zfile[STRING_SIZE];
	char rotfile[STRING_SIZE], ext[8], wm[2];
	char  divfile[STRING_SIZE], pfile[STRING_SIZE];
	FILE *fpx1, *fpy1, *fpz1, *fpx2, *fpy2, *fpp;
	int i,j,k;
	float a=0.0, amp, dh24x, dh24y, dh24z, vyx, vxy, vxx, vyy, vzx, vyz, vxz, vzy, vzz;


	extern float DX, DY, DZ, DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int MYID, POS[4], SNAP_PLANE, LOG;

        float ***vx = v->x;
        float ***vy = v->y;
        float ***vz = v->z;

        float ***sxx = s->xx;
        float ***syy = s->yy;
        float ***szz = s->zz;

	switch(format){
	case 1: 
		sprintf(ext,".su");
		break;
	case 2: 
		sprintf(ext,".asc");
		break;
	case 3: 
		sprintf(ext,".bin");
		break;
	}


	sprintf(xfile,"%s%s.vx.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(yfile,"%s%s.vy.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(zfile,"%s%s.vz.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(divfile,"%s%s.div.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(rotfile,"%s%s.curl.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(pfile,"%s%s.p.%i.%i.%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);

        if (LOG){
	fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);}


	if (nsnap==1) 
		sprintf(wm,"w");
	else 
		sprintf(wm,"a");

	switch(type){
	case 1 :
		fprintf(fp,"\t%s\n", xfile);
		fprintf(fp,"\t%s\n", yfile);
		fprintf(fp,"\t%s\n\n", zfile);
		fpx1=fopen(xfile,wm);
		fpy1=fopen(yfile,wm);
		fpz1=fopen(zfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
				
			
					writedsk(fpx1,vx[j][i][k],format);
					writedsk(fpy1,vy[j][i][k],format);
					writedsk(fpz1,vz[j][i][k],format);
					
					
				}
		fclose(fpx1);
		fclose(fpy1);
		fclose(fpz1);
		break;
	case 2 :
		fprintf(fp,"\t%s\n\n", pfile);
		fpp=fopen(pfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					amp=-sxx[j][i][k]-syy[j][i][k]-szz[j][i][k];
					
				
					writedsk(fpp,amp,format);

				}
		fclose(fpp);
		break;
	case 4 :
		fprintf(fp,"\t%s\n", xfile);
		fprintf(fp,"\t%s\n", yfile);
		fprintf(fp,"\t%s\n\n", zfile);
		fprintf(fp,"\t%s\n\n", pfile);
		fpx1=fopen(xfile,wm);
		fpy1=fopen(yfile,wm);
		fpz1=fopen(zfile,wm);
		fpp=fopen(pfile,wm);
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					amp=-sxx[j][i][k]-syy[j][i][k]-szz[j][i][k];

					writedsk(fpx1,vx[j][i][k],format);
					writedsk(fpy1,vy[j][i][k],format);
					writedsk(fpz1,vz[j][i][k],format);
					writedsk(fpp,amp,format);
					
				}
		fclose(fpx1);
		fclose(fpy1);
		fclose(fpz1);
		fclose(fpp);
	case 3 :
		/* output of the curl of the velocity field according to Dougherty and
		                  Stephen (PAGEOPH, 1988) */
		fprintf(fp,"\t%s\n", divfile);
		fprintf(fp,"\t%s\n\n", rotfile);
		fpx2=fopen(divfile,wm);
		fpy2=fopen(rotfile,wm);
		
		dh24x=1.0/DX;
		dh24y=1.0/DY;
		dh24z=1.0/DZ;
		
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					vxy=(vx[j+1][i][k]-vx[j][i][k])*(dh24y);
					vxz=(vx[j][i][k+1]-vx[j][i][k])*(dh24z);
					vyx=(vy[j][i+1][k]-vy[j][i][k])*(dh24x);
					vyz=(vy[j][i][k+1]-vy[j][i][k])*(dh24z);
					vzx=(vz[j][i+1][k]-vz[j][i][k])*(dh24x);
					vzy=(vz[j+1][i][k]-vz[j][i][k])*(dh24y);
					
					/*amp= absolute value of curl(v)), without sqrt!!!*/
					amp=((vzy-vyz)*(vzy-vyz)+(vxz-vzx)*(vxz-vzx)+(vyx-vxy)*(vyx-vxy)); 
					
					/*note that "y" denotes the vertical coordinate*/
					switch(SNAP_PLANE){
					case 1 : /* energy without sign */
						/* sqrt(Es) with Es = u*amp*amp (second amp removed due to missing sqrt in amp*/
						a=sqrt((u[j][i][k])*amp);
						break;
					case 2 : /* energy with sign true for x-y-plane */
						/*sign(rot(v)x * sqrt(Es) with Es = u*amp*amp (second amp removed due to missing sqrt in amp*/
						a=fsign((vxy-vyx))*sqrt((u[j][i][k])*amp);
						break;
					case 3 : /* energy with sign true for x-z-plane */
						/*sign(rot(v)r * sqrt(Es) with Es = u*amp*amp (second amp removed due to missing sqrt in amp*/
						a=fsign((vxz-vzx))*sqrt((u[j][i][k])*amp);
						break;
					case 4 :/* energy with sign true for y-z-plane */
						/*sign(rot(v)t * sqrt(Es) with Es = u*amp*amp (second amp removed due to missing sqrt in amp*/
						a=fsign((vzy-vyz))*sqrt((u[j][i][k])*amp);
						break;
					case 5 : /*custom force*/ /*not yet working properly*/
						/*sign(rot(v)t * sqrt(Es) with Es = u*amp*amp (second amp removed due to missing sqrt in amp*/
						a=fsign((vxz-vzx))*sqrt((u[j][i][k])*amp);
						break;
					case 6 : /*SH wave*/
						amp=(vxz-vzx)*(vxz-vzx);
						a=fsign((vxz-vzx))*sqrt((u[j][i][k])*amp);
						break;
					case 7 : /*SV wave*/
						amp=(vzy-vyz)*(vzy-vyz)+(vyx-vxy)*(vyx-vxy);
						a=sqrt((u[j][i][k])*amp);
						break;
					
					}
					
					writedsk(fpy2,a,format);
					   

				}



		/* output of the divergence of the velocity field according to Dougherty and
		                  Stephen (PAGEOPH, 1988) */
		for (k=nz1;k<=nz2;k+=idz)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=ny1;j<=ny2;j+=idy){
					vxx=(vx[j][i][k]-vx[j][i-1][k])*(dh24x);
					vyy=(vy[j][i][k]-vy[j-1][i][k])*(dh24y);
					vzz=(vz[j][i][k]-vz[j][i][k-1])*(dh24z);
					
					/*amp= div(v))*/
					amp=(vxx+vyy+vzz);
					
					switch(SNAP_PLANE){
					case 1 : /* energy without sign */
						/* Ep with Ep=pi*amp*amp */
						a=sqrt((pi[j][i][k])*amp*amp);
						break;
					case 2 : /* single force in x */
						/*sign of div(v) * Ep with Ep=pi*amp*amp */
						a=fsign(amp)*sqrt((pi[j][i][k])*amp*amp);
						break;
					case 3 : /* single force in y */
						/*sign of div(v) * Ep with Ep=pi*amp*amp */
						a=fsign(amp)*sqrt((pi[j][i][k])*amp*amp);
						break;
					case 4 : /* single force in z */
						/*sign of div(v) * Ep with Ep=pi*amp*amp */
						a=fsign(amp)*sqrt((pi[j][i][k])*amp*amp);
						break;
					}
					
					writedsk(fpx2,a,format);
				}

		fclose(fpx2);
		fclose(fpy2);
		break;
	}
}


