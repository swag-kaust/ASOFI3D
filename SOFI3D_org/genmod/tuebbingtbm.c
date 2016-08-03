/*
 *   homogeneous half space 
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, DX, DY, DZ, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	/*extern char  MFILE[STRING_SIZE];*/
	char  outfile[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts, ts, tp, sumu, sumpi, d,d2,x,y,z, dsteel;
	int i, j, k, l, ii, jj, kk;
	/*char modfile[STRING_SIZE];*/
	
	/* background */
	const float vp1=5700.0, vs1=3400.0, rho1=2200.0, Qp1=500.0, Qs1=500.0;
	
	/* background reflector */
        const float vp2=4500.0, vs2=2800.0, rho2=1500.0, Qp2=500.0, Qs2=500.0;	
	const float reflect_dist=80.00;	
 
	/* tunnel */
	const float vpt=0.0, vst=1.0e-4, rhot=1.25, Qpt=1000.0, Qst=1000.0, length_tunnel=60.00,

	axis_y=(float)NYG*DY/2.0, axis_z=(float)NZG*DZ/2.0, radius=5.0;
		    
	/*excevation damage zone*/
	const float vpedz=1800.0*sqrt(3.0), vsedz=1800.0, rhoedz=1500.0, 
	Qpedz=500.0, Qsedz=500.0,

	radiusedz=10.0;
	
	/* tuebbing */
	
	const float vpconcrete=4000.0, vsconcrete=vpconcrete/sqrt(3), rhoconcrete=2300.0, 
	Qpconcrete=500.0, Qsconcrete=500.0, 
	vpsteel=5860.0, vssteel=3313.0, rhosteel=7700.0, 
	Qpsteel=500.00, Qssteel=500.0,
	vpbackfill=3000.0, vsbackfill=vpbackfill/sqrt(3), rhobackfill=2000.0, 
	Qpbackfill=500.0, Qsbackfill=500.0,

	thick_backfill=0.4, thick_tueb=0.6, radius_steel=0.2;
	
	/* TBM */
	const float vp_tbm_steel=5050, vs_tbm_steel=vp_tbm_steel*sqrt(3), rho_tbm_steel=7700.0,
	Qp_tbm_steel=500.00, Qs_tbm_steel=500.0,
	vp_tbm_bentonit=1600.0, vs_tbm_bentonit=1.0e-4, rho_tbm_bentonit=1200.0,
	Qp_tbm_bentonit=500.00, Qs_tbm_bentonit=500.0,
	
	thick_schild=0.3, length_tbm=10.0, thick_cutterhead=0.4, spacing2formation=0.0, thick_bentchamber=2,
	radius_tbmaxis=0.5;
	
	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	ws=2.0*PI*FL[1];

	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
					
				x=(float)i*DX;
				y=(float)j*DY;
				z=(float)k*DZ;
				
				/* background */
				Vp=vp1; Vs=vs1; Rho=rho1;
				tp=2.0/Qp1; ts=2.0/Qs1;
				
				/* reflector */
                                if ((d<radius)&&(x>reflect_dist)){
                                        Vp=vp2; Vs=vs2; Rho=rho2;
                                        tp=2.0/Qp2; ts=2.0/Qs2;}
				
				/*excevation damage zone*/	
				/*working front*/
				
				d2=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z))+((x-length_tunnel+radius)*(x-length_tunnel+radius)));		
				
				if ((d2>radius)&&(d2<(radius+radiusedz))){ 
					Vp=vpedz+((vp1-vpedz)*(d2-radius)/radiusedz); 
					
					
					Vs=vsedz+((vs1-vsedz)*(d2-radius)/radiusedz);		
					
					Rho=rhoedz+((rho1-rhoedz)*(d2-radius)/radiusedz);
					
					tp=2.0/Qpedz; 
					ts=2.0/Qsedz;}		
																				
				/*excevation damage zone*/	
				/*around tube*/				
				
				d=sqrt(((y-axis_y)*(y-axis_y))+((z-axis_z)*(z-axis_z)));
				
				if ((d>radius)&&(d<(radius+radiusedz))&&(x<length_tunnel-radius/2)){ 
					Vp=vpedz+((vp1-vpedz)*(d-radius)/radiusedz); 
					
					
					Vs=vsedz+((vs1-vsedz)*(d-radius)/radiusedz);		
					
					Rho=rhoedz+((rho1-rhoedz)*(d-radius)/radiusedz);
					
					tp=2.0/Qpedz; 
					ts=2.0/Qsedz;}	
				
				/* tunnel */
				if ((d<radius)&&(x<length_tunnel)){ 
					Vp=vpt; Vs=vst; Rho=rhot;
					tp=2.0/Qpt; ts=2.0/Qst;} 
				
				/* tuebbing concrete */
				if ((d>radius-thick_tueb-thick_backfill)&&(d<(radius-thick_backfill))&&(x<length_tunnel-radius/2)){ 
					Vp=vpconcrete; 
					Vs=vsconcrete;

					Rho=rhoconcrete;
					
					tp=2.0/Qpconcrete; 
					ts=2.0/Qsconcrete;}
					
				/* concrete behind (backfill)*/
				if ((d>radius-thick_backfill)&&(d<(radius))&&(x<length_tunnel-radius/2)){ 
					Vp=vpbackfill; 
					Vs=vsbackfill;

					Rho=rhobackfill;
					
					tp=2.0/Qpbackfill; 
					ts=2.0/Qsbackfill;}


				/* tuebbing steel within concrete */
				for (l=0;l<360;l=l+15){
				
					dsteel=sqrt((y-axis_y+(radius-thick_backfill-0.5*thick_tueb)*sin(l))*(y-axis_y+(radius-thick_backfill-0.5*thick_tueb)*sin(l))+(z-axis_z+(radius-thick_backfill-0.5*thick_tueb)*cos(l))*(z-axis_z+(radius-thick_backfill-0.5*thick_tueb)*cos(l)));
					if ((dsteel<=radius_steel)&&(x<length_tunnel-radius/2)){
						Vp=vpsteel; 
						Vs=vssteel;

						Rho=rhosteel;
					
						tp=2.0/Qpsteel; 
						ts=2.0/Qssteel;
					}
				}
				
				/* TBM */
				/* remove Tuebbing*/
				if ((d<radius)&&(x<length_tunnel)&&(x>length_tunnel-length_tbm)){ 
					Vp=vpt; 
					Vs=vst;

					Rho=rhot;
					
					tp=2.0/Qpt; 
					ts=2.0/Qst;}			
					
				/*TBM */
				/* Schild tube*/
				if ((d>radius-thick_schild-spacing2formation)&&(d<(radius-spacing2formation))&&(x<length_tunnel-thick_cutterhead-thick_bentchamber)&&(x>length_tunnel-length_tbm)){ 
					Vp=vp_tbm_steel; 
					Vs=vs_tbm_steel;

					Rho=rho_tbm_steel;
					
					tp=2.0/Qp_tbm_steel; 
					ts=2.0/Qs_tbm_steel;}
				/*TBM */
				/* Schild cover tunnel*/
				if ((d<radius-spacing2formation)&&(x<length_tunnel-length_tbm+thick_schild)&&(x>length_tunnel-length_tbm)){ 
					Vp=vp_tbm_steel; 
					Vs=vs_tbm_steel;

					Rho=rho_tbm_steel;
					
					tp=2.0/Qp_tbm_steel; 
					ts=2.0/Qs_tbm_steel;}		
				/*TBM */
				/* Schild cover working face*/
				if ((d<radius-spacing2formation)&&(x<length_tunnel-thick_bentchamber)&&(x>length_tunnel-thick_bentchamber-thick_schild)){ 
					Vp=vp_tbm_steel; 
					Vs=vs_tbm_steel;

					Rho=rho_tbm_steel;
					
					tp=2.0/Qp_tbm_steel; 
					ts=2.0/Qs_tbm_steel;}	
					
				/*TBM */
				/* cutter head */
				if ((d<radius-spacing2formation)&&(x<length_tunnel-spacing2formation)&&(x>length_tunnel-thick_cutterhead-spacing2formation)){ 
					Vp=vp_tbm_steel; 
					Vs=vs_tbm_steel;

					Rho=rho_tbm_steel;
					
					tp=2.0/Qp_tbm_steel; 
					ts=2.0/Qs_tbm_steel;}
				/*TBM */
				/* axis connecting tbm schild body and cutter head */
				if ((d<radius_tbmaxis)&&(x<length_tunnel-spacing2formation-thick_cutterhead)&&(x>length_tunnel-thick_cutterhead-spacing2formation-thick_bentchamber)){ 
					Vp=vp_tbm_steel; 
					Vs=vs_tbm_steel;

					Rho=rho_tbm_steel;
					
					tp=2.0/Qp_tbm_steel; 
					ts=2.0/Qs_tbm_steel;}	
				
				
				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					taus[jj][ii][kk]=ts;
					taup[jj][ii][kk]=tp;
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}	


	
	sprintf(outfile,"model/tuebbingtbm.vs");
	writemod(outfile,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);
	
	
	sprintf(outfile,"model/tuebbingtbm.vp");
	writemod(outfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	/*sprintf(outfile,"model/test.u");
	writemod(outfile,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	*/
	sprintf(outfile,"model/tuebbingtbm.rho");
	writemod(outfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);
	
	/*
	sprintf(outfile,"model/test.pi");
	writemod(outfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(outfile,3);

	*/



	free_vector(pts,1,L);
	
	
}



