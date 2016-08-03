/*------------------------------------------------------------------------
 *   Source excitation using Moment tensor for simulation of earthquakes
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void eqsource(int nt, float *** sxx, float *** syy, float *** szz,
		float *** sxy, float *** syz, float *** sxz,
		float **  srcpos_loc, float ** signals, int nsrc, int * stype,
		float amon, float str, float dip, float rake){

	extern float DT, DX, DY, DZ;
	extern float AMON, STR, DIP, RAKE;
	int i, j, k, l;
	float amp, scale_amp;;
	float m11, m12, m13, m22, m23, m33;

	/* adding source wavelet to stress components
           (moment tensor source) at source points */

	amon=AMON; /*1.0e12;*/
	str=STR*PI/180.0;
	dip=DIP*PI/180.0;
	rake=RAKE*PI/180.0;

	for (l=1;l<=nsrc;l++)
	{
	  
		if(stype[l]==6){
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
			k=(int)srcpos_loc[3][l];

			amp=signals[l][nt];

			scale_amp=DT/(DX*DY*DZ);

			amp=amon*amp*scale_amp;


			m33 = -(sin(dip)*cos(rake)*sin(2.0*str)+
					sin(2.0*dip)*sin(rake)*sin(str)*sin(str));

			m13 = sin(dip)*cos(rake)*cos(2.*str)+
					0.5*(sin(2.*dip)*sin(rake)*sin(2.*str));

			m23 = -(cos(dip)*cos(rake)*cos(str)+
					cos(2.*dip)*sin(rake)*sin(str));

			m11 = sin(dip)*cos(rake)*sin(2.*str)-
					sin(2.*dip)*sin(rake)*cos(str)*cos(str);

			m12 = -(cos(dip)*cos(rake)*sin(str)-
					cos(2.*dip)*sin(rake)*cos(str));

			m22 = sin(2.*dip)*sin(rake);


			sxx[j][i][k]+=amp*m11;
			syy[j][i][k]+=amp*m22;
			szz[j][i][k]+=amp*m33;

			sxy[j][i][k]+=0.25*amp*m12;
			sxy[j][i-1][k]+=0.25*amp*m12;
			sxy[j-1][i][k]+=0.25*amp*m12;
			sxy[j-1][i-1][k]+=0.25*amp*m12;

			syz[j][i][k]+=0.25*amp*m23;
			syz[j][i][k-1]+=0.25*amp*m23;
			syz[j][i-1][k]+=0.25*amp*m23;
			syz[j][i-1][k-1]+=0.25*amp*m23;

			sxz[j][i][k]+=0.25*amp*m13;
			sxz[j-1][i][k]+=0.25*amp*m13;
			sxz[j][i][k-1]+=0.25*amp*m13;
			sxz[j-1][i][k-1]+=0.25*amp*m13;


		}
	}
}
