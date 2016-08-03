/*------------------------------------------------------------------------
 *   Data Conversion Program from Binary to VTK-Format
 *   Daniel Koehn
 *
 *   Raisdorf, the 22nd of february 2004
 *   update: 29.3.2004 - included VTK Binary IO
 *           01.4.2004 - added Decomposition of Data-Set in smaller 
 *                       blocks
 *           28.10.2007 - using this code to cut a plane in x-y-direction out
 *                        of the 3D block
 *
 *  ----------------------------------------------------------------------*/
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

int main( int argc, char *argv[] )
{

/* define input and output files here */
ifstream input("./sofi3D/par/model/random.sofi3D.rho",ios::in | ios::binary);
ofstream output("./sofi3D/par/model/random_shot2_cut.rho",ios::out | ios::binary);
int i,j,k,l,h,NX,NY,NZ;
int IDX, IDY, IDZ, NSNAP;
int CUTPOINT;
float vp, DH;
char test[50];

if(input.fail()){
    cout << "Sorry, can't read that file" << endl;
    exit(1);
    }

/* define geometry (like in *.inp)*/
NX=110;
NY=320;
NZ=440;

IDX=1;
IDY=1;
IDZ=1;

/* number of time steps */
NSNAP=1;

/* define gridpoint in z-direction where the cut plane is located */
CUTPOINT=321;

for(k=1;k<=NSNAP;k++){
cout << "Cutting wavefield time step ... " << k << endl;

    for(l=1;l<=NZ;l+=IDZ){
      for (i=1;i<=NX;i+=IDX){
        for (j=1;j<=NY;j+=IDY){
	
	  input.read((char*)&vp,sizeof(float));
	
	  if(j==1){
	    output.write((char*)&vp,sizeof(float));
	  }
	  
        }
     }
    
   }
}
input.close();
output.close();
return 0;
}
