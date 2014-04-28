#include "integrated.C"

void run(){
int ncent=4;
int ngen=3;
int neta=1;
bool doMC=false;
int centmin[]={0,20,60,100};
int centmax[]={20,60,100,200};

int genreco[]={3,4,5};

double etadijet[]={2};
for(int icent=0;icent<ncent;icent++){
 for(int igen=0;igen<ngen;igen++){
   for(int ieta=0;ieta<neta;ieta++){
    if(!(!doMC && (genreco[igen]==3 || genreco[igen]==4))) continue; 
    integrated(centmin[icent],centmax[icent],genreco[igen],etadijet[ieta]);
   }
 }
}
}