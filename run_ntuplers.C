// #include "ntupler_data.C"
#include "ntupler_test.C"

void run_ntuplers(){
int npt=4;
double ptmin[]={0.5,1,2,4,  8,0.5};
double ptmax[]={  1,2,4,8,100,100};
for(int ipt=3;ipt<npt;ipt++){
 ntupler_test(ptmin[ipt],ptmax[ipt]);
 // ntupler_data(ptmin[ipt],ptmax[ipt]);
}
}