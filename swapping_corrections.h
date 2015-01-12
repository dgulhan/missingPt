class Swap_corr{
 private:
  double           radius;
  dataType         type;
  int              ncent;
  static const int cent_min[5];
  static const int cent_max[5];
  TFile            *f_swap;
  TF1              *fit[5];
  
 public:
  Swap_corr(int table, int radius){
   this->radius=radius;
   this->table=table;
   
   if(!(radius==2 || radius==3 || radius==4 || radius==5)){
    cout<<"second parameter radius has to be 2,3,4 or 5, R=radius/10"<<endl;
   }
   
   switch(type){
    case tPbPb: 
     f_swap = new TFile(Form("SwappingCorrections/swapPlot_v4_R%d.root",radius));
	   ncent=5;
     break;
    case tPP:
	   f_swap = new TFile(Form("SwappingCorrections/swapPlot_PP_v4_R%d.root",radius));
     ncent=1;
 	   break;
    default:
	   cout<<"table can take the values 0 1 2"<<endl;
	   break;
    }
    for(int icent=0; icent<ncent; icent++){
     fit[icent]=(TF1*)f_swap->Get(Form("fPLeading_%d_%d",cent_min[icent],cent_max[icent]));
    }
   }
 
   double get_swap_corr(double Aj, int cent=0){
    int cent_bin = 0;
  
    if(table>0) cent=0;
  
    for(int icent=0; icent<ncent; icent++){
     if(cent>=cent_min[icent] && cent<cent_max[icent]) cent_bin=icent;
    }
  
    double q = fit[cent_bin]->Eval(Aj);
    double p=1-q;
    return (1-2*p);
   } 
};


 const int Swap_corr::cent_min[5]={ 0,20, 60,100,140};
 const int Swap_corr::cent_max[5]={20,60,100,140,200};