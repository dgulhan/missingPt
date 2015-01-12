#include <iostream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

using namespace std;
static const int PI = 3.141592653589;
static const int ETA_CUT = 2;
static const int LEAD_PT_CUT = 120;
static const int SUBLEAD_PT_CUT = 50;
static const int DPHI_CUT = 50;
enum dataType{};

class utilities{
 public:
  static const double get_dphi( double phi1, double phi2) {
   double dphi = phi1 - phi2;
 
   if ( dphi > PI )
    dphi = dphi - 2. * PI;
   if ( dphi <= -PI ) 
    dphi = dphi + 2. * PI;
  
   if ( TMath::Abs(dphi) > PI ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
   }
  
   return dphi;
  } 
};