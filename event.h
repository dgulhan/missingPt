#include "utilities.h"

enum Choice_jetpt_var{ GEN=0, RECO, FRAG, FRAG_RES, REF};

class Jet{
 private:
  double pt;
  double eta;
  double phi;
  double corr_pt;
  double rescorr_pt;
  double ref_pt;
 public:
  Jet(double pt, double eta, double phi, double corr_pt=0, double rescorr_pt=0, double ref_pt=0){
   this->pt=pt;
   this->eta=eta;
   this->phi=phi;
   this->corr_pt=corr_pt;
   this->rescorr_pt=rescorr_pt;
   this->ref_pt=ref_pt;
  }
  void set_corr_pt(double corr_pt){
   this->corr_pt=corr_pt;
  }  
  void set_rescorr_pt(double rescorr_pt){
   this->rescorr_pt=rescorr_pt;
  }  
  void set_ref_pt(double ref_pt){
   this->ref_pt=ref_pt;
  } 
  double get_pt(Choice_jetpt_var choice_pt) const{
   switch(choice_pt){
    case RECO: case GEN:
     return pt;
    case FRAG:
     return corr_pt;
    case FRAG_RES:
     return rescorr_pt;
    case REF:
     return ref_pt;
    default:
     cout<<"Invalid choice of jet pt for sorting"<<endl;
     break;
   }
  } 
  double get_eta(){
   return eta;
  }
  double get_phi(){
   return phi;
  } 
};

struct Jetpt_greater_than{
 private:
  Choice_jetpt_var choice_pt;
 public:
  Jetpt_greater_than(Choice_jetpt_var choice_pt){
   this->choice_pt=choice_pt;
  }
  bool operator()( const Jet& jet1, const Jet& jet2) const {
   return (jet1.get_pt(choice_pt) > jet2.get_pt(choice_pt));
  }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

class Event{
 private:
  bool   good_event;
  double mpt[11]; 
  bool gen_part;
  Choice_jetpt_var choice_pt;
  vector<Jet> jets;
  Jet *lead;
  Jet *sublead;
  Jet *third;
 public:
  void reset(){
   good_event = false;
   gen_part = false;
   for(int i=0; i<11; i++){
    mpt[i]=0;
   }
   lead = NULL;
   sublead = NULL;
   third = NULL;
  }
  Event(){
   reset();
  }
  
  //event selection related functions
  void set_event_sel(bool good_event){
   this->good_event = good_event;
  }
  bool is_good_event(){
   return good_event;
  }
  
  //jet related functions
  void set_choice_pt(Choice_jetpt_var choice_pt){
   this->choice_pt=choice_pt;
  }
  int get_njet(){
   return jets.size();
  }
  void add_jet(Jet jet){
   if(fabs(jet.get_eta()) < ETA_CUT){
    jets.push_back(jet);
   }
  } 
  Jet get_ieth_jet(int i){
   return jets[i];
  }
  void sort_jets(){
   sort( jets.begin(), jets.end(), Jetpt_greater_than(choice_pt));
   lead = new Jet(jets[0].get_pt(RECO), jets[0].get_eta(), jets[0].get_phi(),jets[0].get_pt(FRAG),jets[0].get_pt(FRAG_RES),jets[0].get_pt(REF));
   sublead = new Jet(jets[1].get_pt(RECO), jets[1].get_eta(), jets[1].get_phi(),jets[0].get_pt(FRAG),jets[0].get_pt(FRAG_RES),jets[0].get_pt(REF));
   if(jets.size()>2) third = new Jet(jets[2].get_pt(RECO), jets[2].get_eta(), jets[2].get_phi(),jets[0].get_pt(FRAG),jets[0].get_pt(FRAG_RES),jets[0].get_pt(REF));
  }
  Jet *get_lead(){
   if(get_njet()>0) return lead;
   else{
    cout<<"no leading reco jet"<<endl;
    exit(EXIT_FAILURE);
   }
  }
  Jet *get_sublead(){
   if(get_njet()>1) return sublead;
   else{
    cout<<"no subleading reco jet"<<endl;
    exit(EXIT_FAILURE);
   }
  }
  Jet *get_third(){
   if(get_njet()>2) return third;
   else{
    cout<<"no third reco jet"<<endl;
    exit(EXIT_FAILURE);
   }
  }
  double get_Aj(){
   double pt1, pt2;
   pt1=lead->get_pt(choice_pt);
   pt2=sublead->get_pt(choice_pt);   
   return (pt1-pt2)/(pt1+pt2);
  }
  double get_Aj23(){
   double pt2, pt3;
   pt2=sublead->get_pt(choice_pt);
   pt3=third->get_pt(choice_pt);   
   return (pt2-pt3)/(pt2+pt3);
  }
  double get_dphi12(){
   double phi1, phi2;
   phi1=lead->get_phi();
   phi2=sublead->get_phi();   
   return utilities::get_dphi(phi1,phi2);
  }
  double get_dphi13(){
   double phi1, phi3;
   phi1=lead->get_phi();
   phi3=sublead->get_phi();   
   return utilities::get_dphi(phi1,phi3);
  }
  //MPT related functions
  void set_gen_part(bool gen_part){
   this->gen_part=gen_part;
  }
  void set_mpt(double *mpt){
   for(int i=0;i<11;i++){
    this->mpt[i]=mpt[i];
   }
  }
  bool is_gen_part(){
   return gen_part;
  }
};