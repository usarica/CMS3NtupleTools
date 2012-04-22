// Includes 
#include "CMS2/NtupleMaker/interface/MetFilterMaker.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Namespaces
using namespace edm;
using namespace std;


// Constructor
MetFilterMaker::MetFilterMaker( const ParameterSet& iConfig ) {

  //
  branchprefix_ = "filt";
  aliasprefix_  = branchprefix_;  


  //
  ecalTPInputTag_                   = iConfig.getParameter<InputTag>("ecalTPInputTag"          );
  InputTag trackingFailureInputTag_ = iConfig.getParameter<InputTag>("trackingFailureInputTag" );


  //
  //produces <bool> ( branchprefix_ + "cscBeamHalo"     ).setBranchAlias( aliasprefix_ + "_cscBeamHalo"     );
  //produces <bool> ( branchprefix_ + "hbheNoise"       ).setBranchAlias( aliasprefix_ + "_hbheNoise"       );
  produces <bool> ( branchprefix_ + "ecalBE"          ).setBranchAlias( aliasprefix_ + "_ecalBE"          );
  produces <bool> ( branchprefix_ + "ecalDR"          ).setBranchAlias( aliasprefix_ + "_ecalDR"          );
  produces <bool> ( branchprefix_ + "ecalTP"          ).setBranchAlias( aliasprefix_ + "_ecalTP"          );
  produces <bool> ( branchprefix_ + "greedyMuon"      ).setBranchAlias( aliasprefix_ + "_greedyMuon"      );
  produces <bool> ( branchprefix_ + "hcalLaser"       ).setBranchAlias( aliasprefix_ + "_hcalLaser"       );
  produces <bool> ( branchprefix_ + "inconsistentMuon").setBranchAlias( aliasprefix_ + "_inconsistentMuon");
  produces <bool> ( branchprefix_ + "jetIDFailure"    ).setBranchAlias( aliasprefix_ + "_jetIDFailure"    );
  produces <bool> ( branchprefix_ + "multiEvent"      ).setBranchAlias( aliasprefix_ + "_multiEvent"      );
  produces <bool> ( branchprefix_ + "trackingFailure" ).setBranchAlias( aliasprefix_ + "_trackingFailure" );

} // End Constructor

MetFilterMaker::~MetFilterMaker () {}
void MetFilterMaker::beginJob   () {}
void MetFilterMaker::endJob     () {}

// Producer
void MetFilterMaker::produce( Event& iEvent, const edm::EventSetup& iSetup ) {
  
  //////////////
  // Pointers //
  //////////////

  //auto_ptr <bool> filt_cscBeamHalo      ( new bool(false) );
  //auto_ptr <bool> filt_hbheNoise        ( new bool(false) );
  auto_ptr <bool> filt_ecalBE           ( new bool(false) );
  auto_ptr <bool> filt_ecalDR           ( new bool(false) );
  auto_ptr <bool> filt_ecalTP           ( new bool(false) );
  auto_ptr <bool> filt_greedyMuon       ( new bool(false) );
  auto_ptr <bool> filt_hcalLaser        ( new bool(false) );
  auto_ptr <bool> filt_inconsistentMuon ( new bool(false) );
  auto_ptr <bool> filt_jetIDFailure     ( new bool(false) );
  auto_ptr <bool> filt_multiEvent       ( new bool(false) );
  auto_ptr <bool> filt_trackingFailure  ( new bool(false) );


  ////////////////////////
  // Assign MET Filters //
  ////////////////////////

  Handle<bool> b_ecalBE;
  Handle<bool> b_ecalDR;
  Handle<bool> b_ecalTP;
  Handle<bool> b_greedyMuon;
  Handle<bool> b_hcalLaser;
  Handle<bool> b_inconsistentMuon;
  Handle<bool> b_jetIDFailure;
  Handle<bool> b_multiEvent;
  Handle<bool> b_trackingFailure;

  iEvent.getByLabel( ecalTPInputTag_          , b_ecalTP          );
  iEvent.getByLabel( trackingFailureInputTag_ , b_trackingFailure );



  //
  if( b_ecalTP.isValid()          ){ 
    *filt_ecalTP          = *b_ecalTP;
  } 
  else {
    cout << "Invalid Handle: ecalTP.   File: " << __FILE__ << ", Line: " << __LINE__ << endl; 
  }   



  //
  if( b_trackingFailure.isValid() ){ 
    *filt_trackingFailure = *b_trackingFailure;
  } 
  else {
    cout << "Invalid Handle: trackingFailure" << __FILE__ << ", " << __LINE__ << endl; 
  }


  ///////////
  // Debug //
  ///////////

  //cout << ecalTPInputTag_ << ": " << b_ecalTP << endl;



  //////////////////
  // Write Output //
  //////////////////

  //iEvent.put( filt_cscBeamHalo      , branchprefix_ + "cscBeamHalo"     );
  //iEvent.put( filt_hbheNoise        , branchprefix_ + "hbheNoise"       );
  iEvent.put( filt_ecalBE           , branchprefix_ + "ecalBE"          );
  iEvent.put( filt_ecalDR           , branchprefix_ + "ecalDR"          );
  iEvent.put( filt_ecalTP           , branchprefix_ + "ecalTP"          );
  iEvent.put( filt_greedyMuon       , branchprefix_ + "greedyMuon"      );
  iEvent.put( filt_hcalLaser        , branchprefix_ + "hcalLaser"       );
  iEvent.put( filt_inconsistentMuon , branchprefix_ + "inconsistentMuon");
  iEvent.put( filt_jetIDFailure     , branchprefix_ + "jetIDFailure"    );
  iEvent.put( filt_multiEvent       , branchprefix_ + "multiEvent"      );
  iEvent.put( filt_trackingFailure  , branchprefix_ + "trackingFailure" );

} // End MetFilterMaker::produce()

//define this as a plug-in
DEFINE_FWK_MODULE( MetFilterMaker );
