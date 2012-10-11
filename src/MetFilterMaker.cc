// Includes 
#include "CMS2/NtupleMaker/interface/MetFilterMaker.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Namespaces
using namespace edm;
using namespace std;

 
//
void checkValid( Handle<bool> hbool, InputTag tag ){
    if( !hbool.isValid() ){
        cout << "Error, Invalid Handle: " << tag <<  endl; 
    }
    return;
}


// Constructor
MetFilterMaker::MetFilterMaker( const ParameterSet& iConfig ) {

    //
    branchprefix_ = "filt";
    aliasprefix_  = branchprefix_;  


    //

    ecalBEInputTag_             = iConfig.getParameter<InputTag>("ecalBEInputTag"           );
    ecalDRInputTag_             = iConfig.getParameter<InputTag>("ecalDRInputTag"           );
    ecalTPInputTag_             = iConfig.getParameter<InputTag>("ecalTPInputTag"           );
    greedyMuonInputTag_         = iConfig.getParameter<InputTag>("greedyMuonInputTag"       );
    hcalLaserEventInputTag_     = iConfig.getParameter<InputTag>("hcalLaserEventInputTag"   );
    inconsistentMuonInputTag_   = iConfig.getParameter<InputTag>("inconsistentMuonInputTag" );
    jetIDFailureInputTag_       = iConfig.getParameter<InputTag>("jetIDFailureInputTag"     );
    multiEventFailureInputTag_  = iConfig.getParameter<InputTag>("multiEventFailureInputTag");
    trackingFailureInputTag_    = iConfig.getParameter<InputTag>("trackingFailureInputTag"  );
    eeBadScFilterInputTag_      = iConfig.getParameter<InputTag>("eeBadScFilterInputTag"    );
    ecalLaserCorrFilterInputTag_ = iConfig.getParameter<InputTag>("ecalLaserCorrFilterInputTag");


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
    produces <bool> ( branchprefix_ + "eeBadSc"         ).setBranchAlias( aliasprefix_ + "_eeBadSc"         );
    produces <bool> ( branchprefix_ + "ecalLaser"       ).setBranchAlias( aliasprefix_ + "_ecalLaser"       );
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
    auto_ptr <bool> filt_hcalLaserEvent   ( new bool(false) );
    auto_ptr <bool> filt_inconsistentMuon ( new bool(false) );
    auto_ptr <bool> filt_jetIDFailure     ( new bool(false) );
    auto_ptr <bool> filt_multiEventFailure( new bool(false) );
    auto_ptr <bool> filt_trackingFailure  ( new bool(false) );
    auto_ptr <bool> filt_eeBadSc          ( new bool(false) );
    auto_ptr <bool> filt_ecalLaser        ( new bool(false) );

    ////////////////////////
    // Assign MET Filters //
    ////////////////////////

    Handle<bool> b_ecalBE;
    Handle<bool> b_ecalDR;
    Handle<bool> b_ecalTP;
    Handle<bool> b_greedyMuon;
    Handle<bool> b_hcalLaserEvent;
    Handle<bool> b_inconsistentMuon;
    Handle<bool> b_jetIDFailure;
    Handle<bool> b_multiEventFailure;
    Handle<bool> b_trackingFailure;
    Handle<bool> b_eeBadSc;
    Handle<bool> b_ecalLaser;

    iEvent.getByLabel( ecalBEInputTag_            , b_ecalBE            );
    iEvent.getByLabel( ecalDRInputTag_            , b_ecalDR            );
    iEvent.getByLabel( ecalTPInputTag_            , b_ecalTP            );
    iEvent.getByLabel( greedyMuonInputTag_        , b_greedyMuon        );
    iEvent.getByLabel( hcalLaserEventInputTag_    , b_hcalLaserEvent    );
    iEvent.getByLabel( inconsistentMuonInputTag_  , b_inconsistentMuon  );
    iEvent.getByLabel( jetIDFailureInputTag_      , b_jetIDFailure      );
    iEvent.getByLabel( multiEventFailureInputTag_ , b_multiEventFailure );
    iEvent.getByLabel( trackingFailureInputTag_   , b_trackingFailure   );
    iEvent.getByLabel( eeBadScFilterInputTag_     , b_eeBadSc           );
    iEvent.getByLabel( ecalLaserCorrFilterInputTag_ , b_ecalLaser       );

    checkValid( b_ecalBE            , ecalBEInputTag_            );
    //checkValid( b_ecalDR            , ecalDRInputTag_            );
    checkValid( b_ecalTP            , ecalTPInputTag_            );
    checkValid( b_greedyMuon        , greedyMuonInputTag_        );
    checkValid( b_hcalLaserEvent    , hcalLaserEventInputTag_    );
    checkValid( b_inconsistentMuon  , inconsistentMuonInputTag_  );
    //checkValid( b_jetIDFailure      , jetIDFailureInputTag_      );
    checkValid( b_multiEventFailure , multiEventFailureInputTag_ );
    checkValid( b_trackingFailure   , trackingFailureInputTag_   );
    checkValid( b_eeBadSc           , eeBadScFilterInputTag_     );
    checkValid( b_ecalLaser         , ecalLaserCorrFilterInputTag_ );

    *filt_ecalBE            = *b_ecalBE;
    //*filt_ecalDR            = *b_ecalDR;
    *filt_ecalTP            = *b_ecalTP;
    *filt_greedyMuon        = *b_greedyMuon;
    *filt_hcalLaserEvent    = *b_hcalLaserEvent;
    *filt_inconsistentMuon  = *b_inconsistentMuon;
    //*filt_jetIDFailure      = *b_jetIDFailure;
    *filt_multiEventFailure = *b_multiEventFailure;
    *filt_trackingFailure   = *b_trackingFailure;
    *filt_eeBadSc           = *b_eeBadSc;
    *filt_ecalLaser         = *b_ecalLaser;

    //////////////////
    // Write Output //
    //////////////////

    //iEvent.put( filt_cscBeamHalo      , branchprefix_ + "cscBeamHalo"     );
    //iEvent.put( filt_hbheNoise        , branchprefix_ + "hbheNoise"       );
    iEvent.put( filt_ecalBE           , branchprefix_ + "ecalBE"          );
    iEvent.put( filt_ecalDR           , branchprefix_ + "ecalDR"          );
    iEvent.put( filt_ecalTP           , branchprefix_ + "ecalTP"          );
    iEvent.put( filt_greedyMuon       , branchprefix_ + "greedyMuon"      );
    iEvent.put( filt_hcalLaserEvent   , branchprefix_ + "hcalLaser"       );
    iEvent.put( filt_inconsistentMuon , branchprefix_ + "inconsistentMuon");
    iEvent.put( filt_jetIDFailure     , branchprefix_ + "jetIDFailure"    );
    iEvent.put( filt_multiEventFailure, branchprefix_ + "multiEvent"      );
    iEvent.put( filt_trackingFailure  , branchprefix_ + "trackingFailure" );
    iEvent.put( filt_eeBadSc          , branchprefix_ + "eeBadSc"         );
    iEvent.put( filt_ecalLaser        , branchprefix_ + "ecalLaser"       );

} // End MetFilterMaker::produce()

//define this as a plug-in
DEFINE_FWK_MODULE( MetFilterMaker );
