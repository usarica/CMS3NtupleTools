// Includes 
#include "CMS3/NtupleMaker/interface/MetFilterMaker.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
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
    aliasprefix_     = iConfig.getUntrackedParameter<string>("aliasPrefix");
    branchprefix_    = aliasprefix_;
    filtersInputTag_ = iConfig.getParameter<InputTag> ("filtersInputTag" );
    processName_     = iConfig.getUntrackedParameter<string> ("processName" );

    //
    produces <bool> ( branchprefix_ + "cscBeamHalo"                    ).setBranchAlias( aliasprefix_ + "_cscBeamHalo"                    );
    produces <bool> ( branchprefix_ + "hbheNoise"                      ).setBranchAlias( aliasprefix_ + "_hbheNoise"                      );
    produces <bool> ( branchprefix_ + "ecalTP"                         ).setBranchAlias( aliasprefix_ + "_ecalTP"                         );
    produces <bool> ( branchprefix_ + "hcalLaser"                      ).setBranchAlias( aliasprefix_ + "_hcalLaser"                      );
    produces <bool> ( branchprefix_ + "trackingFailure"                ).setBranchAlias( aliasprefix_ + "_trackingFailure"                );
    produces <bool> ( branchprefix_ + "eeBadSc"                        ).setBranchAlias( aliasprefix_ + "_eeBadSc"                        );
    produces <bool> ( branchprefix_ + "ecalLaser"                      ).setBranchAlias( aliasprefix_ + "_ecalLaser"                      );
    produces <bool> ( branchprefix_ + "metfilter"                      ).setBranchAlias( aliasprefix_ + "_metfilter"                      );
    produces <bool> ( branchprefix_ + "goodVertices"                   ).setBranchAlias( aliasprefix_ + "_goodVertices"                   );
    produces <bool> ( branchprefix_ + "trkPOGFilters"                  ).setBranchAlias( aliasprefix_ + "_trkPOGFilters"                  );
    produces <bool> ( branchprefix_ + "trkPOGlogErrorTooManyClusters"  ).setBranchAlias( aliasprefix_ + "_trkPOG_logErrorTooManyClusters" );
    produces <bool> ( branchprefix_ + "trkPOGmanystripclus53X"	       ).setBranchAlias( aliasprefix_ + "_trkPOG_manystripclus53X"        );
    produces <bool> ( branchprefix_ + "trkPOGtoomanystripclus53X"      ).setBranchAlias( aliasprefix_ + "_trkPOG_toomanystripclus53X"     );

    // For compatibility with CMS2 variable names
    produces <bool> ( "evtcscTightHaloId"                 ).setBranchAlias( "evt_cscTightHaloId"                 );
    produces <bool> ( "evthbheFilter"                     ).setBranchAlias( "evt_hbheFilter"                     );
	    

} // End Constructor

MetFilterMaker::~MetFilterMaker () {}
void MetFilterMaker::beginJob   () {}
void MetFilterMaker::endJob     () {}

// Producer
void MetFilterMaker::produce( Event& iEvent, const edm::EventSetup& iSetup ) {
  
  //////////////
  // Pointers //
  //////////////
  
  auto_ptr <bool> filt_cscBeamHalo                   ( new bool(false) );
  auto_ptr <bool> filt_hbheNoise                     ( new bool(false) );
  auto_ptr <bool> filt_ecalTP                        ( new bool(false) );
  auto_ptr <bool> filt_hcalLaserEvent                ( new bool(false) );
  auto_ptr <bool> filt_trackingFailure               ( new bool(false) );
  auto_ptr <bool> filt_eeBadSc                       ( new bool(false) );
  auto_ptr <bool> filt_ecalLaser                     ( new bool(false) );
  auto_ptr <bool> filt_metfilter                     ( new bool(false) );
  auto_ptr <bool> filt_goodVertices                  ( new bool(false) );
  auto_ptr <bool> filt_trkPOGFilters                 ( new bool(false) );
  auto_ptr <bool> filt_trkPOG_logErrorTooManyClusters( new bool(false) );
  auto_ptr <bool> filt_trkPOG_manystripclus53X       ( new bool(false) );
  auto_ptr <bool> filt_trkPOG_toomanystripclus53X    ( new bool(false) );
  // For compatibility with CMS2 variable names
  auto_ptr <bool> filt_cscTightHaloId                ( new bool(false) );
  auto_ptr <bool> filt_hbheFilter                    ( new bool(false) );

  ////////////////////////
  // Assign MET Filters //
  ////////////////////////
  
  iEvent.getByLabel(edm::InputTag(filtersInputTag_.label(), "", processName_), metFilterResultsH_);
  if (! metFilterResultsH_.isValid())
    throw cms::Exception("MetFilterMaker::produce: error getting TriggerResults_PAT product from Event!");
  
  int idx_cscBeamHalo                     = -1;
  int idx_hbheNoise                       = -1;
  int idx_ecalTP                          = -1;
  int idx_hcalLaserEvent                  = -1;
  int idx_trackingFailure                 = -1;
  int idx_eeBadSc                         = -1;
  int idx_ecalLaser                       = -1;
  int idx_metfilter                       = -1;
  int idx_goodVertices                    = -1;
  int idx_trkPOGFilters                   = -1;
  int idx_trkPOG_logErrorTooManyClusters  = -1;
  int idx_trkPOG_manystripclus53X	  = -1;
  int idx_trkPOG_toomanystripclus53X      = -1; 
  
  
  edm::TriggerNames metFilterNames_ = iEvent.triggerNames(*metFilterResultsH_); 
  for (unsigned int i=0; i<metFilterNames_.size(); i++) {
    if (  metFilterNames_.triggerName(i) == "Flag_CSCTightHaloFilter"	              )  idx_cscBeamHalo                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_HBHENoiseFilter"		      )  idx_hbheNoise                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter" )  idx_ecalTP                         = i;
    if (  metFilterNames_.triggerName(i) == "Flag_hcalLaserEventFilter"               )  idx_hcalLaserEvent                 = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trackingFailureFilter"	      )  idx_trackingFailure                = i;
    if (  metFilterNames_.triggerName(i) == "Flag_eeBadScFilter"		      )  idx_eeBadSc                        = i;
    if (  metFilterNames_.triggerName(i) == "Flag_ecalLaserCorrFilter"	              )  idx_ecalLaser                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_METFilters"                         )  idx_metfilter                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_goodVertices"                       )  idx_goodVertices                   = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOGFilters"                      )  idx_trkPOGFilters                  = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_logErrorTooManyClusters"     )  idx_trkPOG_logErrorTooManyClusters = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_manystripclus53X"            )  idx_trkPOG_manystripclus53X	    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_toomanystripclus53X"         )  idx_trkPOG_toomanystripclus53X     = i;
    //  std::cout << "metFilterName= " << metFilterNames_.triggerName(i) << std::endl;
  }

  
  *filt_cscBeamHalo                          = (idx_cscBeamHalo                    < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHalo                     );
  *filt_hbheNoise                            = (idx_hbheNoise                      < 0) ? false : metFilterResultsH_->accept(idx_hbheNoise                       );
  *filt_ecalTP                               = (idx_ecalTP                         < 0) ? false : metFilterResultsH_->accept(idx_ecalTP                          );
  *filt_hcalLaserEvent                       = (idx_hcalLaserEvent                 < 0) ? false : metFilterResultsH_->accept(idx_hcalLaserEvent                  );
  *filt_trackingFailure                      = (idx_trackingFailure                < 0) ? false : metFilterResultsH_->accept(idx_trackingFailure                 );
  *filt_eeBadSc                              = (idx_eeBadSc                        < 0) ? false : metFilterResultsH_->accept(idx_eeBadSc                         );
  *filt_ecalLaser                            = (idx_ecalLaser                      < 0) ? false : metFilterResultsH_->accept(idx_ecalLaser                       );
  *filt_metfilter                            = (idx_metfilter                      < 0) ? false : metFilterResultsH_->accept(idx_metfilter                       );
  *filt_goodVertices                         = (idx_goodVertices                   < 0) ? false : metFilterResultsH_->accept(idx_goodVertices                    );
  *filt_trkPOGFilters                        = (idx_trkPOGFilters                  < 0) ? false : metFilterResultsH_->accept(idx_trkPOGFilters                   );
  *filt_trkPOG_logErrorTooManyClusters       = (idx_trkPOG_logErrorTooManyClusters < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_logErrorTooManyClusters  );
  *filt_trkPOG_manystripclus53X	             = (idx_trkPOG_manystripclus53X        < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_manystripclus53X	 );
  *filt_trkPOG_toomanystripclus53X           = (idx_trkPOG_toomanystripclus53X     < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_toomanystripclus53X      );
  // For compatibility with CMS2 variable names
  *filt_cscTightHaloId                       = (idx_cscBeamHalo < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHalo                     );
  *filt_hbheFilter                           = (idx_hbheNoise < 0) ? false : metFilterResultsH_->accept(idx_hbheNoise                       );

  //////////////////
  // Write Output //
  //////////////////
  iEvent.put( filt_cscBeamHalo                    , branchprefix_ + "cscBeamHalo"                    );
  iEvent.put( filt_hbheNoise                      , branchprefix_ + "hbheNoise"                      );
  iEvent.put( filt_ecalTP                         , branchprefix_ + "ecalTP"                         );
  iEvent.put( filt_hcalLaserEvent                 , branchprefix_ + "hcalLaser"                      );
  iEvent.put( filt_trackingFailure                , branchprefix_ + "trackingFailure"                );
  iEvent.put( filt_eeBadSc                        , branchprefix_ + "eeBadSc"                        );
  iEvent.put( filt_ecalLaser                      , branchprefix_ + "ecalLaser"                      );
  iEvent.put( filt_metfilter                      , branchprefix_ + "metfilter"                      );
  iEvent.put( filt_goodVertices                   , branchprefix_ + "goodVertices"                   );
  iEvent.put( filt_trkPOGFilters                  , branchprefix_ + "trkPOGFilters"                  );
  iEvent.put( filt_trkPOG_logErrorTooManyClusters , branchprefix_ + "trkPOGlogErrorTooManyClusters"  );
  iEvent.put( filt_trkPOG_manystripclus53X	  , branchprefix_ + "trkPOGmanystripclus53X"	     );
  iEvent.put( filt_trkPOG_toomanystripclus53X     , branchprefix_ + "trkPOGtoomanystripclus53X"      );

  // For compatibility with CMS2 variable names
  iEvent.put( filt_cscTightHaloId                 , "evtcscTightHaloId"                  );
  iEvent.put( filt_hbheFilter                     , "evthbheFilter"                      );

  
} // End MetFilterMaker::produce()

//define this as a plug-in
DEFINE_FWK_MODULE( MetFilterMaker );
