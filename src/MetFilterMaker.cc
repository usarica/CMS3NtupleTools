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
    processName_     = iConfig.getUntrackedParameter<string> ("processName" );
    filtersInputTag_ = iConfig.getParameter<InputTag> ("filtersInputTag" );
    filtersToken = consumes<edm::TriggerResults>(edm::InputTag(filtersInputTag_.label(), "", processName_)
);

    //
    produces <bool> ( branchprefix_ + "cscBeamHalo"                    ).setBranchAlias( aliasprefix_ + "_cscBeamHalo"                    );
    produces <bool> ( branchprefix_ + "globalTightHalo2016"            ).setBranchAlias( aliasprefix_ + "_globalTightHalo2016"            );
    produces <bool> ( branchprefix_ + "globalSuperTightHalo2016"       ).setBranchAlias( aliasprefix_ + "_globalSuperTightHalo2016"       );
    produces <bool> ( branchprefix_ + "cscBeamHalo2015"                ).setBranchAlias( aliasprefix_ + "_cscBeamHalo2015"                );
    produces <bool> ( branchprefix_ + "hbheNoise"                      ).setBranchAlias( aliasprefix_ + "_hbheNoise"                      );
    produces <bool> ( branchprefix_ + "ecalTP"                         ).setBranchAlias( aliasprefix_ + "_ecalTP"                         );
    produces <bool> ( branchprefix_ + "hcalLaser"                      ).setBranchAlias( aliasprefix_ + "_hcalLaser"                      );
    produces <bool> ( branchprefix_ + "trackingFailure"                ).setBranchAlias( aliasprefix_ + "_trackingFailure"                );
    produces <bool> ( branchprefix_ + "chargedHadronTrackResolution"   ).setBranchAlias( aliasprefix_ + "_chargedHadronTrackResolution"   );
    produces <bool> ( branchprefix_ + "eeBadSc"                        ).setBranchAlias( aliasprefix_ + "_eeBadSc"                        );
    produces <bool> ( branchprefix_ + "ecalLaser"                      ).setBranchAlias( aliasprefix_ + "_ecalLaser"                      );
    produces <bool> ( branchprefix_ + "metfilter"                      ).setBranchAlias( aliasprefix_ + "_metfilter"                      );
    produces <bool> ( branchprefix_ + "goodVertices"                   ).setBranchAlias( aliasprefix_ + "_goodVertices"                   );
    produces <bool> ( branchprefix_ + "trkPOGFilters"                  ).setBranchAlias( aliasprefix_ + "_trkPOGFilters"                  );
    produces <bool> ( branchprefix_ + "trkPOGlogErrorTooManyClusters"  ).setBranchAlias( aliasprefix_ + "_trkPOG_logErrorTooManyClusters" );
    produces <bool> ( branchprefix_ + "trkPOGmanystripclus53X"	       ).setBranchAlias( aliasprefix_ + "_trkPOG_manystripclus53X"        );
    produces <bool> ( branchprefix_ + "trkPOGtoomanystripclus53X"      ).setBranchAlias( aliasprefix_ + "_trkPOG_toomanystripclus53X"     );
    produces <bool> ( branchprefix_ + "hbheNoiseIso"                   ).setBranchAlias( aliasprefix_ + "_hbheNoiseIso"                   );
    produces <bool> ( branchprefix_ + "cscBeamHaloTrkMuUnveto"         ).setBranchAlias( aliasprefix_ + "_cscBeamHaloTrkMuUnveto"         );
    produces <bool> ( branchprefix_ + "hcalStrip"                      ).setBranchAlias( aliasprefix_ + "_hcalStrip"                      );
    produces <bool> ( branchprefix_ + "ecalBoundaryEnergy"             ).setBranchAlias( aliasprefix_ + "_ecalBoundaryEnergy"             );
    produces <bool> ( branchprefix_ + "muonBadTrack"                   ).setBranchAlias( aliasprefix_ + "_muonBadTrack"                   );
    produces <bool> ( branchprefix_ + "BadPFMuonFilter"                   ).setBranchAlias( aliasprefix_ + "_BadPFMuonFilter"                   );
    produces <bool> ( branchprefix_ + "BadChargedCandidateFilter"                   ).setBranchAlias( aliasprefix_ + "_BadChargedCandidateFilter"                   );
    produces <bool> ( branchprefix_ + "ecalBadCalibFilter"                   ).setBranchAlias( aliasprefix_ + "_ecalBadCalibFilter"                   );

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
  
  unique_ptr <bool> filt_cscBeamHalo                   ( new bool(false) );
  unique_ptr <bool> filt_cscBeamHalo2015               ( new bool(false) );
  unique_ptr <bool> filt_globalTightHalo2016           ( new bool(false) );
  unique_ptr <bool> filt_globalSuperTightHalo2016      ( new bool(false) );
  unique_ptr <bool> filt_hbheNoise                     ( new bool(false) );
  unique_ptr <bool> filt_ecalTP                        ( new bool(false) );
  unique_ptr <bool> filt_hcalLaserEvent                ( new bool(false) );
  unique_ptr <bool> filt_trackingFailure               ( new bool(false) );
  unique_ptr <bool> filt_chargedHadronTrackResolution  ( new bool(false) );
  unique_ptr <bool> filt_eeBadSc                       ( new bool(false) );
  unique_ptr <bool> filt_ecalLaser                     ( new bool(false) );
  unique_ptr <bool> filt_metfilter                     ( new bool(false) );
  unique_ptr <bool> filt_goodVertices                  ( new bool(false) );
  unique_ptr <bool> filt_trkPOGFilters                 ( new bool(false) );
  unique_ptr <bool> filt_trkPOG_logErrorTooManyClusters( new bool(false) );
  unique_ptr <bool> filt_trkPOG_manystripclus53X       ( new bool(false) );
  unique_ptr <bool> filt_trkPOG_toomanystripclus53X    ( new bool(false) );
  unique_ptr <bool> filt_hbheNoiseIso                  ( new bool(false) );
  unique_ptr <bool> filt_cscBeamHaloTrkMuUnveto        ( new bool(false) );
  unique_ptr <bool> filt_hcalStrip                     ( new bool(false) );
  unique_ptr <bool> filt_ecalBoundaryEnergy            ( new bool(false) );
  unique_ptr <bool> filt_muonBadTrack                  ( new bool(false) );
  unique_ptr <bool> filt_BadPFMuonFilter                  ( new bool(false) );
  unique_ptr <bool> filt_BadChargedCandidateFilter                  ( new bool(false) );
  unique_ptr <bool> filt_ecalBadCalibFilter                  ( new bool(false) );



  // For compatibility with CMS2 variable names
  unique_ptr <bool> filt_cscTightHaloId                ( new bool(false) );
  unique_ptr <bool> filt_hbheFilter                    ( new bool(false) );

  ////////////////////////
  // Assign MET Filters //
  ////////////////////////
  
  iEvent.getByToken(filtersToken, metFilterResultsH_);
  if (! metFilterResultsH_.isValid())
    throw cms::Exception("MetFilterMaker::produce: error getting TriggerResults_PAT product from Event!");
  
  int idx_cscBeamHalo                     = -1;
  int idx_cscBeamHalo2015                 = -1;
  int idx_globalTightHalo2016             = -1;
  int idx_globalSuperTightHalo2016        = -1;
  int idx_hbheNoise                       = -1;
  int idx_ecalTP                          = -1;
  int idx_hcalLaserEvent                  = -1;
  int idx_trackingFailure                 = -1;
  int idx_chargedHadronTrackResolution    = -1;
  int idx_eeBadSc                         = -1;
  int idx_ecalLaser                       = -1;
  int idx_metfilter                       = -1;
  int idx_goodVertices                    = -1;
  int idx_trkPOGFilters                   = -1;
  int idx_trkPOG_logErrorTooManyClusters  = -1;
  int idx_trkPOG_manystripclus53X	      = -1;
  int idx_trkPOG_toomanystripclus53X      = -1; 
  int idx_hbheNoiseIso                    = -1;
  int idx_cscBeamHaloTrkMuUnveto          = -1;
  int idx_hcalStrip                       = -1;
  int idx_ecalBoundaryEnergy              = -1;
  int idx_muonBadTrack                    = -1;
  int idx_BadPFMuonFilter                 = -1;
  int idx_BadChargedCandidateFilter       = -1;
  int idx_ecalBadCalibFilter              = -1;



  
  
  edm::TriggerNames metFilterNames_ = iEvent.triggerNames(*metFilterResultsH_); 
  for (unsigned int i=0; i<metFilterNames_.size(); i++) {
    if (  metFilterNames_.triggerName(i) == "Flag_CSCTightHaloFilter"	              )  idx_cscBeamHalo                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_CSCTightHalo2015Filter"	              )  idx_cscBeamHalo2015                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_globalTightHalo2016Filter"	              )  idx_globalTightHalo2016                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_globalSuperTightHalo2016Filter"	              )  idx_globalSuperTightHalo2016                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_HBHENoiseFilter"		      )  idx_hbheNoise                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter" )  idx_ecalTP                         = i;
    if (  metFilterNames_.triggerName(i) == "Flag_hcalLaserEventFilter"               )  idx_hcalLaserEvent                 = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trackingFailureFilter"	      )  idx_trackingFailure                = i;
    if (  metFilterNames_.triggerName(i) == "Flag_chargedHadronTrackResolutionFilter"	      )  idx_chargedHadronTrackResolution                = i;
    if (  metFilterNames_.triggerName(i) == "Flag_eeBadScFilter"		      )  idx_eeBadSc                        = i;
    if (  metFilterNames_.triggerName(i) == "Flag_ecalLaserCorrFilter"	              )  idx_ecalLaser                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_METFilters"                         )  idx_metfilter                      = i;
    if (  metFilterNames_.triggerName(i) == "Flag_goodVertices"                       )  idx_goodVertices                   = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOGFilters"                      )  idx_trkPOGFilters                  = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_logErrorTooManyClusters"     )  idx_trkPOG_logErrorTooManyClusters = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_manystripclus53X"            )  idx_trkPOG_manystripclus53X	    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_trkPOG_toomanystripclus53X"         )  idx_trkPOG_toomanystripclus53X     = i;
    if (  metFilterNames_.triggerName(i) == "Flag_HBHENoiseIsoFilter"                 ) idx_hbheNoiseIso                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_CSCTightHaloTrkMuUnvetoFilter"      ) idx_cscBeamHaloTrkMuUnveto          = i;
    if (  metFilterNames_.triggerName(i) == "Flag_HcalStripHaloFilter"                ) idx_hcalStrip                       = i;
    if (  metFilterNames_.triggerName(i) == "Flag_EcalDeadCellBoundaryEnergyFilter"   ) idx_ecalBoundaryEnergy              = i;
    if (  metFilterNames_.triggerName(i) == "Flag_muonBadTrackFilter"                 ) idx_muonBadTrack                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_BadPFMuonFilter"                 ) idx_BadPFMuonFilter                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_BadChargedCandidateFilter"                 ) idx_BadChargedCandidateFilter                    = i;
    if (  metFilterNames_.triggerName(i) == "Flag_ecalBadCalibFilter"                 ) idx_ecalBadCalibFilter                    = i;
  }

  
  *filt_cscBeamHalo                          = (idx_cscBeamHalo                    < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHalo                     );
  *filt_cscBeamHalo2015                      = (idx_cscBeamHalo2015                < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHalo2015                 );
  *filt_globalTightHalo2016                  = (idx_globalTightHalo2016            < 0) ? false : metFilterResultsH_->accept(idx_globalTightHalo2016             );
  *filt_globalSuperTightHalo2016             = (idx_globalSuperTightHalo2016       < 0) ? false : metFilterResultsH_->accept(idx_globalSuperTightHalo2016        );
  *filt_hbheNoise                            = (idx_hbheNoise                      < 0) ? false : metFilterResultsH_->accept(idx_hbheNoise                       );
  *filt_ecalTP                               = (idx_ecalTP                         < 0) ? false : metFilterResultsH_->accept(idx_ecalTP                          );
  *filt_hcalLaserEvent                       = (idx_hcalLaserEvent                 < 0) ? false : metFilterResultsH_->accept(idx_hcalLaserEvent                  );
  *filt_trackingFailure                      = (idx_trackingFailure                < 0) ? false : metFilterResultsH_->accept(idx_trackingFailure                 );
  *filt_chargedHadronTrackResolution         = (idx_chargedHadronTrackResolution   < 0) ? false : metFilterResultsH_->accept(idx_chargedHadronTrackResolution    );
  *filt_eeBadSc                              = (idx_eeBadSc                        < 0) ? false : metFilterResultsH_->accept(idx_eeBadSc                         );
  *filt_ecalLaser                            = (idx_ecalLaser                      < 0) ? false : metFilterResultsH_->accept(idx_ecalLaser                       );
  *filt_metfilter                            = (idx_metfilter                      < 0) ? false : metFilterResultsH_->accept(idx_metfilter                       );
  *filt_goodVertices                         = (idx_goodVertices                   < 0) ? false : metFilterResultsH_->accept(idx_goodVertices                    );
  *filt_trkPOGFilters                        = (idx_trkPOGFilters                  < 0) ? false : metFilterResultsH_->accept(idx_trkPOGFilters                   );
  *filt_trkPOG_logErrorTooManyClusters       = (idx_trkPOG_logErrorTooManyClusters < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_logErrorTooManyClusters  );
  *filt_trkPOG_manystripclus53X	             = (idx_trkPOG_manystripclus53X        < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_manystripclus53X	     );
  *filt_trkPOG_toomanystripclus53X           = (idx_trkPOG_toomanystripclus53X     < 0) ? false : metFilterResultsH_->accept(idx_trkPOG_toomanystripclus53X      );
  *filt_hbheNoiseIso                         = (idx_hbheNoiseIso                   < 0) ? false : metFilterResultsH_->accept(idx_hbheNoiseIso                    );
  *filt_cscBeamHaloTrkMuUnveto               = (idx_cscBeamHaloTrkMuUnveto         < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHaloTrkMuUnveto          );
  *filt_hcalStrip                            = (idx_hcalStrip                      < 0) ? false : metFilterResultsH_->accept(idx_hcalStrip                       );
  *filt_ecalBoundaryEnergy                   = (idx_ecalBoundaryEnergy             < 0) ? false : metFilterResultsH_->accept(idx_ecalBoundaryEnergy              );
  *filt_muonBadTrack                         = (idx_muonBadTrack                   < 0) ? false : metFilterResultsH_->accept(idx_muonBadTrack                    );
  *filt_BadPFMuonFilter                         = (idx_BadPFMuonFilter                   < 0) ? false : metFilterResultsH_->accept(idx_BadPFMuonFilter                    );
  *filt_BadChargedCandidateFilter                         = (idx_BadChargedCandidateFilter                   < 0) ? false : metFilterResultsH_->accept(idx_BadChargedCandidateFilter                    );
  *filt_ecalBadCalibFilter                         = (idx_ecalBadCalibFilter                   < 0) ? false : metFilterResultsH_->accept(idx_ecalBadCalibFilter                    );

  // For compatibility with CMS2 variable names
  *filt_cscTightHaloId                       = (idx_cscBeamHalo < 0) ? false : metFilterResultsH_->accept(idx_cscBeamHalo                     );
  *filt_hbheFilter                           = (idx_hbheNoise < 0) ? false : metFilterResultsH_->accept(idx_hbheNoise                       );

  //////////////////
  // Write Output //
  //////////////////

  iEvent.put(std::move( filt_cscBeamHalo                    ), branchprefix_ + "cscBeamHalo"                    );
  iEvent.put(std::move( filt_cscBeamHalo2015                ), branchprefix_ + "cscBeamHalo2015"                );
  iEvent.put(std::move( filt_globalTightHalo2016            ), branchprefix_ + "globalTightHalo2016"            );
  iEvent.put(std::move( filt_globalSuperTightHalo2016       ), branchprefix_ + "globalSuperTightHalo2016"       );
  iEvent.put(std::move( filt_hbheNoise                      ), branchprefix_ + "hbheNoise"                      );
  iEvent.put(std::move( filt_ecalTP                         ), branchprefix_ + "ecalTP"                         );
  iEvent.put(std::move( filt_hcalLaserEvent                 ), branchprefix_ + "hcalLaser"                      );
  iEvent.put(std::move( filt_trackingFailure                ), branchprefix_ + "trackingFailure"                );
  iEvent.put(std::move( filt_chargedHadronTrackResolution   ), branchprefix_ + "chargedHadronTrackResolution"   );
  iEvent.put(std::move( filt_eeBadSc                        ), branchprefix_ + "eeBadSc"                        );
  iEvent.put(std::move( filt_ecalLaser                      ), branchprefix_ + "ecalLaser"                      );
  iEvent.put(std::move( filt_metfilter                      ), branchprefix_ + "metfilter"                      );
  iEvent.put(std::move( filt_goodVertices                   ), branchprefix_ + "goodVertices"                   );
  iEvent.put(std::move( filt_trkPOGFilters                  ), branchprefix_ + "trkPOGFilters"                  );
  iEvent.put(std::move( filt_trkPOG_logErrorTooManyClusters ), branchprefix_ + "trkPOGlogErrorTooManyClusters"  );
  iEvent.put(std::move( filt_trkPOG_manystripclus53X	      ), branchprefix_ + "trkPOGmanystripclus53X"	     );
  iEvent.put(std::move( filt_trkPOG_toomanystripclus53X     ), branchprefix_ + "trkPOGtoomanystripclus53X"      );
  iEvent.put(std::move( filt_hbheNoiseIso                   ), branchprefix_ + "hbheNoiseIso"                   );
  iEvent.put(std::move( filt_cscBeamHaloTrkMuUnveto         ), branchprefix_ + "cscBeamHaloTrkMuUnveto"         );
  iEvent.put(std::move( filt_hcalStrip                      ), branchprefix_ + "hcalStrip"                      );
  iEvent.put(std::move( filt_ecalBoundaryEnergy             ), branchprefix_ + "ecalBoundaryEnergy"             );
  iEvent.put(std::move( filt_muonBadTrack                   ), branchprefix_ + "muonBadTrack"                   );
  iEvent.put(std::move( filt_BadPFMuonFilter                ), branchprefix_ + "BadPFMuonFilter"                );
  iEvent.put(std::move( filt_BadChargedCandidateFilter      ), branchprefix_ + "BadChargedCandidateFilter"      );
  iEvent.put(std::move( filt_ecalBadCalibFilter             ), branchprefix_ + "ecalBadCalibFilter"             );

  // For compatibility with CMS2 variable names
  iEvent.put(std::move( filt_cscTightHaloId                 ), "evtcscTightHaloId"                  );
  iEvent.put(std::move( filt_hbheFilter                     ), "evthbheFilter"                      );
  
} // End MetFilterMaker::produce()

//define this as a plug-in
DEFINE_FWK_MODULE( MetFilterMaker );
