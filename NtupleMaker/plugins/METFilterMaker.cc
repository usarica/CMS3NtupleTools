#include <CMS3/NtupleMaker/interface/plugins/METFilterMaker.h>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <CMS3/NtupleMaker/interface/METFilterInfo.h>
#include <IvyFramework/IvyDataTools/interface/HelperFunctionsCore.h>


using namespace edm;
using namespace std;

 
METFilterMaker::METFilterMaker(const ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),

  doEcalFilterUpdate(iConfig.getParameter<bool>("doEcalFilterUpdate")),

  filtersInputTag_(iConfig.getParameter<InputTag>("filtersInputTag"))
{
  filtersTokenRECO = consumes<edm::TriggerResults>(edm::InputTag(filtersInputTag_.label(), "", "RECO"));
  filtersTokenPAT = consumes<edm::TriggerResults>(edm::InputTag(filtersInputTag_.label(), "", "PAT"));
  ecalBadCalibFilterUpdate_token = consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

  produces<METFilterInfo>();
}

METFilterMaker::~METFilterMaker(){}

void METFilterMaker::beginJob(){}
void METFilterMaker::endJob(){}

void METFilterMaker::produce(Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique<METFilterInfo>();

  edm::Handle<edm::TriggerResults> metFilterResultsH_;
  iEvent.getByToken(filtersTokenPAT, metFilterResultsH_);
  if (!metFilterResultsH_.isValid()) iEvent.getByToken(filtersTokenRECO, metFilterResultsH_);
  if (!metFilterResultsH_.isValid()) throw cms::Exception("METFilterMaker::produce: error getting TriggerResults PAT or RECO product from Event!");
  edm::TriggerNames metFilterNames_ = iEvent.triggerNames(*metFilterResultsH_);
  for (unsigned int i=0; i<metFilterNames_.size(); i++){
    std::string trigname = metFilterNames_.triggerName(i);
    if (trigname.find("Flag_") == std::string::npos) continue;

    HelperFunctions::replaceString<std::string, const char*>(trigname, "Flag_", "");
    result->flag_accept_map[trigname] = metFilterResultsH_->accept(i);
    /*
    if (metFilterNames_.triggerName(i) == "Flag_CSCTightHaloFilter")  idx_cscBeamHalo                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_CSCTightHalo2015Filter")  idx_cscBeamHalo2015                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_globalTightHalo2016Filter")  idx_globalTightHalo2016                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_globalSuperTightHalo2016Filter")  idx_globalSuperTightHalo2016                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_HBHENoiseFilter")  idx_hbheNoise                      = i;
    if (metFilterNames_.triggerName(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter")  idx_ecalTP                         = i;
    if (metFilterNames_.triggerName(i) == "Flag_hcalLaserEventFilter")  idx_hcalLaserEvent                 = i;
    if (metFilterNames_.triggerName(i) == "Flag_trackingFailureFilter")  idx_trackingFailure                = i;
    if (metFilterNames_.triggerName(i) == "Flag_chargedHadronTrackResolutionFilter")  idx_chargedHadronTrackResolution                = i;
    if (metFilterNames_.triggerName(i) == "Flag_eeBadScFilter")  idx_eeBadSc                        = i;
    if (metFilterNames_.triggerName(i) == "Flag_ecalLaserCorrFilter")  idx_ecalLaser                      = i;
    if (metFilterNames_.triggerName(i) == "Flag_METFilters")  idx_metfilter                      = i;
    if (metFilterNames_.triggerName(i) == "Flag_goodVertices")  idx_goodVertices                   = i;
    if (metFilterNames_.triggerName(i) == "Flag_trkPOGFilters")  idx_trkPOGFilters                  = i;
    if (metFilterNames_.triggerName(i) == "Flag_trkPOG_logErrorTooManyClusters")  idx_trkPOG_logErrorTooManyClusters = i;
    if (metFilterNames_.triggerName(i) == "Flag_trkPOG_manystripclus53X")  idx_trkPOG_manystripclus53X	    = i;
    if (metFilterNames_.triggerName(i) == "Flag_trkPOG_toomanystripclus53X")  idx_trkPOG_toomanystripclus53X     = i;
    if (metFilterNames_.triggerName(i) == "Flag_HBHENoiseIsoFilter") idx_hbheNoiseIso                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_CSCTightHaloTrkMuUnvetoFilter") idx_cscBeamHaloTrkMuUnveto          = i;
    if (metFilterNames_.triggerName(i) == "Flag_HcalStripHaloFilter") idx_hcalStrip                       = i;
    if (metFilterNames_.triggerName(i) == "Flag_EcalDeadCellBoundaryEnergyFilter") idx_ecalBoundaryEnergy              = i;
    if (metFilterNames_.triggerName(i) == "Flag_muonBadTrackFilter") idx_muonBadTrack                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_BadPFMuonFilter") idx_BadPFMuonFilter                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_BadChargedCandidateFilter") idx_BadChargedCandidateFilter                    = i;
    if (metFilterNames_.triggerName(i) == "Flag_ecalBadCalibFilter") idx_ecalBadCalibFilter                    = i;
    */
  }

  if (doEcalFilterUpdate){
    edm::Handle<bool> passEcalBadCalibFilterUpdate;
    iEvent.getByToken(ecalBadCalibFilterUpdate_token, passEcalBadCalibFilterUpdate);
    result->flag_accept_map["ecalBadCalibFilterUpdated"] = *passEcalBadCalibFilterUpdate;
  }
  else result->flag_accept_map["ecalBadCalibFilterUpdated"] = true;

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(METFilterMaker);
