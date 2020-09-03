#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/MuonReco/interface/MuonSelectors.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>

#include <CMS3/NtupleMaker/interface/MuonSelectionHelpers.h>


namespace MuonSelectionHelpers{

  float muonEffArea(pat::Muon const& obj, int const& year){
    double eta = std::abs(obj.eta());
    float ea=-1;
    if (year==2016){
      // From https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt
      if (eta<=0.8) ea = 0.0735;
      else if (eta<=1.3) ea = 0.0619;
      else if (eta<=2.) ea = 0.0465;
      else if (eta<=2.2) ea = 0.0433;
      else ea = 0.0577;
    }
    else if (year==2017 || 2018){
      // From https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt
      if (eta<=0.8) ea = 0.0566;
      else if (eta<=1.3) ea = 0.0562;
      else if (eta<=2.) ea = 0.0363;
      else if (eta<=2.2) ea = 0.0119;
      else ea = 0.0064;
    }
    else cms::Exception("UnknownYear") << "MuonSelectionHelpers::muonEffArea: Year " << year << " is not implemented!" << std::endl;

    return ea;
  }

  float muonPFIsoComb(pat::Muon const& obj, int const& year, MuonSelectionHelpers::IsolationType const& type, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr, double* sum_neutral_EAcorr_nofsr){
    reco::MuonPFIsolation const* pfStruct = nullptr;

    double dR = 0;
    if (type==PFIso03){
      pfStruct = &(obj.pfIsolationR03());
      dR = 0.3;
    }
    else if (type==PFIso04){
      pfStruct = &(obj.pfIsolationR04());
      dR = 0.4;
    }
    else cms::Exception("UnknownIsoDR") << "MuonSelectionHelpers::muonPFIsoComb: Type " << type << " is not implemented!" << std::endl;

    double ch = pfStruct->sumChargedHadronPt;
    double nh = pfStruct->sumNeutralHadronEt;
    double em = pfStruct->sumPhotonEt;
    double db = pfStruct->sumPUPt;

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - 0.5*db;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    double ea = MuonSelectionHelpers::muonEffArea(obj, year);
    ea *= std::pow(dR / 0.3, 2);
    double sum_neutral_EAcorr_nofsr_val = nh + em - rho * ea;
    if (sum_neutral_EAcorr_nofsr) *sum_neutral_EAcorr_nofsr = sum_neutral_EAcorr_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }
  float muonMiniIsoComb(pat::Muon const& obj, int const& year, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one

    pat::PFIsolation const& miniiso = obj.miniPFIsolation();
    double ch = miniiso.chargedHadronIso();
    double nh = miniiso.neutralHadronIso();
    double em = miniiso.photonIso();

    double ea = MuonSelectionHelpers::muonEffArea(obj, year);
    double dR = 10. / std::min(std::max(uncorr_pt, 50.), 200.);
    ea *= std::pow(dR / 0.3, 2);

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - rho * ea;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }

  bool testUpdatedHighPtMuon(pat::Muon const& obj, reco::Vertex const& vtx, int const& year){
    if (year<2016 || year>2018) return ((obj.selectors() & reco::Muon::CutBasedIdGlobalHighPt) == reco::Muon::CutBasedIdGlobalHighPt);
    //
    // Below is a direct copy from CMSSW_10_4_X/DataFormats/MuonReco/src/MuonSelectors.cc::isHighPtMuon
    //
#if CMSSW_VERSION_MAJOR>10 || (CMSSW_VERSION_MAJOR==10 && CMSSW_VERSION_MINOR>=4)
    return muon::isHighPtMuon(obj);
#else
    if (!obj.isGlobalMuon()) return false;

    bool muValHits = (
      obj.globalTrack()->hitPattern().numberOfValidMuonHits()>0
      ||
      obj.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0
      );

    bool muMatchedSt = obj.numberOfMatchedStations()>1;
    if (!muMatchedSt){
      if (obj.isTrackerMuon() && obj.numberOfMatchedStations()==1){
        if (
          obj.expectedNnumberOfMatchedStations()<2
          ||
          !(obj.stationMask()==1 || obj.stationMask()==16)
          ||
          obj.numberOfMatchedRPCLayers()>2
          )
          muMatchedSt = true;
      }
    }

    bool muID = muValHits && muMatchedSt;

    bool hits =
      obj.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
      &&
      obj.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

    bool momQuality = obj.tunePMuonBestTrack()->ptError()/obj.tunePMuonBestTrack()->pt() < 0.3;

    bool ip = fabs(obj.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(obj.innerTrack()->dz(vtx.position())) < 0.5;

    return muID && hits && momQuality && ip;
#endif
  }

  bool testGoodMETPFMuon(pat::PackedCandidate const& pfcand){
    // The following selection requirements come from process.basicJetsForMetModifiedMET [of type EDProducer("PATJetCleanerForType1MET")]
    //   skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    //   skipMuons = cms.bool(True),
    return (pfcand.isGlobalMuon() || pfcand.isStandAloneMuon());
  }

  bool testLooseTriggerId(pat::Muon const& obj, int const& /*year*/){
#if CMSSW_VERSION_MAJOR>=10 
    return muon::isLooseTriggerMuon(obj);
#else
    bool pass_looseTriggerId = false;
    const reco::TrackRef innerTrack = obj.innerTrack();
    if (innerTrack.isNonnull()){
      bool tk_id = muon::isGoodMuon(obj, muon::TMOneStationTight);
      if (tk_id){
        bool layer_requirements =
          innerTrack->hitPattern().trackerLayersWithMeasurement() > 5 &&
          innerTrack->hitPattern().pixelLayersWithMeasurement() > 0;
        bool match_requirements = (obj.expectedNnumberOfMatchedStations()<2) || (obj.numberOfMatchedStations()>1) || (obj.pt()<8.); // This pT has to be the uncorrected pT.
        pass_looseTriggerId = (layer_requirements && match_requirements);
      }
    }
    return pass_looseTriggerId;
#endif
  }
  bool testTightCharge(pat::Muon const& obj, int const& year){
    //const reco::TrackRef chosenTrack = obj.innerTrack(); // Used by SNT
    const reco::TrackRef chosenTrack = obj.muonBestTrack();
    if (chosenTrack.isNonnull()){
      double trk_pt_err = chosenTrack->ptError();
      double trk_pt = chosenTrack->pt();
      return (trk_pt_err/trk_pt < tightCharge_pt_err_rel_thr);
    }
    else return false;
  }
  bool testMuonTiming(pat::Muon const& obj, int const& /*year*/){
    // Cut suggestions from Piotr for out-of-time muons from https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf
    // reco::Muon::InTimeMuon selector bit flag also stores the same info, but it is only set for 10.2.X+.
    // The following code is a copy of MuonSelectors.cc::outOfTimeMuon. Why not call that directly? Because the header does not include this function!
    auto const& combinedTime = obj.time();
    auto const& rpcTime = obj.rpcTime();
    bool combinedTimeIsOk = (combinedTime.nDof>7);
    bool rpcTimeIsOk = (rpcTime.nDof>1 && std::abs(rpcTime.timeAtIpInOutErr)<0.001);
    bool outOfTime = false;
    if (rpcTimeIsOk) outOfTime = (
      std::abs(rpcTime.timeAtIpInOut)>10.
      &&
      !(combinedTimeIsOk && std::abs(combinedTime.timeAtIpInOut)<10.)
      );
    else outOfTime = (combinedTimeIsOk && (combinedTime.timeAtIpInOut>20. || combinedTime.timeAtIpInOut<-45.));
    return !outOfTime;
  }

  bool testProbeMuonForTnP(pat::Muon const& obj, int const& /*year*/){
    return (obj.innerTrack().isNonnull());
  }
  bool testProbeMuonSTAForTnP(pat::Muon const& obj, int const& /*year*/){
    return (obj.outerTrack().isNonnull());
  }

  bool testSkimMuon(pat::Muon const& obj, int const& /*year*/){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    bool passAnyPOGBit = (obj.userInt("POG_selector_bits")!=0);
    bool passTnP = ((obj.userInt("is_probeForTnP")+obj.userInt("is_probeForTnP_STA"))>0);
    return (
      (passTnP || passAnyPOGBit) &&
      eta<selection_skim_eta && (
        uncorr_pt>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_scale_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_scale_totalDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_smear_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_smear_totalDn")>=selection_skim_pt
        )
      );
  }

}
