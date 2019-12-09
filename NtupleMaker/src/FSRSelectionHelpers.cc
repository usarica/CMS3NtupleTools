#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/FSRSelectionHelpers.h>


namespace FSRSelectionHelpers{

  bool testSCVeto(pat::PackedCandidate const* pfcand, pat::Electron const* electron){
    bool SCVeto = false;
    auto electron_pfcands = electron->associatedPackedPFCandidates();
    if (!electron_pfcands.empty()){
      for (auto electron_pfcand:electron_pfcands){
        if (pfcand == &(*electron_pfcand)){
          SCVeto = true;
          break;
        }
      }
    }
    /*
    // Prior to 2016, this was used:
    double dR = reco::deltaR(*(electron->superCluster()), *pfcand);
    if (
      (std::abs(pfcand->eta() - electron->superCluster()->eta())<0.05 && std::abs(reco::deltaPhi(pfcand->phi(), electron->superCluster()->phi()))<2.)
      ||
      dR<0.15
      ){
      SCVeto=true;
      break;
    }
    */
    return SCVeto;
  }
  bool testSCVeto(pat::PackedCandidate const* pfcand, pat::Photon const* photon){
    bool SCVeto = false;
    auto photon_pfcands = photon->associatedPackedPFCandidates();
    if (!photon_pfcands.empty()){
      for (auto photon_pfcand:photon_pfcands){
        if (pfcand == &(*photon_pfcand)){
          SCVeto = true;
          break;
        }
      }
    }
    return SCVeto;
  }

  bool testSkimFSR_PtEta(pat::PackedCandidate const& obj, int const& /*year*/){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double abs_eta = std::abs(obj.eta()); // Has to be the uncorrected one
    return (uncorr_pt>=selection_skim_fsr_pt && abs_eta<selection_skim_fsr_eta);
  }
  bool testSkimFSR_Iso(pat::PackedCandidate const& obj, int const& /*year*/, double const& fsrIso){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double relFSRIso = fsrIso/uncorr_pt;
    return (relFSRIso<selection_skim_fsr_reliso);
  }
  bool testSkimFSR_MinDeltaR(pat::PackedCandidate const& obj, int const& /*year*/, double const& mindr){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double min_rel_mindr = mindr / std::pow(uncorr_pt, 2);
    return (min_rel_mindr<selection_skim_fsr_min_reldr);
  }

}
