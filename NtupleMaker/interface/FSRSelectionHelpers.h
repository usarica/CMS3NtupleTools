#ifndef CMS3_FSRSELECTIONHELPERS_H
#define CMS3_FSRSELECTIONHELPERS_H

#include <cmath>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include "CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h"


namespace FSRSelectionHelpers{
  // Photon iso dR value, fixed in miniAOD
  constexpr double selection_iso_deltaR = PhotonSelectionHelpers::selection_iso_deltaR;

  // FSR skim
  constexpr double selection_match_fsr_deltaR = 0.5;
  constexpr double selection_skim_fsr_pt = 2.;
  constexpr double selection_skim_fsr_eta = 2.4;
  constexpr double selection_skim_fsr_min_reldr = 0.012; // min. DR/pt**2 threshold
  constexpr double selection_skim_fsr_reliso = 1.8; // fsrIso / pT threshold

  template<typename PFCandIterable> float fsrIso(pat::PackedCandidate const& obj, int const& year, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end);

  // Test if the supercluster veto is to be applied
  bool testSCVeto(pat::PackedCandidate const* pfcand, pat::Electron const* electron);
  bool testSCVeto(pat::PackedCandidate const* pfcand, pat::Photon const* photon);

  // FSR preselection comparator
  template<typename T> bool fsrPreselectionComparator(T const* obj, pat::PackedCandidate const* pfcand);
  template<> bool fsrPreselectionComparator(pat::Electron const* obj, pat::PackedCandidate const* pfcand);

  bool testSkimFSR_PtEta(pat::PackedCandidate const& obj, int const& year);
  bool testSkimFSR_Iso(pat::PackedCandidate const& obj, int const& year, double const& fsrIso);
  bool testSkimFSR_MinDeltaR(pat::PackedCandidate const& obj, int const& year, double const& mindr);

}

template<typename T> bool FSRSelectionHelpers::fsrPreselectionComparator(T const* obj, pat::PackedCandidate const* pfcand){
  if (pfcand->pdgId()!=22 || pfcand->pt()<FSRSelectionHelpers::selection_skim_fsr_pt || std::abs(pfcand->eta())>=FSRSelectionHelpers::selection_skim_fsr_eta) return false;
  return (reco::deltaR(pfcand->p4(), obj->p4())<FSRSelectionHelpers::selection_match_fsr_deltaR);
}
template<> bool FSRSelectionHelpers::fsrPreselectionComparator(pat::Electron const* obj, pat::PackedCandidate const* pfcand){
  if (pfcand->pdgId()!=22 || pfcand->pt()<FSRSelectionHelpers::selection_skim_fsr_pt || std::abs(pfcand->eta())>=FSRSelectionHelpers::selection_skim_fsr_eta) return false;
  return (!testSCVeto(pfcand, obj) && reco::deltaR(pfcand->p4(), obj->p4())<FSRSelectionHelpers::selection_match_fsr_deltaR);
}

template<typename PFCandIterable> float FSRSelectionHelpers::fsrIso(pat::PackedCandidate const& obj, int const& /*year*/, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end){
  constexpr double cut_deltaR = selection_iso_deltaR;

  constexpr double cut_deltaRself_ch = 0.0001;
  constexpr double cut_pt_ch = 0.2;
  double sum_ch = 0;

  constexpr double cut_deltaRself_ne = 0.01;
  constexpr double cut_pt_ne = 0.5;
  double sum_ne = 0;

  for (PFCandIterable it_pfcands = pfcands_begin; it_pfcands!=pfcands_end; it_pfcands++){
    pat::PackedCandidate const* pfcand;
    CMS3ObjectHelpers::getObjectPointer(it_pfcands, pfcand);
    if (!pfcand) continue;
    if (&obj==pfcand) continue; // Obviously don't include self

    double dr = reco::deltaR(obj.p4(), pfcand->p4());
    if (dr>=cut_deltaR) continue;

    unsigned int abs_id = std::abs(pfcand->pdgId());
    int charge = pfcand->charge();
    double pt = pfcand->pt();
    if (charge!=0){
      // Charged hadrons
      if (abs_id==211 && dr>cut_deltaRself_ch && pt>cut_pt_ch) sum_ch += pt;
    }
    else{
      // Neutral particles
      if ((abs_id==22 || abs_id==130) && dr>cut_deltaRself_ne && pt>cut_pt_ne) sum_ne += pt;
    }
  }

  return (sum_ch + sum_ne);
}


#endif
