#ifndef CMS3_PHOTONSELECTIONHELPERS_H
#define CMS3_PHOTONSELECTIONHELPERS_H

#include <cmath>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include "CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h"


namespace PhotonSelectionHelpers{
  enum EffectiveAreaType{
    PhotonEA_ch,
    PhotonEA_nh,
    PhotonEA_em
  };

  // Skim selection
  constexpr double selection_skim_pt = 5.;
  constexpr double selection_skim_eta = 5.0;
  constexpr double selection_iso_deltaR = 0.3; // This is a constant since the only isolation available uses dR<0.3

  // FSR skim
  constexpr double selection_match_fsr_deltaR = 0.5;
  constexpr double selection_skim_fsr_pt = 2.;
  constexpr double selection_skim_fsr_min_reldr = 0.012; // min. DR/pt**2 threshold
  constexpr double selection_skim_fsr_reliso = 1.8; // fsrIso / pT threshold

  float photonEffArea(pat::Photon const& obj, int const& year, PhotonSelectionHelpers::EffectiveAreaType const& eatype); // For PF and mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/photons_cff.py EAFile_MiniIso entries

  template<typename PFCandIterable> float photonFSRIso(pat::Photon const& obj, int const& year, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end);

  float photonPFIsoChargedHadron(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoNeutralHadron(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoEM(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoComb(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta

  bool testSkimPhoton(pat::Photon const& obj, int const& year);

}

template<typename PFCandIterable> float PhotonSelectionHelpers::photonFSRIso(pat::Photon const& obj, int const& /*year*/, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end){
  constexpr double cut_deltaR = selection_iso_deltaR;

  constexpr double cut_deltaRself_ch = 0.0001;
  constexpr double cut_pt_ch = 0.2;
  double sum_ch=0;

  constexpr double cut_deltaRself_ne = 0.01;
  constexpr double cut_pt_ne = 0.5;
  double sum_nh=0;

  for (PFCandIterable it_pfcands = pfcands_begin; it_pfcands!=pfcands_end; it_pfcands++){
    pat::PackedCandidate const* pfcand;
    CMS3ObjectHelpers::getObjectPointer(it_pfcands, pfcand);
    if (!pfcand) continue;

    double dr = reco::deltaR(obj.p4(), pfcand->p4());
    if (dr>=cut_deltaR) continue;

    int id = std::abs(pfcand->pdgId());
    int charge = pfcand->charge();
    double pt = pfcand->pt();
    if (charge!=0){
      // Charged hadrons
      if (dr>cut_deltaRself_ch && pt>cut_pt_ch && id==211) sum_ch += pt;
    }
    else{
      // Neutral particles
      if (dr>cut_deltaRself_ne && pt>cut_pt_ne && (id==22|| id==130)) sum_nh += pt;
    }
  }

  return (sum_ch + sum_nh);
}


#endif
