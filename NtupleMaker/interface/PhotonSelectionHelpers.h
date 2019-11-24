#ifndef CMS3_PHOTONSELECTIONHELPERS_H
#define CMS3_PHOTONSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Photon.h>


namespace PhotonSelectionHelpers{
  enum EffectiveAreaType{
    PhotonEA_ch,
    PhotonEA_nh,
    PhotonEA_em
  };

  // Skim selection
  constexpr double selection_skim_pt = 70.;
  constexpr double selection_skim_eta = 5.0;

  float photonEffArea(pat::Photon const& obj, int const& year, PhotonSelectionHelpers::EffectiveAreaType const& eatype); // For PF and mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/photons_cff.py EAFile_MiniIso entries

  float photonPFIsoCharged(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoComb(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta

  bool testSkimPhoton(pat::Photon const& obj, int const& year);

}


#endif
