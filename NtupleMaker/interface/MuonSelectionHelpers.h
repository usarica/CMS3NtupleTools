#ifndef CMS3_MUONSELECTIONHELPERS_H
#define CMS3_MUONSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Muon.h>


namespace MuonSelectionHelpers{
  enum IsolationType{
    PFIso03,
    PFIso04,
    MiniIso
  };

  // Skim selection
  constexpr double selection_skim_pt = 5.;
  constexpr double selection_skim_eta = 2.4;

  float muonEffArea(pat::Muon const& obj, int const& year); // For mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/muons_cff.py EAFile_MiniIso entries

  float muonPFIsoComb(pat::Muon const& obj, int const& year, MuonSelectionHelpers::IsolationType const& type, double const& fsr); // Absolute PF iso. value, uses delta beta correction instead of rho
  float muonMiniIsoComb(pat::Muon const& obj, int const& year, double const& rho, double const& fsr); // Absolute mini. iso. value

  bool testSkimMuon(pat::Muon const& obj, int const& year);

}


#endif
