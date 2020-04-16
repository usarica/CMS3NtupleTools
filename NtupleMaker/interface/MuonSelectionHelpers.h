#ifndef CMS3_MUONSELECTIONHELPERS_H
#define CMS3_MUONSELECTIONHELPERS_H

#include <cassert>
#include <DataFormats/PatCandidates/interface/Muon.h>


namespace MuonSelectionHelpers{
  enum IsolationType{
    PFIso03,
    PFIso04,
    MiniIso
  };

  // Skim selection
  constexpr double selection_skim_pt = 3.;
  constexpr double selection_skim_eta = 2.4;

  constexpr double tightCharge_pt_err_rel_thr = 0.2;


  template<typename T> constexpr unsigned int getPOGSelectorBitPosition(T const&);

  float muonEffArea(pat::Muon const& obj, int const& year); // For mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/muons_cff.py EAFile_MiniIso entries

  float muonPFIsoComb(pat::Muon const& obj, int const& year, MuonSelectionHelpers::IsolationType const& type, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr, double* sum_neutral_EAcorr_nofsr); // Absolute PF iso. value, uses delta beta correction instead of EA*rho, but sum_neutral_EAcorr_nofsr return the EA-corrected value for optional use.
  float muonMiniIsoComb(pat::Muon const& obj, int const& year, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr); // Absolute mini. iso. value

  bool testLooseTriggerId(pat::Muon const& obj, int const& year); // Test loose trigger id
  bool testTightCharge(pat::Muon const& obj, int const& year); // Test error on charge determination via track pT error
  bool testMuonTiming(pat::Muon const& obj, int const& year); // Test muon timing from RPC and combined measurements

  bool testSkimMuon(pat::Muon const& obj, int const& year);

}

template<typename T> constexpr unsigned int MuonSelectionHelpers::getPOGSelectorBitPosition(T const& mask){
  unsigned int i = 0;
  while (!((mask >> i) & 1)){
    i++;
    if (i>64){
      assert(0);
      break;
    }
  }
  return i;
}


#endif
