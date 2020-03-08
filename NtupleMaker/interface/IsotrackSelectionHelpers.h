#ifndef CMS3_ISOTRACKSELECTIONHELPERS_H
#define CMS3_ISOTRACKSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>
#include <CMS3/NtupleMaker/interface/IsotrackInfo.h>


namespace IsotrackSelectionHelpers{
  enum IsolationType{
    PFIso03,
    MiniIso
  };

  // Skim selection
  constexpr double selection_skim_pt = 5.;
  constexpr double selection_skim_hadron_pt = 10.;
  constexpr double selection_skim_hadron_eta = 2.5;
  // No eta cut on leptons or general tracks

  float isotrackPFIsoComb(pat::IsolatedTrack const& obj, int const& year, IsotrackSelectionHelpers::IsolationType const& type, double const& fsr); // Absolute PF iso. value, uses delta beta correction instead of rho
  float isotrackMiniIsoComb(pat::IsolatedTrack const& obj, int const& year, double const& fsr); // Absolute mini. iso. value

  bool testSkimIsotrack(IsotrackInfo const& obj, int const& year);

}


#endif
