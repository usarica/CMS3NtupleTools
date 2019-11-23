#ifndef CMS3_AK8JETSELECTIONHELPERS_H
#define CMS3_AK8JETSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Jet.h>


namespace AK8JetSelectionHelpers{
  enum AK8JetType{
    AK8PFCHS,
    AK8PFPUPPI
  };

  // Skim selection
  constexpr double selection_skim_pt = 100.;
  constexpr double selection_skim_eta = 4.7;

  bool testSkimAK8Jet(pat::Jet const& obj, int const& year);
}


#endif
