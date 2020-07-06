#ifndef CMS3_AK4JETSELECTIONHELPERS_H
#define CMS3_AK4JETSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Jet.h>


namespace AK4JetSelectionHelpers{
  enum AK4JetType{
    AK4PFCHS,
    AK4PFPUPPI
  };

  // Skim selection
  constexpr double selection_skim_pt = 20.;
  constexpr double selection_skim_eta = 4.7;

  double getUncorrectedJetEnergy(pat::Jet const& obj);
  double getUncorrectedJetPt(pat::Jet const& obj);

  bool testSkimAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  bool testLooseAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  bool testTightAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  bool testLeptonVetoAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  // Bad muon id requires accessing MET, so it should not be done here

}


#endif
