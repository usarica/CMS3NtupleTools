#ifndef CMS3_AK8JETSELECTIONHELPERS_H
#define CMS3_AK8JETSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Jet.h>


namespace AK8JetSelectionHelpers{
  enum AK8JetType{
    AK8PFCHS,
    AK8PFPUPPI
  };

  // Cone radius
  constexpr double ConeRadiusConstant = 0.8;

  // Skim selection
  constexpr double selection_skim_pt = 100.;
  constexpr double selection_skim_eta = 4.7;

  double getUncorrectedJetEnergy(pat::Jet const& obj);
  double getUncorrectedJetPt(pat::Jet const& obj);
  double getUncorrectedJetMass(pat::Jet const& obj);

  bool testLooseAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type);
  bool testTightAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type);
  bool testLeptonVetoAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type);
  // No need for PU jet id

  bool testSkimAK8Jet(pat::Jet const& obj, int const& year);

}


#endif
