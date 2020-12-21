#ifndef CMS3_AK4JETSELECTIONHELPERS_H
#define CMS3_AK4JETSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Jet.h>


namespace AK4JetSelectionHelpers{
  enum AK4JetType{
    AK4PFCHS,
    AK4PFPUPPI
  };

  // Cone radius
  constexpr double ConeRadiusConstant = 0.4;

  // MET pT threshold
  constexpr double selection_METJERC_pt = 15.; // Applied on corrected no-mu p4

  // Skim selection
  constexpr double selection_skim_pt = 20.;
  constexpr double selection_skim_eta = 4.7;

  double getUncorrectedJetEnergy(pat::Jet const& obj);
  double getUncorrectedJetPt(pat::Jet const& obj);
  double getUncorrectedJetMass(pat::Jet const& obj);

  bool testLooseAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  bool testTightAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);
  bool testLeptonVetoAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);

  bool testAK4JetMETSafety(pat::Jet const& obj);

  bool testAK4JetMETFixSafety_NoPt(double const& eta, int const& year);
  bool testAK4JetMETFixSafety_NoPt(pat::Jet const& obj, int const& year);
  bool testAK4JetMETFixSafety(double const& uncorrected_pt, double const& eta, int const& year);
  bool testAK4JetMETFixSafety(pat::Jet const& obj, int const& year, bool isPFJetMakerOutput);

  bool testSkimAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type);

}


#endif
