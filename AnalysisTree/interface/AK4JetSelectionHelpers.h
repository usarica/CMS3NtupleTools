#ifndef AK4JETSELECTIONHELPERS_H
#define AK4JETSELECTIONHELPERS_H

#include "AK4JetObject.h"


namespace AK4JetSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

    kLoosePUJetId,
    kMediumPUJetId,
    kTightPUJetId,

    kLoosePUJetId_default,
    kMediumPUJetId_default,
    kTightPUJetId_default,

    kTightLeptonVetoId,

    kNoisyJet,

    kLooseId,
    kLooseKin,

    kTightId,
    kTightKin,

    kBtaggableKin,
    kBtaggable_NoPUJetId,
    kBtaggable,

    // These two flags are distinct. Do not apply a direct AND of them.
    kMETJERCSafe_Base,
    kMETFixRecoverySuitable,

    kPreselectionLoose,
    kPreselectionTight,

    nSelectionBits
  };

  constexpr float ptThr_gen = 15.;
  constexpr float ptThr_METJERCSafety = 15.; // This is not exactly a cut on the pt of the jet itself, but on the muon-subtracted version.
  constexpr float ptThr_btag = 20.;
  constexpr float ptThr_PUId = 50.; // Upper bound
  constexpr float ptThr_skim_loose = ptThr_btag;
  constexpr float ptThr_skim_tight = 30.;

  constexpr float etaThr_gen = 5.2;
  constexpr float etaThr_skim_loose = 4.7;
  constexpr float etaThr_skim_tight = 4.7;
  constexpr float etaThr_PUId = 5.;

  // Preselection loose = AND of bits below
  const SelectionBits bit_preselectionLoose_id = kTightId;
  const SelectionBits bit_preselectionLoose_kin = kTightKin;

  // Preselection tight = AND of bits below && PU jet id
  const SelectionBits bit_preselectionTight_id = kTightId;
  const SelectionBits bit_preselectionTight_kin = kTightKin;

  // Control additional flags to preselection tight id
  void setPUIdWP(SelectionBits flag);
  void setApplyTightLeptonVetoIdToJets(bool flag);

  bool const& getApplyTightLeptonVetoIdToJetsFlag();

  void setSelectionBits(AK4JetObject& part, bool resetIDs = true, bool resetKinematics = true);
}


#endif
