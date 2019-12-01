#ifndef AK4JETSELECTIONHELPERS_H
#define AK4JETSELECTIONHELPERS_H

#include "AK4JetObject.h"


namespace AK4JetSelectionHelpers{
  constexpr float ptThr_gen = 15.;
  constexpr float ptThr_skim_loose = 30.;
  constexpr float ptThr_skim_tight = 30.;

  constexpr float etaThr_gen = 5.2;
  constexpr float etaThr_skim_loose = 4.7;
  constexpr float etaThr_skim_tight = 4.7;

  enum SelectionBits{
    kGenPtEta,

    kPUJetId,

    kLooseId,
    kLooseKin,
    kTightId,
    kTightKin,

    kPreselection,

    nSelectionBits
  };
  const SelectionBits bit_preselection_id = kTightId;
  const SelectionBits bit_preselection_kin = kTightKin;

  bool testPtEtaGen(AK4JetObject const& part);

  bool testPUJetId(AK4JetObject const& part);

  bool testLooseId(AK4JetObject const& part);
  bool testLooseKin(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightKin(AK4JetObject const& part);

  bool testPreselection(AK4JetObject const& part);

  void setSelectionBits(AK4JetObject& part);
}


#endif
