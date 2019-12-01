#ifndef AK8JETSELECTIONHELPERS_H
#define AK8JETSELECTIONHELPERS_H

#include "AK8JetObject.h"


namespace AK8JetSelectionHelpers{
  constexpr float ptThr_gen = 150.;
  constexpr float ptThr_skim_loose = 175.;
  constexpr float ptThr_skim_tight = 200.;

  constexpr float etaThr_gen = 5.2;
  constexpr float etaThr_skim_loose = 4.7;
  constexpr float etaThr_skim_tight = 4.7;

  enum SelectionBits{
    kGenPtEta,

    kLooseId,
    kLooseKin,
    kTightId,
    kTightKin,

    kPreselection,

    nSelectionBits
  };
  const SelectionBits bit_preselection_id = kTightId;
  const SelectionBits bit_preselection_kin = kTightKin;

  bool testPtEtaGen(AK8JetObject const& part);

  bool testLooseId(AK8JetObject const& part);
  bool testLooseKin(AK8JetObject const& part);

  bool testTightId(AK8JetObject const& part);
  bool testTightKin(AK8JetObject const& part);

  bool testPreselection(AK8JetObject const& part);

  void setSelectionBits(AK8JetObject& part);
}


#endif
