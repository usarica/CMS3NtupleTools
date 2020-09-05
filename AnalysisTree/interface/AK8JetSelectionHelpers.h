#ifndef AK8JETSELECTIONHELPERS_H
#define AK8JETSELECTIONHELPERS_H

#include "AK8JetObject.h"


namespace AK8JetSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

    kLooseId,
    kLooseKin,
    kTightId,
    kTightKin,

    kPreselectionLoose,
    kPreselectionTight,

    nSelectionBits
  };

  constexpr float ptThr_gen = 150.;
  constexpr float ptThr_skim_HEMcheck = 30.; // Literally any ak8 jet
  constexpr float ptThr_skim_loose = 175.;
  constexpr float ptThr_skim_tight = 200.;

  constexpr float etaThr_gen = 5.2;
  constexpr float etaThr_skim_loose = 4.7;
  constexpr float etaThr_skim_tight = 4.7;

  const SelectionBits bit_preselectionLoose_id = kTightId;
  const SelectionBits bit_preselectionLoose_kin = kTightKin;

  const SelectionBits bit_preselectionTight_id = kTightId;
  const SelectionBits bit_preselectionTight_kin = kTightKin;

  void setSelectionBits(AK8JetObject& part);
}


#endif
