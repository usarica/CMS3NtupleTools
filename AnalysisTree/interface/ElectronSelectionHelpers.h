#ifndef ELECTRONSELECTIONHELPERS_H
#define ELECTRONSELECTIONHELPERS_H

#include "ElectronObject.h"


namespace ElectronSelectionHelpers{
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 5.;
  constexpr float ptThr_skim_medium = 5.;
  constexpr float ptThr_skim_tight = 5.;

  // Last ECAL crystal in barrel is at |eta|=1.4442
  constexpr float etaThr_gen = 2.5;
  constexpr float etaThr_skim_veto = 2.5;
  constexpr float etaThr_skim_loose = 2.5;
  constexpr float etaThr_skim_medium = 2.5;
  constexpr float etaThr_skim_tight = 2.5;

  constexpr float isoThr_veto = 0.1;
  constexpr float isoThr_loose = 0.1;
  constexpr float isoThr_medium = 0.1;
  constexpr float isoThr_tight = 0.1;

  enum SelectionBits{
    kGenPtEta,
    kVetoID,
    kVetoIDIso,
    kLooseID,
    kLooseIDIso,
    kMediumID,
    kMediumIDIso,
    kTightID,
    kTightIDIso,

    kSkimPtEta,
    kPreselection,

    nSelectionBits
  };
  const SelectionBits bit_preselection_idiso = kTightID;
  const SelectionBits bit_preselection_idisoreco = kTightIDIso;

  float absMiniIso_DR0p3(ElectronObject const& part);
  float relMiniIso_DR0p3(ElectronObject const& part);

  float absPFIso_DR0p3(ElectronObject const& part);
  float relPFIso_DR0p3(ElectronObject const& part);

  float absPFIso_DR0p4(ElectronObject const& part);
  float relPFIso_DR0p4(ElectronObject const& part);

  bool testPtEtaGen(ElectronObject const& part);

  bool testVetoId(ElectronObject const& part);
  bool testVetoIdIso(ElectronObject const& part);

  bool testLooseId(ElectronObject const& part);
  bool testLooseIdIso(ElectronObject const& part);

  bool testMediumId(ElectronObject const& part);
  bool testMediumIdIso(ElectronObject const& part);

  bool testTightId(ElectronObject const& part);
  bool testTightIdIso(ElectronObject const& part);

  bool testPtEtaSkim(ElectronObject const& part);
  bool testPreselection(ElectronObject const& part);

  void setSelectionBits(ElectronObject& part);

}


#endif
