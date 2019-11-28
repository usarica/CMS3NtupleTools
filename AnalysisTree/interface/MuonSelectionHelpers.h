#ifndef MUONSELECTIONHELPERS_H
#define MUONSELECTIONHELPERS_H

#include "MuonObject.h"


namespace MuonSelectionHelpers{
  constexpr float ptThr_gen = 7.;
  constexpr float ptThr_skim_veto = 7.;
  constexpr float ptThr_skim_loose = 7.;
  constexpr float ptThr_skim_medium = 7.;
  constexpr float ptThr_skim_tight = 7.;

  constexpr float etaThr_gen = 2.4;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_medium = 2.4;
  constexpr float etaThr_skim_tight = 2.4;

  constexpr float isoThr_veto = 0.1;
  constexpr float isoThr_loose = 0.1;
  constexpr float isoThr_medium = 0.1;
  constexpr float isoThr_tight = 0.1;

  enum SelectionBits{
    kGenPtEta,

    kValidMuonSystemTime,

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

  float absMiniIso_DR0p3(MuonObject const& part);
  float relMiniIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p3(MuonObject const& part);
  float relPFIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p4(MuonObject const& part);
  float relPFIso_DR0p4(MuonObject const& part);

  bool testPtEtaGen(MuonObject const& part);

  bool testMuonSystemTime(MuonObject const& part);

  bool testVetoId(MuonObject const& part);
  bool testVetoIdIso(MuonObject const& part);

  bool testLooseId(MuonObject const& part);
  bool testLooseIdIso(MuonObject const& part);

  bool testMediumId(MuonObject const& part);
  bool testMediumIdIso(MuonObject const& part);

  bool testTightId(MuonObject const& part);
  bool testTightIdIso(MuonObject const& part);

  bool testPtEtaSkim(MuonObject const& part);
  bool testPreselection(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

}


#endif
