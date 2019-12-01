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

    kVetoId,
    kVetoIso,
    kVetoKin,
    kLooseId,
    kLooseIso,
    kLooseKin,
    kMediumId,
    kMediumIso,
    kMediumKin,
    kTightId,
    kTightIso,
    kTightKin,

    kPreselection,

    nSelectionBits
  };
  const SelectionBits bit_preselection_id = kTightId;
  const SelectionBits bit_preselection_iso = kTightIso;
  const SelectionBits bit_preselection_kin = kTightKin;
  const SelectionBits bit_preselection_time = kValidMuonSystemTime;

  float absMiniIso_DR0p3(MuonObject const& part);
  float relMiniIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p3(MuonObject const& part);
  float relPFIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p4(MuonObject const& part);
  float relPFIso_DR0p4(MuonObject const& part);

  bool testPtEtaGen(MuonObject const& part);

  bool testMuonSystemTime(MuonObject const& part);

  bool testVetoId(MuonObject const& part);
  bool testVetoIso(MuonObject const& part);
  bool testVetoKin(MuonObject const& part);

  bool testLooseId(MuonObject const& part);
  bool testLooseIso(MuonObject const& part);
  bool testLooseKin(MuonObject const& part);

  bool testMediumId(MuonObject const& part);
  bool testMediumIso(MuonObject const& part);
  bool testMediumKin(MuonObject const& part);

  bool testTightId(MuonObject const& part);
  bool testTightIso(MuonObject const& part);
  bool testTightKin(MuonObject const& part);

  bool testPreselection(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

}


#endif
