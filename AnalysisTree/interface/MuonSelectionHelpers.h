#ifndef MUONSELECTIONHELPERS_H
#define MUONSELECTIONHELPERS_H

#include "MuonObject.h"


namespace MuonSelectionHelpers{
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

    kSoftId,
    kSoftIso,
    kSoftKin,

    kPreselection,

    nSelectionBits
  };
  enum MuonId{
    kCutBasedId_MuonPOG
  };
  enum MuonIso{
    kMiniIsoDR0p3,
    kPFIsoDR0p3,
    kPFIsoDR0p4
  };

  // Kinematic pT thresholds
  constexpr float ptThr_gen = 7.;
  constexpr float ptThr_skim_veto = 7.;
  constexpr float ptThr_skim_loose = 7.;
  constexpr float ptThr_skim_medium = 7.;
  constexpr float ptThr_skim_tight = 7.;
  constexpr float ptThr_skim_soft = 3.;

  // Kinematic eta thresholds
  constexpr float etaThr_gen = 2.4;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_medium = 2.4;
  constexpr float etaThr_skim_tight = 2.4;
  constexpr float etaThr_skim_soft = 2.4;

  // Isolation thresholds
  constexpr float isoThr_veto = 0.2;
  constexpr float isoThr_loose = 0.2;
  constexpr float isoThr_medium = 0.1;
  constexpr float isoThr_tight = 0.1;
  constexpr float isoThr_soft = 0.2;

  constexpr SelectionBits bit_preselection_id = kTightId;
  constexpr SelectionBits bit_preselection_iso = kTightIso;
  constexpr SelectionBits bit_preselection_kin = kTightKin;
  constexpr SelectionBits bit_preselection_time = kValidMuonSystemTime;
  constexpr MuonId idType_preselection = kCutBasedId_MuonPOG;
  constexpr MuonIso isoType_preselection = kPFIsoDR0p3;

  float absMiniIso_DR0p3(MuonObject const& part);
  float relMiniIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p3(MuonObject const& part);
  float relPFIso_DR0p3(MuonObject const& part);

  float absPFIso_DR0p4(MuonObject const& part);
  float relPFIso_DR0p4(MuonObject const& part);

  float computeIso(MuonObject const& part);

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

  bool testSoftId(MuonObject const& part);
  bool testSoftIso(MuonObject const& part);
  bool testSoftKin(MuonObject const& part);

  bool testPreselection(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

}


#endif
