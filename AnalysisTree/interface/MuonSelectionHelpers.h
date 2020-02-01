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

    kPreselectionVeto,
    kPreselectionLoose,
    kPreselectionAccept,

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
  constexpr float isoThr_medium = 0.15;
  constexpr float isoThr_tight = 0.15;
  constexpr float isoThr_soft = 0.2;

  // Determine how the physics objects are vetoed/cleaned/accepted
  constexpr MuonId idType_preselection = kCutBasedId_MuonPOG;
  constexpr MuonIso isoType_preselection = kPFIsoDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kMediumId;
  constexpr SelectionBits bit_preselectionVeto_iso = kLooseIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kLooseKin;

  constexpr SelectionBits bit_preselectionLoose_id = kMediumId;
  constexpr SelectionBits bit_preselectionLoose_iso = kLooseIso;
  constexpr SelectionBits bit_preselectionLoose_kin = kLooseKin;

  constexpr SelectionBits bit_preselectionAccept_id = kMediumId;
  constexpr SelectionBits bit_preselectionAccept_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionAccept_kin = kMediumKin;

  constexpr SelectionBits bit_preselection_time = kValidMuonSystemTime;

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

  bool testPreselectionVeto(MuonObject const& part);
  bool testPreselectionLoose(MuonObject const& part);
  bool testPreselectionAccept(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

}


#endif
