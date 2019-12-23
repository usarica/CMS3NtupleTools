#ifndef ELECTRONSELECTIONHELPERS_H
#define ELECTRONSELECTIONHELPERS_H

#include "ElectronObject.h"


namespace ElectronSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

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
  enum ElectronId{
    kCutBasedId_Fall17V1,
    kCutBasedId_Fall17V2,
    kMVAId_Fall17V2_NoIso,
    kMVAId_Fall17V2_Iso,
    kMVAId_HZZRun2Legacy_Iso
  };
  enum ElectronIso{
    kMiniIsoDR0p3,
    kPFIsoDR0p3,
    kPFIsoDR0p4,
    kMVAIso_Fall17V2_Iso,
    kMVAIso_HZZRun2Legacy_Iso
  };

  // Kinematic pT thresholds
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 5.;
  constexpr float ptThr_skim_medium = 5.;
  constexpr float ptThr_skim_tight = 5.;

  // Kinematic eta thresholds
  // Last ECAL crystal in barrel is at |eta|=1.4442
  constexpr float etaThr_gen = 2.5;
  constexpr float etaThr_skim_veto = 2.5;
  constexpr float etaThr_skim_loose = 2.5;
  constexpr float etaThr_skim_medium = 2.5;
  constexpr float etaThr_skim_tight = 2.5;

  // Isolation thresholds
  constexpr float isoThr_veto = 0.1;
  constexpr float isoThr_loose = 0.1;
  constexpr float isoThr_medium = 0.1;
  constexpr float isoThr_tight = 0.1;

  constexpr SelectionBits bit_preselection_id = kTightId;
  constexpr SelectionBits bit_preselection_iso = kTightIso;
  constexpr SelectionBits bit_preselection_kin = kTightKin;
  constexpr ElectronId idType_preselection = kMVAId_Fall17V2_NoIso;
  constexpr ElectronIso isoType_preselection = kPFIsoDR0p3;

  float absMiniIso_DR0p3(ElectronObject const& part);
  float relMiniIso_DR0p3(ElectronObject const& part);

  float absPFIso_DR0p3(ElectronObject const& part);
  float relPFIso_DR0p3(ElectronObject const& part);

  float absPFIso_DR0p4(ElectronObject const& part);
  float relPFIso_DR0p4(ElectronObject const& part);

  float computeIso(ElectronObject const& part);

  bool testPtEtaGen(ElectronObject const& part);

  bool testVetoId(ElectronObject const& part);
  bool testVetoIso(ElectronObject const& part);
  bool testVetoKin(ElectronObject const& part);

  bool testLooseId(ElectronObject const& part);
  bool testLooseIso(ElectronObject const& part);
  bool testLooseKin(ElectronObject const& part);

  bool testMediumId(ElectronObject const& part);
  bool testMediumIso(ElectronObject const& part);
  bool testMediumKin(ElectronObject const& part);

  bool testTightId(ElectronObject const& part);
  bool testTightIso(ElectronObject const& part);
  bool testTightKin(ElectronObject const& part);

  bool testPreselection(ElectronObject const& part);

  void setSelectionBits(ElectronObject& part);

}


#endif
