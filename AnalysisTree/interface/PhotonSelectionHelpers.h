#ifndef PHOTONSELECTIONHELPERS_H
#define PHOTONSELECTIONHELPERS_H

#include "PhotonObject.h"


namespace PhotonSelectionHelpers{
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
  enum PhotonId{
    kCutBasedId_Fall17V2,
    kMVAId_Fall17V2
  };
  enum PhotonIso{
    kPFIsoDR0p3
  };

  // Kinematic pT thresholds
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_veto = 25.;
  constexpr float ptThr_skim_loose = 25.;
  constexpr float ptThr_skim_medium = 50.;
  constexpr float ptThr_skim_tight = 50.;

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
  constexpr PhotonId idType_preselection = kCutBasedId_Fall17V2; // FIXME: Do not use MVA id in photons. The conversion veto and pixel seed falgs are not yet recorded in the trees.
  constexpr PhotonIso isoType_preselection = kPFIsoDR0p3;

  float absPFIso_DR0p3(PhotonObject const& part);
  float relPFIso_DR0p3(PhotonObject const& part);

  bool testPtEtaGen(PhotonObject const& part);

  bool testVetoId(PhotonObject const& part);
  bool testVetoIso(PhotonObject const& part);
  bool testVetoKin(PhotonObject const& part);

  bool testLooseId(PhotonObject const& part);
  bool testLooseIso(PhotonObject const& part);
  bool testLooseKin(PhotonObject const& part);

  bool testMediumId(PhotonObject const& part);
  bool testMediumIso(PhotonObject const& part);
  bool testMediumKin(PhotonObject const& part);

  bool testTightId(PhotonObject const& part);
  bool testTightIso(PhotonObject const& part);
  bool testTightKin(PhotonObject const& part);

  bool testPreselection(PhotonObject const& part);

  void setSelectionBits(PhotonObject& part);

}


#endif
