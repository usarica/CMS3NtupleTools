#ifndef PHOTONSELECTIONHELPERS_H
#define PHOTONSELECTIONHELPERS_H

#include "PhotonObject.h"


namespace PhotonSelectionHelpers{
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_veto = 25.;
  constexpr float ptThr_skim_loose = 25.;
  constexpr float ptThr_skim_medium = 50.;
  constexpr float ptThr_skim_tight = 50.;

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
