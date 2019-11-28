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

  float absPFIso_DR0p3(PhotonObject const& part);
  float relPFIso_DR0p3(PhotonObject const& part);

  bool testPtEtaGen(PhotonObject const& part);

  bool testVetoId(PhotonObject const& part);
  bool testVetoIdIso(PhotonObject const& part);

  bool testLooseId(PhotonObject const& part);
  bool testLooseIdIso(PhotonObject const& part);

  bool testMediumId(PhotonObject const& part);
  bool testMediumIdIso(PhotonObject const& part);

  bool testTightId(PhotonObject const& part);
  bool testTightIdIso(PhotonObject const& part);

  bool testPtEtaSkim(PhotonObject const& part);
  bool testPreselection(PhotonObject const& part);

  void setSelectionBits(PhotonObject& part);

}


#endif
