#ifndef PHOTONSELECTIONHELPERS_H
#define PHOTONSELECTIONHELPERS_H

#include "PhotonObject.h"


namespace PhotonSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

    kConversionSafe,

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

    kPreselectionVeto,
    kPreselectionLoose,
    kPreselectionTight,

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
  constexpr float ptThr_skim_veto = 30.;
  constexpr float ptThr_skim_loose = 30.;
  constexpr float ptThr_skim_medium = 30.;
  constexpr float ptThr_skim_tight = 30.;

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

  constexpr PhotonId idType_preselection = kCutBasedId_Fall17V2;
  constexpr PhotonIso isoType_preselection = kPFIsoDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kMediumId;
  constexpr SelectionBits bit_preselectionVeto_iso = kLooseIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kLooseKin;

  constexpr SelectionBits bit_preselectionLoose_id = kMediumId;
  constexpr SelectionBits bit_preselectionLoose_iso = kLooseIso;
  constexpr SelectionBits bit_preselectionLoose_kin = kLooseKin;

  constexpr SelectionBits bit_preselectionTight_id = kMediumId;
  constexpr SelectionBits bit_preselectionTight_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionTight_kin = kMediumKin;

  constexpr SelectionBits bit_preselection_conversion = nSelectionBits; // kConversionSafe (enable in loose and tight preselection) or nSelectionBits (disable)

  float absPFIso_DR0p3(PhotonObject const& part);
  float relPFIso_DR0p3(PhotonObject const& part);

  void setSelectionBits(PhotonObject& part);

}


#endif
