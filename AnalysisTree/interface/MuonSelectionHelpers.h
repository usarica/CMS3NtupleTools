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

    kProbeId, // TnP
    kProbeSTAId, // TnP STA

    kFakeableBaseIso,
    kFakeableBase,
    kFakeable,

    kPreselectionVeto,
    kPreselectionLoose,
    kPreselectionTight,

    nSelectionBits
  };
  enum MuonId{
    kCutBasedId_MuonPOG
  };
  enum MuonIso{
    kMiniIso,
    kPFIsoDR0p3,
    kPFIsoDR0p3_EACorrected,
    kPFIsoDR0p4,
    kPFIsoDR0p4_EACorrected
  };

  // Kinematic pT thresholds
  constexpr float ptThr_gen = 3.;
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 5.;
  constexpr float ptThr_skim_medium = 5.;
  constexpr float ptThr_skim_tight = 5.;
  constexpr float ptThr_skim_soft = 3.;

  // Kinematic eta thresholds
  // Barrel region ends at |eta|=1.2
  constexpr float etaThr_gen = 2.4;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_medium = 2.4;
  constexpr float etaThr_skim_tight = 2.4;
  constexpr float etaThr_skim_soft = 2.4;

  // Isolation thresholds
  constexpr float isoThr_veto = 0.15;
  constexpr float isoThr_loose = 0.15;
  constexpr float isoThr_medium = 0.15;
  constexpr float isoThr_tight = 0.15;
  constexpr float isoThr_soft = 0.15;
  constexpr float isoThr_fakeable = 0.4;

  // Determine how the physics objects are vetoed/cleaned/accepted
  constexpr MuonId idType_preselection = kCutBasedId_MuonPOG;
  constexpr MuonIso isoType_preselection = kPFIsoDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kMediumId;
  constexpr SelectionBits bit_preselectionVeto_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kMediumKin;

  constexpr SelectionBits bit_preselectionLoose_id = kMediumId;
  constexpr SelectionBits bit_preselectionLoose_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionLoose_kin = kMediumKin;

  constexpr SelectionBits bit_preselectionTight_id = kMediumId;
  constexpr SelectionBits bit_preselectionTight_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionTight_kin = kMediumKin;

  constexpr SelectionBits bit_preselection_time = nSelectionBits; // kValidMuonSystemTime (enable in loose and tight preselection) or nSelectionBits (disable)

  float getIsolationDRmax(MuonObject const& part);

  float absPFIso_DR0p3(MuonObject const& part);
  float relPFIso_DR0p3(MuonObject const& part);

  float absPFIso_EACorr_DR0p3(MuonObject const& part);
  float relPFIso_EACorr_DR0p3(MuonObject const& part);

  float absPFIso_DR0p4(MuonObject const& part);
  float relPFIso_DR0p4(MuonObject const& part);

  float absPFIso_EACorr_DR0p4(MuonObject const& part);
  float relPFIso_EACorr_DR0p4(MuonObject const& part);

  float absMiniIso(MuonObject const& part);
  float relMiniIso(MuonObject const& part);

  float computeIso(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

  void doRequireTrackerIsolationInFakeable(float const& isothr);
}


#endif
