#ifndef ELECTRONSELECTIONHELPERS_H
#define ELECTRONSELECTIONHELPERS_H

#include <vector>
#include <CMS3/Dictionaries/interface/ElectronTriggerCutEnums.h>
#include "ElectronObject.h"


namespace ElectronSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

    kConversionSafe,
    kInTimeSeed,

    kSpikeSafe,

    kPFElectronId,
    kPFElectronPreferable,
    kPFMETSafe,

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

    kProbeId,

    kFakeableBaseIso,
    kFakeableBase,
    kFakeable,

    kPreselectionVeto,
    kPreselectionLoose_NoIso,
    kPreselectionLoose,
    kPreselectionTight,

    nSelectionBits
  };
  enum ElectronId{
#if ELECTRONS_HAVE_FALL17V1_CUTBASED == 1
    kCutBasedId_Fall17V1,
#endif
    kCutBasedId_Fall17V2,
    kMVAId_Fall17V2_NoIso,
    kMVAId_Fall17V2_Iso,
    kMVAId_HZZRun2Legacy_Iso
  };
  enum ElectronIso{
    kMiniIso,
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
  // Gap region is between 1.4442 and 1.56, crossing is at 1.479. See ECALGeometrySpecifications.h.
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
  constexpr float isoThr_fakeable = 0.4;

  // Seed time threshold for cosmics and other stuff
  // In ns. A bit tight, but oh well...
  constexpr float seedTimeThr = 1.;

  // Spike cleanup
  constexpr float full5x5_sigmaIEtaIEtaThr = 0.001;
  constexpr float full5x5_sigmaIPhiIPhiThr = 0.001;

  // PF electron id delta R matching threshold
  // Hopefully all values are already ~0 in reco.
  constexpr float mindRThr_electron_pfelectron = 0.04;

  constexpr ElectronId idType_preselection = kMVAId_Fall17V2_NoIso;
  constexpr ElectronIso isoType_preselection = kPFIsoDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kMediumId;
  constexpr SelectionBits bit_preselectionVeto_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kMediumKin;

  constexpr SelectionBits bit_preselectionLoose_id = kMediumId;
  constexpr SelectionBits bit_preselectionLoose_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionLoose_kin = kMediumKin;

  constexpr SelectionBits bit_preselectionTight_id = kMediumId;
  constexpr SelectionBits bit_preselectionTight_iso = kMediumIso;
  constexpr SelectionBits bit_preselectionTight_kin = kMediumKin;

  float getIsolationDRmax(ElectronObject const& part);

  float absPFIso_DR0p3(ElectronObject const& part);
  float relPFIso_DR0p3(ElectronObject const& part);

  float absPFIso_DR0p4(ElectronObject const& part);
  float relPFIso_DR0p4(ElectronObject const& part);

  float absMiniIso(ElectronObject const& part);
  float relMiniIso(ElectronObject const& part);

  float computeIso(ElectronObject const& part);

  void setSelectionBits(ElectronObject& part);

  // For trigger emulation selection
  void clearTriggerEmulationBits();
  void setTriggerEmulationV1Bits(std::vector<ElectronTriggerCutEnums::ElectronTriggerCutTypes> const& bitlist);
  void setTriggerEmulationV2Bits(std::vector<ElectronTriggerCutEnums::ElectronTriggerCutTypes> const& bitlist);

  // User functions to disable or enable other selection features for 'loose' and 'tight' preselection
  void setApplyConversionSafety(bool flag);
  void setApplySeedTimeVeto(bool flag);
  void setApplySpikeVeto(bool flag);
  // Notice these two are separated!
  // Notice also that if PF id is to be applied, there should be an external check
  // for the kPFElectronPreferable flag as well when photon overlaps are present.
  void setApplyPFId(bool flag);
  void setApplyPFMETSafety(bool flag);
  // T&P and fake rate functions
  void setAllowProbeIdInLooseSelection(bool flag);
  void setAllowFakeableInLooseSelection(bool flag);

  // Get functions to read the state of the flags
  bool getApplyConversionSafety();
  bool getApplySeedTimeVeto();
  bool getApplySpikeVeto();
  bool getApplyPFId();
  bool getApplyPFMETSafety();
  bool getAllowProbeIdInLooseSelection();
  bool getAllowFakeableInLooseSelection();

}


#endif
