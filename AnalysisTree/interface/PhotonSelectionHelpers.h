#ifndef PHOTONSELECTIONHELPERS_H
#define PHOTONSELECTIONHELPERS_H

#include "PhotonObject.h"


namespace PhotonSelectionHelpers{
  enum SelectionBits{
    kGenPtEta,

    kConversionSafe,
    kInTimeSeed,
    kBeamHaloSafe,

    kSpikeSafe,

    kPFPhotonId,
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

    kPreselectionVeto,
    kPreselectionLoose,
    kPreselectionTight,

    kSFTampon,

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
  constexpr float ptThr_skim_veto = 20.;
  constexpr float ptThr_skim_loose = 20.;
  constexpr float ptThr_skim_medium = 20.;
  constexpr float ptThr_skim_tight = 20.;

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

  // Seed time threshold for cosmics and other stuff
  // In ns. Beware that 2016, 2017 vs 2018 have different distributions.
  constexpr float seedTimeThr_DnUp = 2.;
  constexpr float seedTimeThr_Up_2018 = 1.;

  // MIP total energy threshold for beam halo safety (meaningful for barrel photons only, endcap photons have MIPTotalEnergy set to 0)
  constexpr float mipTotalEnergyThr = 4.9;

  // Spike cleanup
  constexpr float full5x5_sigmaIEtaIEtaThr = 0.001;
  constexpr float full5x5_sigmaIPhiIPhiThr = 0.001;

  // PF photon id delta R matching threshold
  constexpr float mindRThr_photon_pfphoton = 0.04;

  constexpr PhotonId idType_preselection = kCutBasedId_Fall17V2;
  constexpr PhotonIso isoType_preselection = kPFIsoDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kTightId;
  constexpr SelectionBits bit_preselectionVeto_iso = kTightIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kTightKin;

  constexpr SelectionBits bit_preselectionLoose_id = kTightId;
  constexpr SelectionBits bit_preselectionLoose_iso = kTightIso;
  constexpr SelectionBits bit_preselectionLoose_kin = kTightKin;

  constexpr SelectionBits bit_preselectionTight_id = kTightId;
  constexpr SelectionBits bit_preselectionTight_iso = kTightIso;
  constexpr SelectionBits bit_preselectionTight_kin = kTightKin;

  constexpr SelectionBits bit_SFTampon_id = kMediumId;
  constexpr SelectionBits bit_SFTampon_iso = kMediumIso;
  constexpr SelectionBits bit_SFTampon_kin = bit_preselectionTight_kin;

  float getIsolationDRmax(PhotonObject const& part);

  float absPFIso_DR0p3(PhotonObject const& part);
  float relPFIso_DR0p3(PhotonObject const& part);

  void setSelectionBits(PhotonObject& part);

  // User functions to disable or enable selection features for 'loose' and 'tight' preselection
  void setApplyConversionSafety(bool flag);
  void setApplySeedTimeVeto(bool flag);
  void setApplyBeamHaloVeto(bool flag);
  void setApplySpikeVeto(bool flag);
  // Notice these two are separated!
  // Notice also that if PF id is to be applied, there should be an external check
  // for the ElectronSelectionHelpers::kPFElectronPreferable flag as well for any overlapping electron.
  void setApplyPFId(bool flag);
  void setApplyPFMETSafety(bool flag);

  // Get functions to read the state of the flags
  bool getApplyConversionSafety();
  bool getApplySeedTimeVeto();
  bool getApplyBeamHaloVeto();
  bool getApplySpikeVeto();
  bool getApplyPFId();
  bool getApplyPFMETSafety();

}


#endif
