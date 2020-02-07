#ifndef ISOTRACKSELECTIONHELPERS_H
#define ISOTRACKSELECTIONHELPERS_H

#include "IsotrackObject.h"


namespace IsotrackSelectionHelpers{
  enum SelectionBits{
    kVetoId,
    kVetoIso,
    kVetoKin,

    kPreselectionVeto,

    nSelectionBits
  };
  enum IsotrackIso{
    kMiniIsoCharged,
    kMiniIsoComb,
    kPFIsoChargedDR0p3,
    kPFIsoCombDR0p3
  };

  // Kinematic pT thresholds
  constexpr float ptThr_lepton_veto = 5.;
  constexpr float ptThr_hadron_veto = 10.;

  // Kinematic eta thresholds
  constexpr float etaThr_lepton_veto = 2.4;
  constexpr float etaThr_hadron_veto = 2.4;

  // Isolation thresholds
  constexpr float relIsoThr_lepton_veto = 0.2;
  constexpr float absIsoThr_lepton_veto = 5.;
  constexpr float relIsoThr_hadron_veto = 0.2;
  constexpr float absIsoThr_hadron_veto = 5.;

  // Determine how the physics objects are vetoed/cleaned/accepted
  constexpr IsotrackIso isoType_preselection = kPFIsoChargedDR0p3;

  constexpr SelectionBits bit_preselectionVeto_id = kVetoId;
  constexpr SelectionBits bit_preselectionVeto_iso = kVetoIso;
  constexpr SelectionBits bit_preselectionVeto_kin = kVetoKin;

  float getIsolationDRmax(IsotrackObject const& part);

  float absMiniIsoCharged(IsotrackObject const& part);
  float relMiniIsoCharged(IsotrackObject const& part);

  float absPFIsoCharged_DR0p3(IsotrackObject const& part);
  float relPFIsoCharged_DR0p3(IsotrackObject const& part);

  float absMiniIsoComb(IsotrackObject const& part);
  float relMiniIsoComb(IsotrackObject const& part);

  float absPFIsoComb_DR0p3(IsotrackObject const& part);
  float relPFIsoComb_DR0p3(IsotrackObject const& part);

  float computeIso(IsotrackObject const& part);

  void setSelectionBits(IsotrackObject& part);

}


#endif
