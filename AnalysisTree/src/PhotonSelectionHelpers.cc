#include <cassert>
#include "PhotonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


float PhotonSelectionHelpers::absPFIso_DR0p3(PhotonObject const& part){ return part.extras.pfIso_comb; }
float PhotonSelectionHelpers::relPFIso_DR0p3(PhotonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

#define ID_PASS_VETO id_cutBased_Fall17V2_Loose_Bits
#define ID_PASS_LOOSE id_cutBased_Fall17V2_Loose_Bits
#define ID_PASS_MEDIUM id_cutBased_Fall17V2_Medium_Bits
#define ID_PASS_TIGHT id_cutBased_Fall17V2_Tight_Bits

bool PhotonSelectionHelpers::testVetoId(PhotonObject const& part){
  return (part.extras.ID_PASS_VETO == 127);
}
bool PhotonSelectionHelpers::testLooseId(PhotonObject const& part){
  return (part.extras.ID_PASS_LOOSE == 127);
}
bool PhotonSelectionHelpers::testMediumId(PhotonObject const& part){
  return (part.extras.ID_PASS_MEDIUM == 127);
}
bool PhotonSelectionHelpers::testTightId(PhotonObject const& part){
  return (part.extras.ID_PASS_TIGHT == 127);
}

#define ISO_FCN PhotonSelectionHelpers::relPFIso_DR0p3

// Cut-based id applies isolation veto as well, so no need to apply it again...
bool PhotonSelectionHelpers::testVetoIso(PhotonObject const& part){
  return (true/* && ISO_FCN(part)<isoThr_veto*/);
}
bool PhotonSelectionHelpers::testLooseIso(PhotonObject const& part){
  return (true/* && ISO_FCN(part)<isoThr_loose*/);
}
bool PhotonSelectionHelpers::testMediumIso(PhotonObject const& part){
  return (true/* && ISO_FCN(part)<isoThr_medium*/);
}
bool PhotonSelectionHelpers::testTightIso(PhotonObject const& part){
  return (true/* && ISO_FCN(part)<isoThr_tight*/);
}

bool PhotonSelectionHelpers::testVetoKin(PhotonObject const& part){
  return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
}
bool PhotonSelectionHelpers::testLooseKin(PhotonObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool PhotonSelectionHelpers::testMediumKin(PhotonObject const& part){
  return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
}
bool PhotonSelectionHelpers::testTightKin(PhotonObject const& part){
  return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
}

#ifdef ID_PASS_VETO
#undef ID_PASS_VETO
#endif
#ifdef ID_PASS_LOOSE
#undef ID_PASS_LOOSE
#endif
#ifdef ID_PASS_MEDIUM
#undef ID_PASS_MEDIUM
#endif
#ifdef ID_PASS_TIGHT
#undef ID_PASS_TIGHT
#endif
#ifdef ISO_FCN
#undef ISO_FCN
#endif

bool PhotonSelectionHelpers::testPtEtaGen(PhotonObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool PhotonSelectionHelpers::testPreselection(PhotonObject const& part){
  return (
    (
    (bit_preselection_iso == kVetoIso && testVetoIso(part))
      ||
      (bit_preselection_iso == kLooseIso && testLooseIso(part))
      ||
      (bit_preselection_iso == kMediumIso && testMediumIso(part))
      ||
      (bit_preselection_iso == kTightIso && testTightIso(part))
      )
    &&
    (
    (bit_preselection_id == kVetoId && testVetoId(part))
      ||
      (bit_preselection_id == kLooseId && testLooseId(part))
      ||
      (bit_preselection_id == kMediumId && testMediumId(part))
      ||
      (bit_preselection_id == kTightId && testTightId(part))
      )
    &&
    (
    (bit_preselection_kin == kVetoKin && testVetoKin(part))
      ||
      (bit_preselection_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselection_kin == kMediumKin && testMediumKin(part))
      ||
      (bit_preselection_kin == kTightKin && testTightKin(part))
      )
    );
}
void PhotonSelectionHelpers::setSelectionBits(PhotonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);

  if (testVetoId(part)) part.setSelectionBit(kVetoId);
  if (testVetoIso(part)) part.setSelectionBit(kVetoIso);
  if (testVetoKin(part)) part.setSelectionBit(kVetoKin);

  if (testLooseId(part)) part.setSelectionBit(kLooseId);
  if (testLooseIso(part)) part.setSelectionBit(kLooseIso);
  if (testLooseKin(part)) part.setSelectionBit(kLooseKin);

  if (testMediumId(part)) part.setSelectionBit(kMediumId);
  if (testMediumIso(part)) part.setSelectionBit(kMediumIso);
  if (testMediumKin(part)) part.setSelectionBit(kMediumKin);

  if (testTightId(part)) part.setSelectionBit(kTightId);
  if (testTightIso(part)) part.setSelectionBit(kTightIso);
  if (testTightKin(part)) part.setSelectionBit(kTightKin);

  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
