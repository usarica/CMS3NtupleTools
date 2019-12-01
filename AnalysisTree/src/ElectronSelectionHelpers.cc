#include <cassert>
#include "ElectronSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


float ElectronSelectionHelpers::absMiniIso_DR0p3(ElectronObject const& part){ return part.extras.miniIso_comb_nofsr; }
float ElectronSelectionHelpers::relMiniIso_DR0p3(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIso_DR0p3(part)/pt : 0.f); }

float ElectronSelectionHelpers::absPFIso_DR0p3(ElectronObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float ElectronSelectionHelpers::relPFIso_DR0p3(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

float ElectronSelectionHelpers::absPFIso_DR0p4(ElectronObject const& part){ return part.extras.pfIso04_comb_nofsr; }
float ElectronSelectionHelpers::relPFIso_DR0p4(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p4(part)/pt : 0.f); }

#define ID_PASS_VETO id_MVA_Fall17V2_NoIso_pass_wpLoose
#define ID_PASS_LOOSE id_MVA_Fall17V2_NoIso_pass_wpLoose
#define ID_PASS_MEDIUM id_MVA_Fall17V2_NoIso_pass_wp90
#define ID_PASS_TIGHT id_MVA_Fall17V2_NoIso_pass_wp80

bool ElectronSelectionHelpers::testVetoId(ElectronObject const& part){
  return part.extras.ID_PASS_VETO;
}
bool ElectronSelectionHelpers::testLooseId(ElectronObject const& part){
  return part.extras.ID_PASS_LOOSE;
}
bool ElectronSelectionHelpers::testMediumId(ElectronObject const& part){
  return part.extras.ID_PASS_MEDIUM;
}
bool ElectronSelectionHelpers::testTightId(ElectronObject const& part){
  return part.extras.ID_PASS_TIGHT;
}

#define ISO_FCN ElectronSelectionHelpers::relPFIso_DR0p3

bool ElectronSelectionHelpers::testVetoIso(ElectronObject const& part){
  return (ISO_FCN(part)<isoThr_veto);
}
bool ElectronSelectionHelpers::testLooseIso(ElectronObject const& part){
  return (ISO_FCN(part)<isoThr_loose);
}
bool ElectronSelectionHelpers::testMediumIso(ElectronObject const& part){
  return (ISO_FCN(part)<isoThr_medium);
}
bool ElectronSelectionHelpers::testTightIso(ElectronObject const& part){
  return (ISO_FCN(part)<isoThr_tight);
}

bool ElectronSelectionHelpers::testVetoKin(ElectronObject const& part){
  return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
}
bool ElectronSelectionHelpers::testLooseKin(ElectronObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool ElectronSelectionHelpers::testMediumKin(ElectronObject const& part){
  return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
}
bool ElectronSelectionHelpers::testTightKin(ElectronObject const& part){
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

bool ElectronSelectionHelpers::testPtEtaGen(ElectronObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool ElectronSelectionHelpers::testPreselection(ElectronObject const& part){
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
void ElectronSelectionHelpers::setSelectionBits(ElectronObject& part){
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
