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

bool ElectronSelectionHelpers::testVetoIdIso(ElectronObject const& part){
  return (testVetoId(part) && ISO_FCN(part)<isoThr_veto);
}
bool ElectronSelectionHelpers::testLooseIdIso(ElectronObject const& part){
  return (testLooseId(part) && ISO_FCN(part)<isoThr_loose);
}
bool ElectronSelectionHelpers::testMediumIdIso(ElectronObject const& part){
  return (testMediumId(part) && ISO_FCN(part)<isoThr_medium);
}
bool ElectronSelectionHelpers::testTightIdIso(ElectronObject const& part){
  return (testTightId(part) && ISO_FCN(part)<isoThr_tight);
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
bool ElectronSelectionHelpers::testPtEtaSkim(ElectronObject const& part){
  // pT and eta skim cut
  if (testTightIdIso(part)) return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
  else if (testMediumIdIso(part)) return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
  else if (testLooseIdIso(part)) return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
  else if (testVetoIdIso(part)) return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
  else return false;
}
bool ElectronSelectionHelpers::testPreselection(ElectronObject const& part){
  return (
    (
    (bit_preselection_idisoreco == kVetoIDIso && testVetoIdIso(part))
      ||
      (bit_preselection_idisoreco == kLooseIDIso && testLooseIdIso(part))
      ||
      (bit_preselection_idisoreco == kMediumIDIso && testMediumIdIso(part))
      ||
      (bit_preselection_idisoreco == kTightIDIso && testTightIdIso(part))
      ) && testPtEtaSkim(part)
    );
}
void ElectronSelectionHelpers::setSelectionBits(ElectronObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);
  if (testVetoId(part)) part.setSelectionBit(kVetoID);
  if (testVetoIdIso(part)) part.setSelectionBit(kVetoIDIso);
  if (testLooseId(part)) part.setSelectionBit(kLooseID);
  if (testLooseIdIso(part)) part.setSelectionBit(kLooseIDIso);
  if (testMediumId(part)) part.setSelectionBit(kMediumID);
  if (testMediumIdIso(part)) part.setSelectionBit(kMediumIDIso);
  if (testTightId(part)) part.setSelectionBit(kTightID);
  if (testTightIdIso(part)) part.setSelectionBit(kTightIDIso);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
