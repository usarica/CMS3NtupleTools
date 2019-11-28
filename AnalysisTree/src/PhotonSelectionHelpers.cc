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
bool PhotonSelectionHelpers::testVetoIdIso(PhotonObject const& part){
  return (testVetoId(part)/* && ISO_FCN(part)<isoThr_veto*/);
}
bool PhotonSelectionHelpers::testLooseIdIso(PhotonObject const& part){
  return (testLooseId(part)/* && ISO_FCN(part)<isoThr_loose*/);
}
bool PhotonSelectionHelpers::testMediumIdIso(PhotonObject const& part){
  return (testMediumId(part)/* && ISO_FCN(part)<isoThr_medium*/);
}
bool PhotonSelectionHelpers::testTightIdIso(PhotonObject const& part){
  return (testTightId(part)/* && ISO_FCN(part)<isoThr_tight*/);
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
bool PhotonSelectionHelpers::testPtEtaSkim(PhotonObject const& part){
  // pT and eta skim cut
  if (testTightIdIso(part)) return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
  else if (testMediumIdIso(part)) return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
  else if (testLooseIdIso(part)) return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
  else if (testVetoIdIso(part)) return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
  else return false;
}
bool PhotonSelectionHelpers::testPreselection(PhotonObject const& part){
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
void PhotonSelectionHelpers::setSelectionBits(PhotonObject& part){
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
