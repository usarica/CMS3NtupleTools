#include <cassert>
#include "PhotonSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


float PhotonSelectionHelpers::absPFIso_DR0p3(PhotonObject const& part){ return part.extras.pfIso_comb; }
float PhotonSelectionHelpers::relPFIso_DR0p3(PhotonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

/*
From https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Applying_Individual_Cuts_of_a_Se
For cut-based selection, the bit map is the following:
0: Min. pT cut
1: SC eta multi. range
2: Single tower H/E
3: Full 5x5 sigmaIetaIeta
4: Iso_ch 
5: Iso_nh
6: Iso_em
*/
#define TEST_CUTBASED_BIT(ibit) (HelperFunctions::test_bit(ibit, 0) && HelperFunctions::test_bit(ibit, 1) && HelperFunctions::test_bit(ibit, 2) && HelperFunctions::test_bit(ibit, 3))
bool PhotonSelectionHelpers::testVetoId(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2:
    return part.extras.id_MVA_Fall17V2_pass_wp90;
  default:
    MELAerr << "PhotonSelectionHelpers::testVetoId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testLooseId(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2:
    return part.extras.id_MVA_Fall17V2_pass_wp90;
  default:
    MELAerr << "PhotonSelectionHelpers::testLooseId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testMediumId(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Medium_Bits);
  case kMVAId_Fall17V2:
    return part.extras.id_MVA_Fall17V2_pass_wp80;
  default:
    MELAerr << "PhotonSelectionHelpers::testMediumId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testTightId(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Tight_Bits);
  case kMVAId_Fall17V2:
    return part.extras.id_MVA_Fall17V2_pass_wp80;
  default:
    MELAerr << "PhotonSelectionHelpers::testTightId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
#undef TEST_CUTBASED_BIT

// Cut-based id applies isolation veto as well, so no need to calculate it again...
// MVA id uses isolation as a variable (see slide 4 in https://indico.cern.ch/event/697079/contributions/2968123/attachments/1632966/2604131/PhotonID_EGM_13.04.2018.pdf#search=Kuntal%20Mondal),
// so isolation test should just return true here.
#define TEST_CUTBASED_BIT(ibit) (HelperFunctions::test_bit(ibit, 4) && HelperFunctions::test_bit(ibit, 5) && HelperFunctions::test_bit(ibit, 6))
bool PhotonSelectionHelpers::testVetoIso(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2:
    return true;
  default:
    MELAerr << "PhotonSelectionHelpers::testVetoIso: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testLooseIso(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2:
    return true;
  default:
    MELAerr << "PhotonSelectionHelpers::testLooseIso: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testMediumIso(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Medium_Bits);
  case kMVAId_Fall17V2:
    return true;
  default:
    MELAerr << "PhotonSelectionHelpers::testMediumIso: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool PhotonSelectionHelpers::testTightIso(PhotonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Tight_Bits);
  case kMVAId_Fall17V2:
    return true;
  default:
    MELAerr << "PhotonSelectionHelpers::testTightIso: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
#undef TEST_CUTBASED_BIT

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
