#include <cassert>
#include "ElectronSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


float ElectronSelectionHelpers::absMiniIso_DR0p3(ElectronObject const& part){ return part.extras.miniIso_comb_nofsr; }
float ElectronSelectionHelpers::relMiniIso_DR0p3(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIso_DR0p3(part)/pt : 0.f); }

float ElectronSelectionHelpers::absPFIso_DR0p3(ElectronObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float ElectronSelectionHelpers::relPFIso_DR0p3(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

float ElectronSelectionHelpers::absPFIso_DR0p4(ElectronObject const& part){ return part.extras.pfIso04_comb_nofsr; }
float ElectronSelectionHelpers::relPFIso_DR0p4(ElectronObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p4(part)/pt : 0.f); }

float ElectronSelectionHelpers::computeIso(ElectronObject const& part){
  // If id is an MVA id with iso., return lowest possible iso;
  // otherwise, return the output of the corresponding function.
  if (
    idType_preselection == kMVAId_Fall17V2_Iso
    ||
    idType_preselection == kMVAId_HZZRun2Legacy_Iso
    ) return 0.f;
  else if (isoType_preselection == kMiniIsoDR0p3) return relMiniIso_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p3) return relPFIso_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p4) return relPFIso_DR0p4(part);
  else MELAerr << "ElectronSelectionHelpers::computeIso: Isolation " << isoType_preselection << " with id " << idType_preselection << " is not implemented." << endl;
  return 999.f;
}

/*
From https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Applying_Individual_Cuts_of_a_Se
For cut-based selection, the bit map is the following:
0: Min. pT cut
1: SC eta multi. range
2: dEtaIn seed
3: dPhiIn
4: Full 5x5 sigmaIetaIeta
5: H/E
6: 1/E - 1/p
7: Eff. area PF iso.
8: Conversion veto
9: Missing hits
*/
#define TEST_CUTBASED_BIT(ibit) (ibit == 1023 || ibit == 895)
bool ElectronSelectionHelpers::testVetoId(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Veto_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Veto_Bits);
  case kMVAId_Fall17V2_NoIso:
    return part.extras.id_MVA_Fall17V2_NoIso_pass_wpLoose;
  case kMVAId_Fall17V2_Iso:
    return part.extras.id_MVA_Fall17V2_Iso_pass_wpLoose;
  case kMVAId_HZZRun2Legacy_Iso:
    return part.extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
  default:
    MELAerr << "ElectronSelectionHelpers::testVetoId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool ElectronSelectionHelpers::testLooseId(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Loose_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2_NoIso:
    return part.extras.id_MVA_Fall17V2_NoIso_pass_wpLoose;
  case kMVAId_Fall17V2_Iso:
    return part.extras.id_MVA_Fall17V2_Iso_pass_wpLoose;
  case kMVAId_HZZRun2Legacy_Iso:
    return part.extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
  default:
    MELAerr << "ElectronSelectionHelpers::testLooseId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool ElectronSelectionHelpers::testMediumId(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Medium_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Medium_Bits);
  case kMVAId_Fall17V2_NoIso:
    return part.extras.id_MVA_Fall17V2_NoIso_pass_wp90;
  case kMVAId_Fall17V2_Iso:
    return part.extras.id_MVA_Fall17V2_Iso_pass_wp90;
  case kMVAId_HZZRun2Legacy_Iso:
    return part.extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
  default:
    MELAerr << "ElectronSelectionHelpers::testMediumId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool ElectronSelectionHelpers::testTightId(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Tight_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Tight_Bits);
  case kMVAId_Fall17V2_NoIso:
    return part.extras.id_MVA_Fall17V2_NoIso_pass_wp80;
  case kMVAId_Fall17V2_Iso:
    return part.extras.id_MVA_Fall17V2_Iso_pass_wp80;
  case kMVAId_HZZRun2Legacy_Iso:
    return part.extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
  default:
    MELAerr << "ElectronSelectionHelpers::testTightId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
#undef TEST_CUTBASED_BIT

#define TEST_CUTBASED_BIT(ibit) HelperFunctions::test_bit(ibit, 7)
bool ElectronSelectionHelpers::testVetoIso(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Veto_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Veto_Bits);
  case kMVAId_Fall17V2_Iso:
  case kMVAId_HZZRun2Legacy_Iso:
    return true;
  default:
    return (computeIso(part)<isoThr_veto);
  };
}
bool ElectronSelectionHelpers::testLooseIso(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Loose_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Loose_Bits);
  case kMVAId_Fall17V2_Iso:
  case kMVAId_HZZRun2Legacy_Iso:
    return true;
  default:
    return (computeIso(part)<isoThr_loose);
  };
}
bool ElectronSelectionHelpers::testMediumIso(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Medium_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Medium_Bits);
  case kMVAId_Fall17V2_Iso:
  case kMVAId_HZZRun2Legacy_Iso:
    return true;
  default:
    return (computeIso(part)<isoThr_medium);
  };
}
bool ElectronSelectionHelpers::testTightIso(ElectronObject const& part){
  switch (idType_preselection){
  case kCutBasedId_Fall17V1:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V1_Tight_Bits);
  case kCutBasedId_Fall17V2:
    return TEST_CUTBASED_BIT(part.extras.id_cutBased_Fall17V2_Tight_Bits);
  case kMVAId_Fall17V2_Iso:
  case kMVAId_HZZRun2Legacy_Iso:
    return true;
  default:
    return (computeIso(part)<isoThr_tight);
  };
}
#undef TEST_CUTBASED_BIT

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
