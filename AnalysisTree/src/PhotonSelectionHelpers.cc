#include <cassert>

#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>

#include "PhotonSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


namespace PhotonSelectionHelpers{
  bool testPtEtaGen(PhotonObject const& part);

  bool testConversionSafe(PhotonObject const& part);

  bool testInTimeSeed(PhotonObject const& part);
  bool testBeamHaloSafe(PhotonObject const& part);
  bool testSpikeSafe(PhotonObject const& part);

  bool testPFPhotonId(PhotonObject const& part);
  bool testPFMETSafe(PhotonObject const& part);

  bool testVetoId(PhotonObject const& part);
  bool testVetoIso(PhotonObject const& part);
  bool testVetoKin(PhotonObject const& part);

  bool testLooseId(PhotonObject const& part);
  bool testLooseIso(PhotonObject const& part);
  bool testLooseKin(PhotonObject const& part);

  bool testMediumId(PhotonObject const& part);
  bool testMediumIso(PhotonObject const& part);
  bool testMediumKin(PhotonObject const& part);

  bool testTightId(PhotonObject const& part);
  bool testTightIso(PhotonObject const& part);
  bool testTightKin(PhotonObject const& part);

  bool testPreselectionVeto(PhotonObject const& part);
  bool testPreselectionLoose(PhotonObject const& part);
  bool testPreselectionTight(PhotonObject const& part);

  bool testSFTampon(PhotonObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


float PhotonSelectionHelpers::getIsolationDRmax(PhotonObject const& /*part*/){
  if (isoType_preselection == kPFIsoDR0p3) return 0.3;
  else{
    MELAerr << "PhotonSelectionHelpers::getIsolationDRmax: Isolation type " << isoType_preselection << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

float PhotonSelectionHelpers::absPFIso_DR0p3(PhotonObject const& part){ return part.extras.pfIso_comb; }
float PhotonSelectionHelpers::relPFIso_DR0p3(PhotonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

bool PhotonSelectionHelpers::testConversionSafe(PhotonObject const& part){ return (!part.extras.hasPixelSeed && part.extras.passElectronVeto); }

bool PhotonSelectionHelpers::testInTimeSeed(PhotonObject const& part){ return std::abs(part.extras.seedTime)<seedTimeThr; }
bool PhotonSelectionHelpers::testBeamHaloSafe(PhotonObject const& part){ return part.extras.MIPTotalEnergy<mipTotalEnergyThr; }
bool PhotonSelectionHelpers::testSpikeSafe(PhotonObject const& part){
  return part.extras.full5x5_sigmaIEtaIEta<full5x5_sigmaIEtaIEtaThr && part.extras.full5x5_sigmaIPhiIPhi<full5x5_sigmaIPhiIPhiThr;
}

bool PhotonSelectionHelpers::testPFPhotonId(PhotonObject const& part){
  auto const& ibit = part.extras.id_egamma_pfPhoton_Bits;
  constexpr bool testBadHCAL = true;
  return (
    part.extras.n_associated_pfphotons==1
    &&
    HelperFunctions::test_bit(ibit, ISEGAMMAPFPHOTON_BASE)
    &&
    (!testBadHCAL || HelperFunctions::test_bit(ibit, ISEGAMMAPFPHOTON_BASE_BADHCALMITIGATED))
    &&
    part.extras.min_dR_photon_pfphoton_associated<mindRThr_photon_pfphoton
    );
}
bool PhotonSelectionHelpers::testPFMETSafe(PhotonObject const& part){ return HelperFunctions::test_bit(part.extras.id_egamma_pfPhoton_Bits, ISEGAMMAPFPHOTON_METSAFE); }

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
bool PhotonSelectionHelpers::testPreselectionVeto(PhotonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionVeto_id)
    &&
    part.testSelectionBit(bit_preselectionVeto_iso)
    &&
    part.testSelectionBit(bit_preselectionVeto_kin)
    );
}
bool PhotonSelectionHelpers::testPreselectionLoose(PhotonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_iso)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    &&
    (bit_preselection_conversion != kConversionSafe || part.testSelectionBit(bit_preselection_conversion))
    );
}
bool PhotonSelectionHelpers::testPreselectionTight(PhotonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(bit_preselectionTight_iso)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    &&
    (bit_preselection_conversion != kConversionSafe || part.testSelectionBit(bit_preselection_conversion))
    );
}
bool PhotonSelectionHelpers::testSFTampon(PhotonObject const& part){
  static_assert(bit_SFTampon_id <= bit_preselectionTight_id);
  static_assert(bit_SFTampon_iso <= bit_preselectionTight_iso);
  static_assert(bit_SFTampon_kin <= bit_preselectionTight_kin);
  return (
    part.testSelectionBit(bit_SFTampon_id)
    &&
    part.testSelectionBit(bit_SFTampon_iso)
    &&
    part.testSelectionBit(bit_SFTampon_kin)
    &&
    (bit_preselection_conversion != kConversionSafe || part.testSelectionBit(bit_preselection_conversion))
    );
}
void PhotonSelectionHelpers::setSelectionBits(PhotonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  part.setSelectionBit(kConversionSafe, testConversionSafe(part));

  part.setSelectionBit(kInTimeSeed, testInTimeSeed(part));
  part.setSelectionBit(kBeamHaloSafe, testBeamHaloSafe(part));
  part.setSelectionBit(kSpikeSafe, testSpikeSafe(part));

  part.setSelectionBit(kPFPhotonId, testPFPhotonId(part));
  part.setSelectionBit(kPFMETSafe, testPFMETSafe(part));

  part.setSelectionBit(kVetoId, testVetoId(part));
  part.setSelectionBit(kVetoIso, testVetoIso(part));
  part.setSelectionBit(kVetoKin, testVetoKin(part));

  part.setSelectionBit(kLooseId, testLooseId(part));
  part.setSelectionBit(kLooseIso, testLooseIso(part));
  part.setSelectionBit(kLooseKin, testLooseKin(part));

  part.setSelectionBit(kMediumId, testMediumId(part));
  part.setSelectionBit(kMediumIso, testMediumIso(part));
  part.setSelectionBit(kMediumKin, testMediumKin(part));

  part.setSelectionBit(kTightId, testTightId(part));
  part.setSelectionBit(kTightIso, testTightIso(part));
  part.setSelectionBit(kTightKin, testTightKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionVeto, testPreselectionVeto(part));
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));

  part.setSelectionBit(kSFTampon, testSFTampon(part));
}
