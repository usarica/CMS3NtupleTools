#include <cassert>
#include <cmath>

#include <DataFormats/MuonReco/interface/Muon.h>

#include <CMS3/Dictionaries/interface/MuonEnums.h>

#include "MuonSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


// These are functions hidden from the user
namespace MuonSelectionHelpers{
  bool allowProbeIdInLooseSelection = false;
  bool allowFakeableInLooseSelection = false;
  float isoThr_fakeable_trkIso = -1;

  bool testPtEtaGen(MuonObject const& part);

  bool testMuonSystemTime(MuonObject const& part);

  bool testVetoId(MuonObject const& part);
  bool testVetoIso(MuonObject const& part);
  bool testVetoKin(MuonObject const& part);

  bool testLooseId(MuonObject const& part);
  bool testLooseIso(MuonObject const& part);
  bool testLooseKin(MuonObject const& part);

  bool testMediumId(MuonObject const& part);
  bool testMediumIso(MuonObject const& part);
  bool testMediumKin(MuonObject const& part);

  bool testTightId(MuonObject const& part);
  bool testTightIso(MuonObject const& part);
  bool testTightKin(MuonObject const& part);

  bool testSoftId(MuonObject const& part);
  bool testSoftIso(MuonObject const& part);
  bool testSoftKin(MuonObject const& part);

  bool testProbeId(MuonObject const& part);
  bool testProbeSTAId(MuonObject const& part);

  bool testFakeableBaseIso(MuonObject const& part);
  bool testFakeableBase(MuonObject const& part);
  bool testFakeable(MuonObject const& part);

  bool testPreselectionVeto(MuonObject const& part);
  bool testPreselectionLoose_NoIso(MuonObject const& part);
  bool testPreselectionLoose(MuonObject const& part);
  bool testPreselectionTight(MuonObject const& part);
}


using namespace std;
using namespace reco;
using namespace MuonEnums;
using namespace MELAStreamHelpers;


void MuonSelectionHelpers::setAllowProbeIdInLooseSelection(bool flag){ allowProbeIdInLooseSelection = flag; }
void MuonSelectionHelpers::setAllowFakeableInLooseSelection(bool flag){ allowFakeableInLooseSelection = flag; }
void MuonSelectionHelpers::doRequireTrackerIsolationInFakeable(float const& isothr){ isoThr_fakeable_trkIso = isothr; }

bool MuonSelectionHelpers::getAllowProbeIdInLooseSelection(){ return allowProbeIdInLooseSelection; }
bool MuonSelectionHelpers::getAllowFakeableInLooseSelection(){ return allowFakeableInLooseSelection; }

bool MuonSelectionHelpers::testGoodMETPFMuon(PFCandidateObject const& part){
  // The following selection requirements should be identical to NtupleMaker::MuonSelectionHelpers:
  return (part.isGlobalMuon() || part.isStandaloneMuon());
}

float MuonSelectionHelpers::absPFIso_DR0p3(MuonObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_EACorr_DR0p3(MuonObject const& part){ return (part.extras.pfIso03_sum_charged_nofsr + std::max(0.f, part.extras.pfIso03_sum_neutral_EAcorr_nofsr)); }
float MuonSelectionHelpers::relPFIso_EACorr_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_EACorr_DR0p3(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_DR0p4(MuonObject const& part){ return part.extras.pfIso04_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p4(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p4(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_EACorr_DR0p4(MuonObject const& part){ return (part.extras.pfIso04_sum_charged_nofsr + std::max(0.f, part.extras.pfIso04_sum_neutral_EAcorr_nofsr)); }
float MuonSelectionHelpers::relPFIso_EACorr_DR0p4(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_EACorr_DR0p4(part)/pt : 0.f); }

float MuonSelectionHelpers::absMiniIso(MuonObject const& part){ return part.extras.miniIso_comb_nofsr; }
float MuonSelectionHelpers::relMiniIso(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIso(part)/pt : 0.f); }

float MuonSelectionHelpers::getIsolationDRmax(MuonObject const& part){
  if (isoType_preselection == kPFIsoDR0p3 || isoType_preselection == kPFIsoDR0p3_EACorrected) return 0.3;
  else if (isoType_preselection == kPFIsoDR0p4 || isoType_preselection == kPFIsoDR0p4_EACorrected) return 0.4;
  else if (isoType_preselection == kMiniIso) return (10. / std::min(std::max(part.uncorrected_pt(), 50.), 200.));
  else{
    MELAerr << "MuonSelectionHelpers::getIsolationDRmax: Isolation type " << isoType_preselection << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

float MuonSelectionHelpers::computeIso(MuonObject const& part){
  if (isoType_preselection == kPFIsoDR0p3) return relPFIso_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p3_EACorrected) return relPFIso_EACorr_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p4) return relPFIso_DR0p4(part);
  else if (isoType_preselection == kPFIsoDR0p4_EACorrected) return relPFIso_EACorr_DR0p4(part);
  else if (isoType_preselection == kMiniIso) return relMiniIso(part);
  else MELAerr << "MuonSelectionHelpers::computeIso: Isolation " << isoType_preselection << " with id " << idType_preselection << " is not implemented." << endl;
  return 999.f;
}

#define ID_CUTBASED_MUONTIME reco::Muon::InTimeMuon
#define ID_CUTBASED_VETO reco::Muon::CutBasedIdLoose
#define ID_CUTBASED_LOOSE reco::Muon::CutBasedIdLoose
#define ID_CUTBASED_MEDIUM reco::Muon::CutBasedIdMediumPrompt
#define ID_CUTBASED_TIGHT reco::Muon::CutBasedIdTight
bool MuonSelectionHelpers::testMuonSystemTime(MuonObject const& part){
  return ((part.extras.POG_selector_bits & ID_CUTBASED_MUONTIME) == ID_CUTBASED_MUONTIME);
}
bool MuonSelectionHelpers::testVetoId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & ID_CUTBASED_VETO) == ID_CUTBASED_VETO);
  case kCutBasedId_H4l:
    return (
      HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_Minimal) && HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_SIP3D)
      &&
      (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_PFID) || (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_HighPt) && part.pt()>200.f))
      );
  default:
    MELAerr << "MuonSelectionHelpers::testVetoId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testLooseId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & ID_CUTBASED_LOOSE) == ID_CUTBASED_LOOSE);
  case kCutBasedId_H4l:
    return (
      HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_Minimal) && HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_SIP3D)
      &&
      (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_PFID) || (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_HighPt) && part.pt()>200.f))
      );
  default:
    MELAerr << "MuonSelectionHelpers::testLooseId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testMediumId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & ID_CUTBASED_MEDIUM) == ID_CUTBASED_MEDIUM);
  case kCutBasedId_H4l:
    return (
      HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_Minimal) && HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_SIP3D)
      &&
      (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_PFID) || (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_HighPt) && part.pt()>200.f))
      );
  default:
    MELAerr << "MuonSelectionHelpers::testMediumId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testTightId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & ID_CUTBASED_TIGHT) == ID_CUTBASED_TIGHT);
  case kCutBasedId_H4l:
    return (
      HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_Minimal) && HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_SIP3D)
      &&
      (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_PFID) || (HelperFunctions::test_bit(part.extras.id_cutBased_H4l_Bits, kH4lSelection_HighPt) && part.pt()>200.f))
      );
  default:
    MELAerr << "MuonSelectionHelpers::testTightId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testSoftId(MuonObject const& part){ return ((part.extras.POG_selector_bits & Muon::SoftCutBasedId) == Muon::SoftCutBasedId); }
#undef ID_CUTBASED_VETO
#undef ID_CUTBASED_LOOSE
#undef ID_CUTBASED_MEDIUM
#undef ID_CUTBASED_TIGHT

bool MuonSelectionHelpers::testVetoIso(MuonObject const& part){ return (computeIso(part)<isoThr_veto); }
bool MuonSelectionHelpers::testLooseIso(MuonObject const& part){ return (computeIso(part)<isoThr_loose); }
bool MuonSelectionHelpers::testMediumIso(MuonObject const& part){ return (computeIso(part)<isoThr_medium); }
bool MuonSelectionHelpers::testTightIso(MuonObject const& part){ return (computeIso(part)<isoThr_tight); }
bool MuonSelectionHelpers::testSoftIso(MuonObject const& part){ return (computeIso(part)<isoThr_soft); }
bool MuonSelectionHelpers::testFakeableBaseIso(MuonObject const& part){ return (computeIso(part)<isoThr_fakeable); }

bool MuonSelectionHelpers::testVetoKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_veto && std::abs(part.eta())<etaThr_skim_veto);
}
bool MuonSelectionHelpers::testLooseKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_loose && std::abs(part.eta())<etaThr_skim_loose);
}
bool MuonSelectionHelpers::testMediumKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_medium && std::abs(part.eta())<etaThr_skim_medium);
}
bool MuonSelectionHelpers::testTightKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_tight && std::abs(part.eta())<etaThr_skim_tight);
}
bool MuonSelectionHelpers::testSoftKin(MuonObject const& part){ return (part.pt()>=ptThr_skim_soft && std::abs(part.eta())<etaThr_skim_soft); }

bool MuonSelectionHelpers::testProbeId(MuonObject const& part){ return part.extras.is_probeForTnP; }
bool MuonSelectionHelpers::testProbeSTAId(MuonObject const& part){ return part.extras.is_probeForTnP_STA; }

bool MuonSelectionHelpers::testPtEtaGen(MuonObject const& part){
  return (part.pt()>=ptThr_gen && std::abs(part.eta())<etaThr_gen);
}
bool MuonSelectionHelpers::testFakeableBase(MuonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(kFakeableBaseIso)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    &&
    (bit_preselection_time != kValidMuonSystemTime || part.testSelectionBit(bit_preselection_time))
    );
}
bool MuonSelectionHelpers::testFakeable(MuonObject const& part){
  return (
    part.testSelectionBit(kFakeableBase)
    &&
    (isoThr_fakeable_trkIso<0.f || part.extras.trkIso03_trackerSumPt<isoThr_fakeable_trkIso)
    );
}
bool MuonSelectionHelpers::testPreselectionVeto(MuonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionVeto_id)
    &&
    part.testSelectionBit(bit_preselectionVeto_iso)
    &&
    part.testSelectionBit(bit_preselectionVeto_kin)
    );
}
bool MuonSelectionHelpers::testPreselectionLoose_NoIso(MuonObject const& part){
  bool const isProbe = (!allowProbeIdInLooseSelection ? false : part.testSelectionBit(kProbeId));
  bool const isFakeable = (!allowFakeableInLooseSelection ? false : part.testSelectionBit(kFakeable));
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    &&
    (bit_preselection_time != kValidMuonSystemTime || part.testSelectionBit(bit_preselection_time))
    )
    ||
    isFakeable
    ||
    isProbe;
}
bool MuonSelectionHelpers::testPreselectionLoose(MuonObject const& part){
  bool const isProbe = (!allowProbeIdInLooseSelection ? false : part.testSelectionBit(kProbeId));
  bool const isFakeable = (!allowFakeableInLooseSelection ? false : part.testSelectionBit(kFakeable));
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_iso)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    &&
    (bit_preselection_time != kValidMuonSystemTime || part.testSelectionBit(bit_preselection_time))
    )
    ||
    isFakeable
    ||
    isProbe;
}
bool MuonSelectionHelpers::testPreselectionTight(MuonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(bit_preselectionTight_iso)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    &&
    (bit_preselection_time != kValidMuonSystemTime || part.testSelectionBit(bit_preselection_time))
    );
}
void MuonSelectionHelpers::setSelectionBits(MuonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  part.setSelectionBit(kValidMuonSystemTime, testMuonSystemTime(part));

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

  part.setSelectionBit(kSoftId, testSoftId(part));
  part.setSelectionBit(kSoftIso, testSoftIso(part));
  part.setSelectionBit(kSoftKin, testSoftKin(part));

  part.setSelectionBit(kProbeId, testProbeId(part));
  part.setSelectionBit(kProbeSTAId, testProbeSTAId(part));

  part.setSelectionBit(kFakeableBaseIso, testFakeableBaseIso(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kFakeableBase, testFakeableBase(part));
  part.setSelectionBit(kFakeable, testFakeable(part));
  part.setSelectionBit(kPreselectionVeto, testPreselectionVeto(part));
  part.setSelectionBit(kPreselectionLoose_NoIso, testPreselectionLoose_NoIso(part));
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
