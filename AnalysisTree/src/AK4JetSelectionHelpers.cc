#include <cassert>
#include <cmath>

#include "SamplesCore.h"
#include "AK4JetSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


namespace AK4JetSelectionHelpers{
  SelectionBits PUIdWP = kTightPUJetId;
  bool applyTightLeptonVetoIdToJets = false;

  bool testLoosePUJetId(AK4JetObject const& part);
  bool testMediumPUJetId(AK4JetObject const& part);
  bool testTightPUJetId(AK4JetObject const& part);

  bool testLoosePUJetId_Default(AK4JetObject const& part);
  bool testMediumPUJetId_Default(AK4JetObject const& part);
  bool testTightPUJetId_Default(AK4JetObject const& part);

  bool testTightLeptonVetoId(AK4JetObject const& part);

  bool testPtEtaGen(AK4JetObject const& part);

  bool testNoisyJet(AK4JetObject const& part);

  bool testMETJERCSafe_Base(AK4JetObject const& part);

  bool testLooseId(AK4JetObject const& part);
  bool testLooseKin(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightKin(AK4JetObject const& part);

  bool testBtaggableKin(AK4JetObject const& part);
  bool testBtaggable_NoPUJetId(AK4JetObject const& part);
  bool testBtaggable(AK4JetObject const& part);

  bool testMETFixRecoverySuitable(AK4JetObject const& part);

  bool testPreselectionLoose(AK4JetObject const& part);
  bool testPreselectionTight(AK4JetObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


bool const& AK4JetSelectionHelpers::getApplyTightLeptonVetoIdToJetsFlag(){ return applyTightLeptonVetoIdToJets; }

bool AK4JetSelectionHelpers::testLoosePUJetId(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId, 2);
}
bool AK4JetSelectionHelpers::testMediumPUJetId(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId, 1);
}
bool AK4JetSelectionHelpers::testTightPUJetId(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId, 0);
}

bool AK4JetSelectionHelpers::testLoosePUJetId_Default(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId_default, 2);
}
bool AK4JetSelectionHelpers::testMediumPUJetId_Default(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId_default, 1);
}
bool AK4JetSelectionHelpers::testTightPUJetId_Default(AK4JetObject const& part){
  return part.pt()>=ptThr_PUId || std::abs(part.eta())>=etaThr_PUId || HelperFunctions::test_bit(part.extras.pileupJetId_default, 0);
}

bool AK4JetSelectionHelpers::testTightLeptonVetoId(AK4JetObject const& part){ return part.extras.pass_leptonVetoId; }

bool AK4JetSelectionHelpers::testNoisyJet(AK4JetObject const& part){
  int const& year = SampleHelpers::getDataYear();
  if (year!=2017 && year!=2018) return false;
  float const abs_eta = std::abs(part.eta());
  float const thr_NEMF = (year==2017 ? 0.5 : 0.7);
  return (abs_eta>=2.65 && abs_eta<=3.139 && part.extras.NEMF>=thr_NEMF);
}

bool AK4JetSelectionHelpers::testMETJERCSafe_Base(AK4JetObject const& part){
  bool const hasMETJERCSafeEM = (part.extras.CEMF+part.extras.NEMF<=0.9);
  // Do not apply pT cut or muon subtraction here; they depend on the p4-preservation and JER application.
  // eta is safe to apply because the corrections do not change eta.
  return hasMETJERCSafeEM && std::abs(part.eta())<=9.9;
}

bool AK4JetSelectionHelpers::testLooseId(AK4JetObject const& part){ return part.extras.pass_looseId; }
bool AK4JetSelectionHelpers::testTightId(AK4JetObject const& part){ return part.extras.pass_tightId; }

bool AK4JetSelectionHelpers::testPtEtaGen(AK4JetObject const& part){
  return (part.pt()>=ptThr_gen && std::abs(part.eta())<etaThr_gen);
}
bool AK4JetSelectionHelpers::testLooseKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && std::abs(part.eta())<etaThr_skim_loose);
}
bool AK4JetSelectionHelpers::testTightKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && std::abs(part.eta())<etaThr_skim_tight);
}

bool AK4JetSelectionHelpers::testPreselectionLoose(AK4JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    );
}
bool AK4JetSelectionHelpers::testPreselectionTight(AK4JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    &&
    (PUIdWP==nSelectionBits || part.testSelectionBit(PUIdWP))
    &&
    (!applyTightLeptonVetoIdToJets || part.testSelectionBit(kTightLeptonVetoId))
    );
}

bool AK4JetSelectionHelpers::testBtaggableKin(AK4JetObject const& part){
  return (
    part.pt()>=ptThr_btag && std::abs(part.eta())<(SampleHelpers::theDataYear<=2016 ? 2.4f : 2.5f) // For jets out of tracker acceptance, b-tagging is meanningless.
    );
}
bool AK4JetSelectionHelpers::testBtaggable_NoPUJetId(AK4JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    (!applyTightLeptonVetoIdToJets || part.testSelectionBit(kTightLeptonVetoId))
    &&
    part.testSelectionBit(kBtaggableKin)
    );
}
bool AK4JetSelectionHelpers::testBtaggable(AK4JetObject const& part){
  return (
    part.testSelectionBit(kBtaggable_NoPUJetId)
    &&
    (PUIdWP==nSelectionBits || part.testSelectionBit(PUIdWP))    );
}

bool AK4JetSelectionHelpers::testMETFixRecoverySuitable(AK4JetObject const& part){
  if (SampleHelpers::getDataYear()!=2017) return false;

  bool const useDefaultPUJetId = (PUIdWP>=kLoosePUJetId_default && PUIdWP<=kTightPUJetId_default);
  float const abs_eta = std::abs(part.eta());
  float const uncorrected_pt = part.uncorrected_p4().Pt();
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    &&
    part.testSelectionBit((!useDefaultPUJetId ? kTightPUJetId : kTightPUJetId_default)) // Must use tight PU jet id
    &&
    (abs_eta>=2.65 && abs_eta<=3.139 && uncorrected_pt<=50. && !part.testSelectionBit(kNoisyJet))
    );
}

void AK4JetSelectionHelpers::setPUIdWP(SelectionBits flag){
  if (
    flag==kLoosePUJetId || flag==kMediumPUJetId || flag==kTightPUJetId
    ||
    flag==kLoosePUJetId_default || flag==kMediumPUJetId_default || flag==kTightPUJetId_default
    ) PUIdWP = flag;
  else PUIdWP = nSelectionBits;
}
void AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(bool flag){ applyTightLeptonVetoIdToJets = flag; }

void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part, bool resetIDs, bool resetKinematics){
  static_assert(std::numeric_limits<ParticleObject::SelectionBitsType_t>::digits >= nSelectionBits);

  if (resetKinematics) part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  if (resetIDs) part.setSelectionBit(kLoosePUJetId, testLoosePUJetId(part));
  if (resetIDs) part.setSelectionBit(kMediumPUJetId, testMediumPUJetId(part));
  if (resetIDs) part.setSelectionBit(kTightPUJetId, testTightPUJetId(part));

  if (resetIDs) part.setSelectionBit(kLoosePUJetId_default, testLoosePUJetId_Default(part));
  if (resetIDs) part.setSelectionBit(kMediumPUJetId_default, testMediumPUJetId_Default(part));
  if (resetIDs) part.setSelectionBit(kTightPUJetId_default, testTightPUJetId_Default(part));

  if (resetIDs) part.setSelectionBit(kTightLeptonVetoId, testTightLeptonVetoId(part));

  if (resetIDs) part.setSelectionBit(kNoisyJet, testNoisyJet(part));

  if (resetIDs) part.setSelectionBit(kMETJERCSafe_Base, testMETJERCSafe_Base(part));

  if (resetIDs) part.setSelectionBit(kLooseId, testLooseId(part));
  if (resetKinematics) part.setSelectionBit(kLooseKin, testLooseKin(part));

  if (resetIDs) part.setSelectionBit(kTightId, testTightId(part));
  if (resetKinematics) part.setSelectionBit(kTightKin, testTightKin(part));

  // The functions below test the bits set in the steps above.
  // The b-taggable function tests pT and eta, but they should pertain to values before any external jet modifications. Hence the use of 'resetIDs'.
  // The rest should always be re-tested.
  if (resetKinematics) part.setSelectionBit(kBtaggableKin, testBtaggableKin(part));
  if (resetIDs){
    part.setSelectionBit(kBtaggable_NoPUJetId, testBtaggable_NoPUJetId(part)); // Tests kBtaggableKin
    part.setSelectionBit(kBtaggable, testBtaggable(part)); // Tests kBtaggable_NoPUJetId
  }

  // This function uses the PU jet id bits above.
  if (resetIDs) part.setSelectionBit(kMETFixRecoverySuitable, testMETFixRecoverySuitable(part));

  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
