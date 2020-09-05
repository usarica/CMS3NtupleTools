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

  bool testLooseId(AK4JetObject const& part);
  bool testLooseKin(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightKin(AK4JetObject const& part);

  bool testBtaggable(AK4JetObject const& part);

  bool testPreselectionLoose(AK4JetObject const& part);
  bool testPreselectionTight(AK4JetObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


bool AK4JetSelectionHelpers::testLoosePUJetId(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId, 2);
}
bool AK4JetSelectionHelpers::testMediumPUJetId(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId, 1);
}
bool AK4JetSelectionHelpers::testTightPUJetId(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId, 0);
}

bool AK4JetSelectionHelpers::testLoosePUJetId_Default(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId_default, 2);
}
bool AK4JetSelectionHelpers::testMediumPUJetId_Default(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId_default, 1);
}
bool AK4JetSelectionHelpers::testTightPUJetId_Default(AK4JetObject const& part){
  return HelperFunctions::test_bit(part.extras.pileupJetId_default, 0);
}

bool AK4JetSelectionHelpers::testTightLeptonVetoId(AK4JetObject const& part){ return part.extras.pass_leptonVetoId; }

bool AK4JetSelectionHelpers::testLooseId(AK4JetObject const& part){ return part.extras.pass_looseId; }
bool AK4JetSelectionHelpers::testTightId(AK4JetObject const& part){ return part.extras.pass_tightId; }

bool AK4JetSelectionHelpers::testLooseKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && std::abs(part.eta())<etaThr_skim_loose);
}
bool AK4JetSelectionHelpers::testTightKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && std::abs(part.eta())<etaThr_skim_tight);
}

bool AK4JetSelectionHelpers::testPtEtaGen(AK4JetObject const& part){
  return (part.pt()>=ptThr_gen && std::abs(part.eta())<etaThr_gen);
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

bool AK4JetSelectionHelpers::testBtaggable(AK4JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.pt()>=ptThr_skim_btag && std::abs(part.eta())<(SampleHelpers::theDataYear<=2016 ? 2.4f : 2.5f)
    &&
    (PUIdWP==nSelectionBits || part.testSelectionBit(PUIdWP))
    &&
    (!applyTightLeptonVetoIdToJets || part.testSelectionBit(kTightLeptonVetoId))
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

void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  part.setSelectionBit(kLoosePUJetId, testLoosePUJetId(part));
  part.setSelectionBit(kMediumPUJetId, testMediumPUJetId(part));
  part.setSelectionBit(kTightPUJetId, testTightPUJetId(part));

  part.setSelectionBit(kLoosePUJetId_default, testLoosePUJetId_Default(part));
  part.setSelectionBit(kMediumPUJetId_default, testMediumPUJetId_Default(part));
  part.setSelectionBit(kTightPUJetId_default, testTightPUJetId_Default(part));

  part.setSelectionBit(kTightLeptonVetoId, testTightLeptonVetoId(part));

  part.setSelectionBit(kLooseId, testLooseId(part));
  part.setSelectionBit(kLooseKin, testLooseKin(part));

  part.setSelectionBit(kTightId, testTightId(part));
  part.setSelectionBit(kTightKin, testTightKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kBtaggable, testBtaggable(part));
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
