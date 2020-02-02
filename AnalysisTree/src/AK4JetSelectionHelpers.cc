#include <cassert>
#include "AK4JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


namespace AK4JetSelectionHelpers{
  bool testPtEtaGen(AK4JetObject const& part);

  bool testPUJetId(AK4JetObject const& part);

  bool testLooseId(AK4JetObject const& part);
  bool testLooseKin(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightKin(AK4JetObject const& part);

  bool testPreselectionLoose(AK4JetObject const& part);
  bool testPreselectionTight(AK4JetObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


bool AK4JetSelectionHelpers::testLooseId(AK4JetObject const& part){ return part.extras.pass_looseId; }
bool AK4JetSelectionHelpers::testTightId(AK4JetObject const& part){ return part.extras.pass_tightId; }
bool AK4JetSelectionHelpers::testPUJetId(AK4JetObject const& part){ return part.extras.pass_puId; }

bool AK4JetSelectionHelpers::testLooseKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool AK4JetSelectionHelpers::testTightKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
}

bool AK4JetSelectionHelpers::testPtEtaGen(AK4JetObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
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
    part.testSelectionBit(kPUJetId)
    );
}
void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  part.setSelectionBit(kPUJetId, testPUJetId(part));

  part.setSelectionBit(kLooseId, testLooseId(part));
  part.setSelectionBit(kLooseKin, testLooseKin(part));

  part.setSelectionBit(kTightId, testTightId(part));
  part.setSelectionBit(kTightKin, testTightKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
