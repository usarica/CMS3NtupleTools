#include <cassert>
#include <cmath>

#include "AK8JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


namespace AK8JetSelectionHelpers{
  bool testPtEtaGen(AK8JetObject const& part);

  bool testLooseId(AK8JetObject const& part);
  bool testLooseKin(AK8JetObject const& part);

  bool testTightId(AK8JetObject const& part);
  bool testTightKin(AK8JetObject const& part);

  bool testPreselectionLoose(AK8JetObject const& part);
  bool testPreselectionTight(AK8JetObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


bool AK8JetSelectionHelpers::testLooseId(AK8JetObject const& part){ return true; }
bool AK8JetSelectionHelpers::testTightId(AK8JetObject const& part){ return true; }

bool AK8JetSelectionHelpers::testLooseKin(AK8JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && std::abs(part.eta())<etaThr_skim_loose);
}
bool AK8JetSelectionHelpers::testTightKin(AK8JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && std::abs(part.eta())<etaThr_skim_tight);
}

bool AK8JetSelectionHelpers::testPtEtaGen(AK8JetObject const& part){
  return (part.pt()>=ptThr_gen && std::abs(part.eta())<etaThr_gen);
}
bool AK8JetSelectionHelpers::testPreselectionLoose(AK8JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    );
}
bool AK8JetSelectionHelpers::testPreselectionTight(AK8JetObject const& part){
  return (
    part.testSelectionBit(bit_preselectionTight_id)
    &&
    part.testSelectionBit(bit_preselectionTight_kin)
    );
}
void AK8JetSelectionHelpers::setSelectionBits(AK8JetObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kGenPtEta, testPtEtaGen(part));

  part.setSelectionBit(kLooseId, testLooseId(part));
  part.setSelectionBit(kLooseKin, testLooseKin(part));

  part.setSelectionBit(kTightId, testTightId(part));
  part.setSelectionBit(kTightKin, testTightKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
