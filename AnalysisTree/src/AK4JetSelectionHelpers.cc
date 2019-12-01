#include <cassert>
#include "AK4JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


bool AK4JetSelectionHelpers::testLooseId(AK4JetObject const& part){
  return part.extras.pass_looseId;
}
bool AK4JetSelectionHelpers::testTightId(AK4JetObject const& part){
  return part.extras.pass_tightId;
}
bool AK4JetSelectionHelpers::testPUJetId(AK4JetObject const& part){
  return part.extras.pass_puId;
}

bool AK4JetSelectionHelpers::testLooseKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool AK4JetSelectionHelpers::testTightKin(AK4JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
}

bool AK4JetSelectionHelpers::testPtEtaGen(AK4JetObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool AK4JetSelectionHelpers::testPreselection(AK4JetObject const& part){
  return (
    (
      (bit_preselection_id == kLooseId && testLooseId(part))
      ||
      (bit_preselection_id == kTightId && testTightId(part))
      )
    &&
    (
      (bit_preselection_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselection_kin == kTightKin && testTightKin(part))
      )
    &&
    testPUJetId(part)
    );
}
void AK4JetSelectionHelpers::setSelectionBits(AK4JetObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);

  if (testPUJetId(part)) part.setSelectionBit(kPUJetId);

  if (testLooseId(part)) part.setSelectionBit(kLooseId);
  if (testLooseKin(part)) part.setSelectionBit(kLooseKin);

  if (testTightId(part)) part.setSelectionBit(kTightId);
  if (testTightKin(part)) part.setSelectionBit(kTightKin);

  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
