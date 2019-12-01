#include <cassert>
#include "AK8JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


bool AK8JetSelectionHelpers::testLooseId(AK8JetObject const& part){
  return true;
}
bool AK8JetSelectionHelpers::testTightId(AK8JetObject const& part){
  return true;
}

bool AK8JetSelectionHelpers::testLooseKin(AK8JetObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool AK8JetSelectionHelpers::testTightKin(AK8JetObject const& part){
  return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
}

bool AK8JetSelectionHelpers::testPtEtaGen(AK8JetObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool AK8JetSelectionHelpers::testPreselection(AK8JetObject const& part){
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
    );
}
void AK8JetSelectionHelpers::setSelectionBits(AK8JetObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);

  if (testLooseId(part)) part.setSelectionBit(kLooseId);
  if (testLooseKin(part)) part.setSelectionBit(kLooseKin);

  if (testTightId(part)) part.setSelectionBit(kTightId);
  if (testTightKin(part)) part.setSelectionBit(kTightKin);

  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
