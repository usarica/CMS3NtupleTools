#include <cassert>
#include "IsotrackSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


// These are functions hidden from the user
namespace IsotrackSelectionHelpers{
  bool testVetoId(IsotrackObject const& part);
  bool testVetoIso(IsotrackObject const& part);
  bool testVetoKin(IsotrackObject const& part);

  bool testPreselectionVeto(IsotrackObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;
using namespace reco;


float IsotrackSelectionHelpers::getIsolationDRmax(IsotrackObject const& part){
  switch (isoType_preselection){
  case kPFIsoChargedDR0p3:
  case kPFIsoCombDR0p3:
    return 0.3;
  case kMiniIsoCharged:
  case kMiniIsoComb:
    return (10. / std::min(std::max(part.pt(), 50.), 200.));
  default:
    MELAerr << "IsotrackSelectionHelpers::getIsolationDRmax: Isolation type " << isoType_preselection << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

float IsotrackSelectionHelpers::absMiniIsoCharged(IsotrackObject const& part){ return part.extras.miniIso_ch; }
float IsotrackSelectionHelpers::relMiniIsoCharged(IsotrackObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIsoCharged(part)/pt : 0.f); }

float IsotrackSelectionHelpers::absPFIsoCharged_DR0p3(IsotrackObject const& part){ return part.extras.pfIso03_ch; }
float IsotrackSelectionHelpers::relPFIsoCharged_DR0p3(IsotrackObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIsoCharged_DR0p3(part)/pt : 0.f); }

float IsotrackSelectionHelpers::absMiniIsoComb(IsotrackObject const& part){ return part.extras.miniIso_comb_nofsr; }
float IsotrackSelectionHelpers::relMiniIsoComb(IsotrackObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIsoComb(part)/pt : 0.f); }

float IsotrackSelectionHelpers::absPFIsoComb_DR0p3(IsotrackObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float IsotrackSelectionHelpers::relPFIsoComb_DR0p3(IsotrackObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIsoComb_DR0p3(part)/pt : 0.f); }

float IsotrackSelectionHelpers::computeIso(IsotrackObject const& part){
  if (isoType_preselection == kPFIsoCombDR0p3) return relPFIsoComb_DR0p3(part);
  else if (isoType_preselection == kMiniIsoComb) return relMiniIsoComb(part);
  else if (isoType_preselection == kPFIsoChargedDR0p3) return relPFIsoCharged_DR0p3(part);
  else if (isoType_preselection == kMiniIsoCharged) return relMiniIsoCharged(part);
  else MELAerr << "IsotrackSelectionHelpers::computeIso: Isolation " << isoType_preselection << " is not implemented." << endl;
  return 999.f;
}

bool IsotrackSelectionHelpers::testVetoId(IsotrackObject const& part){
  auto const& extras = part.extras;
  return (extras.is_pfCand && part.pdgId()!=0 && fabs(extras.dz)<0.1);
}

bool IsotrackSelectionHelpers::testVetoIso(IsotrackObject const& part){
  int const abs_id = std::abs(part.pdgId());
  bool const isLepton = abs_id<15;
  float const pt = part.pt();
  float const absIso = pt * computeIso(part);
  float const iso_thr = std::min(pt*(isLepton ? relIsoThr_lepton_veto : relIsoThr_hadron_veto), (isLepton ? absIsoThr_lepton_veto : absIsoThr_hadron_veto));
  return (absIso < iso_thr);
}

bool IsotrackSelectionHelpers::testVetoKin(IsotrackObject const& part){
  int const abs_id = std::abs(part.pdgId());
  bool const isLepton = abs_id<15;
  return (part.pt()>=(isLepton ? ptThr_lepton_veto : ptThr_hadron_veto) && std::abs(part.eta())<(isLepton ? etaThr_lepton_veto : etaThr_hadron_veto));
}

bool IsotrackSelectionHelpers::testPreselectionVeto(IsotrackObject const& part){
  return (
    part.testSelectionBit(bit_preselectionVeto_id)
    &&
    part.testSelectionBit(bit_preselectionVeto_iso)
    &&
    part.testSelectionBit(bit_preselectionVeto_kin)
    );
}
void IsotrackSelectionHelpers::setSelectionBits(IsotrackObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  part.setSelectionBit(kVetoId, testVetoId(part));
  part.setSelectionBit(kVetoIso, testVetoIso(part));
  part.setSelectionBit(kVetoKin, testVetoKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionVeto, testPreselectionVeto(part));
}
