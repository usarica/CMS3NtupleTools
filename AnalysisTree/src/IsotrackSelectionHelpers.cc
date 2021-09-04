#include <cassert>
#include "IsotrackSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


// These are functions hidden from the user
namespace IsotrackSelectionHelpers{
  bool testVetoId(IsotrackObject const& part);
  bool testVetoIso(IsotrackObject const& part);
  bool testVetoKin(IsotrackObject const& part);

  bool testPreselectionVeto(IsotrackObject const& part);
}


using namespace std;
using namespace IvyStreamHelpers;


float IsotrackSelectionHelpers::getIsolationDRmax(IsotrackObject const& part){
  switch (isoType_preselection){
  case kPFIsoChargedDR0p3:
  case kPFIsoCombDR0p3:
    return 0.3;
  case kMiniIsoCharged:
  case kMiniIsoComb:
    return (10. / std::min(std::max(part.pt(), 50.), 200.));
  default:
    IVYerr << "IsotrackSelectionHelpers::getIsolationDRmax: Isolation type " << isoType_preselection << " is not implemented." << endl;
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

float IsotrackSelectionHelpers::computeAbsIso(IsotrackObject const& part){
  if (isoType_preselection == kPFIsoCombDR0p3) return absPFIsoComb_DR0p3(part);
  else if (isoType_preselection == kMiniIsoComb) return absMiniIsoComb(part);
  else if (isoType_preselection == kPFIsoChargedDR0p3) return absPFIsoCharged_DR0p3(part);
  else if (isoType_preselection == kMiniIsoCharged) return absMiniIsoCharged(part);
  else IVYerr << "IsotrackSelectionHelpers::computeAbsIso: Isolation " << isoType_preselection << " is not implemented." << endl;
  return 999.f;
}

float IsotrackSelectionHelpers::computeRelIso(IsotrackObject const& part){
  if (isoType_preselection == kPFIsoCombDR0p3) return relPFIsoComb_DR0p3(part);
  else if (isoType_preselection == kMiniIsoComb) return relMiniIsoComb(part);
  else if (isoType_preselection == kPFIsoChargedDR0p3) return relPFIsoCharged_DR0p3(part);
  else if (isoType_preselection == kMiniIsoCharged) return relMiniIsoCharged(part);
  else IVYerr << "IsotrackSelectionHelpers::computeRelIso: Isolation " << isoType_preselection << " is not implemented." << endl;
  return 999.f;
}

bool IsotrackSelectionHelpers::testVetoId(IsotrackObject const& part){
  auto const& extras = part.extras;
  cms3_absid_t const abs_id = std::abs(part.pdgId());
  bool const isMuon = abs_id==13;
  bool const isElectron = abs_id==11;
  bool const isChargedHadron = abs_id>100;
  return ((isMuon || isChargedHadron || isElectron) && extras.fromPV && extras.is_pfCand && std::abs(extras.dz)<0.1);
}

bool IsotrackSelectionHelpers::testVetoIso(IsotrackObject const& part){
  cms3_absid_t const abs_id = std::abs(part.pdgId());
  float const pt = part.pt();
  float const absIso = computeAbsIso(part);

  bool const isMuon = abs_id==13;
  bool const isElectron = abs_id==11;
  bool const isChargedHadron = abs_id>100;
  float const relIso_thr = (isMuon ? relIsoThr_muon_veto : (isChargedHadron ? relIsoThr_hadron_veto : (isElectron ? relIsoThr_electron_veto : 0.f)));
  float const absIso_thr = (isMuon ? absIsoThr_muon_veto : (isChargedHadron ? absIsoThr_hadron_veto : (isElectron ? absIsoThr_electron_veto : 0.f)));
  float const iso_thr = std::min(absIso_thr, pt*relIso_thr);
  return (absIso < iso_thr);
}

bool IsotrackSelectionHelpers::testVetoKin(IsotrackObject const& part){
  cms3_absid_t const abs_id = std::abs(part.pdgId());
  bool const isMuon = abs_id==13;
  bool const isElectron = abs_id==11;
  bool const isChargedHadron = abs_id>100;
  if (isMuon) return (part.pt()>=ptThr_muon_veto && std::abs(part.eta())<etaThr_muon_veto);
  else if (isChargedHadron) return (part.pt()>=ptThr_hadron_veto && std::abs(part.eta())<etaThr_hadron_veto);
  else if (isElectron) return (part.pt()>=ptThr_electron_veto && std::abs(part.eta())<etaThr_electron_veto);
  else return false;
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
  static_assert(std::numeric_limits<ParticleObject::SelectionBitsType_t>::digits >= nSelectionBits);

  part.setSelectionBit(kVetoId, testVetoId(part));
  part.setSelectionBit(kVetoIso, testVetoIso(part));
  part.setSelectionBit(kVetoKin, testVetoKin(part));

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionVeto, testPreselectionVeto(part));
}
