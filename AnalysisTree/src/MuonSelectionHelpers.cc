#include <cassert>

#include <DataFormats/MuonReco/interface/Muon.h>

#include "MuonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace reco;


float MuonSelectionHelpers::absMiniIso_DR0p3(MuonObject const& part){ return part.extras.miniIso_comb_nofsr; }
float MuonSelectionHelpers::relMiniIso_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIso_DR0p3(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_DR0p3(MuonObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_DR0p4(MuonObject const& part){ return part.extras.pfIso04_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p4(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p4(part)/pt : 0.f); }

float MuonSelectionHelpers::computeIso(MuonObject const& part){
  if (isoType_preselection == kPFIsoDR0p3) return relPFIso_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p4) return relPFIso_DR0p4(part);
  else if (isoType_preselection == kMiniIsoDR0p3) return relMiniIso_DR0p3(part);
  else MELAerr << "MuonSelectionHelpers::computeIso: Isolation " << isoType_preselection << " with id " << idType_preselection << " is not implemented." << endl;
  return 999.f;
}

bool MuonSelectionHelpers::testMuonSystemTime(MuonObject const& part){
  // Cut suggestions from Piotr for out-of-time muons from https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf:
  float const& cmb = part.extras.time_comb_IPInOut;
  float const& rpc = part.extras.time_rpc_IPInOut;
  //float const& cmberr = part.extras.time_comb_IPInOutError;
  float const& rpcerr = part.extras.time_rpc_IPInOutError;
  int const& cmbndof = part.extras.time_comb_ndof;
  int const& rpcndof = part.extras.time_rpc_ndof;
  bool cmbok = (cmbndof>7);
  // RPC timing stored is the average over all RPC hits
  // The measurements are in multiples of the bunch crossing time since only the bunch crossing id is measured.
  // nDof>=2 ensures at least two measurements, and time error = 0 ensures measurement at the SAME BX!
  bool rpcok =(rpcndof>=2 && rpcerr==0.);
  if (rpcok){
    if ((std::abs(rpc)>10.) && !(cmbok && std::abs(cmb)<10.)) return false;
  }
  else{
    if (cmbok && (cmb>20. || cmb<-45.)) return false;
  }
  return true;
}

#define ID_CUTBASED_VETO CutBasedIdLoose
#define ID_CUTBASED_LOOSE CutBasedIdLoose
#define ID_CUTBASED_MEDIUM CutBasedIdMediumPrompt
#define ID_CUTBASED_TIGHT CutBasedIdTight
bool MuonSelectionHelpers::testVetoId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & Muon::ID_CUTBASED_VETO) == Muon::ID_CUTBASED_VETO);
  default:
    MELAerr << "MuonSelectionHelpers::testVetoId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testLooseId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & Muon::ID_CUTBASED_LOOSE) == Muon::ID_CUTBASED_LOOSE);
  default:
    MELAerr << "MuonSelectionHelpers::testLooseId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testMediumId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & Muon::ID_CUTBASED_MEDIUM) == Muon::ID_CUTBASED_MEDIUM);
  default:
    MELAerr << "MuonSelectionHelpers::testMediumId: Id " << idType_preselection << " is not implemented!" << endl;
    assert(0);
    return false;
  };
}
bool MuonSelectionHelpers::testTightId(MuonObject const& part){
  switch (idType_preselection){
  case kCutBasedId_MuonPOG:
    return ((part.extras.POG_selector_bits & Muon::ID_CUTBASED_TIGHT) == Muon::ID_CUTBASED_TIGHT);
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

bool MuonSelectionHelpers::testVetoKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
}
bool MuonSelectionHelpers::testLooseKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
}
bool MuonSelectionHelpers::testMediumKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
}
bool MuonSelectionHelpers::testTightKin(MuonObject const& part){
  return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
}
bool MuonSelectionHelpers::testSoftKin(MuonObject const& part){ return (part.pt()>=ptThr_skim_soft && fabs(part.eta())<etaThr_skim_soft); }

bool MuonSelectionHelpers::testPtEtaGen(MuonObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool MuonSelectionHelpers::testPreselectionVeto(MuonObject const& part){
  return (
    (
    (bit_preselectionVeto_iso == kVetoIso && testVetoIso(part))
      ||
      (bit_preselectionVeto_iso == kLooseIso && testLooseIso(part))
      ||
      (bit_preselectionVeto_iso == kMediumIso && testMediumIso(part))
      ||
      (bit_preselectionVeto_iso == kTightIso && testTightIso(part))
      )
    &&
    (
    (bit_preselectionVeto_id == kVetoId && testVetoId(part))
      ||
      (bit_preselectionVeto_id == kLooseId && testLooseId(part))
      ||
      (bit_preselectionVeto_id == kMediumId && testMediumId(part))
      ||
      (bit_preselectionVeto_id == kTightId && testTightId(part))
      )
    &&
    (
    (bit_preselectionVeto_kin == kVetoKin && testVetoKin(part))
      ||
      (bit_preselectionVeto_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselectionVeto_kin == kMediumKin && testMediumKin(part))
      ||
      (bit_preselectionVeto_kin == kTightKin && testTightKin(part))
      )
    );
}
bool MuonSelectionHelpers::testPreselectionLoose(MuonObject const& part){
  return (
    (
    (bit_preselectionLoose_iso == kVetoIso && testVetoIso(part))
      ||
      (bit_preselectionLoose_iso == kLooseIso && testLooseIso(part))
      ||
      (bit_preselectionLoose_iso == kMediumIso && testMediumIso(part))
      ||
      (bit_preselectionLoose_iso == kTightIso && testTightIso(part))
      )
    &&
    (
    (bit_preselectionLoose_id == kVetoId && testVetoId(part))
      ||
      (bit_preselectionLoose_id == kLooseId && testLooseId(part))
      ||
      (bit_preselectionLoose_id == kMediumId && testMediumId(part))
      ||
      (bit_preselectionLoose_id == kTightId && testTightId(part))
      )
    &&
    (
    (bit_preselectionLoose_kin == kVetoKin && testVetoKin(part))
      ||
      (bit_preselectionLoose_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselectionLoose_kin == kMediumKin && testMediumKin(part))
      ||
      (bit_preselectionLoose_kin == kTightKin && testTightKin(part))
      )
    );
}
bool MuonSelectionHelpers::testPreselectionAccept(MuonObject const& part){
  return (
    (
    (bit_preselectionAccept_iso == kVetoIso && testVetoIso(part))
      ||
      (bit_preselectionAccept_iso == kLooseIso && testLooseIso(part))
      ||
      (bit_preselectionAccept_iso == kMediumIso && testMediumIso(part))
      ||
      (bit_preselectionAccept_iso == kTightIso && testTightIso(part))
      )
    &&
    (
    (bit_preselectionAccept_id == kVetoId && testVetoId(part))
      ||
      (bit_preselectionAccept_id == kLooseId && testLooseId(part))
      ||
      (bit_preselectionAccept_id == kMediumId && testMediumId(part))
      ||
      (bit_preselectionAccept_id == kTightId && testTightId(part))
      )
    &&
    (
    (bit_preselectionAccept_kin == kVetoKin && testVetoKin(part))
      ||
      (bit_preselectionAccept_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselectionAccept_kin == kMediumKin && testMediumKin(part))
      ||
      (bit_preselectionAccept_kin == kTightKin && testTightKin(part))
      )
    );
}
void MuonSelectionHelpers::setSelectionBits(MuonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);

  if (testMuonSystemTime(part)) part.setSelectionBit(kValidMuonSystemTime);

  if (testVetoId(part)) part.setSelectionBit(kVetoId);
  if (testVetoIso(part)) part.setSelectionBit(kVetoIso);
  if (testVetoKin(part)) part.setSelectionBit(kVetoKin);

  if (testLooseId(part)) part.setSelectionBit(kLooseId);
  if (testLooseIso(part)) part.setSelectionBit(kLooseIso);
  if (testLooseKin(part)) part.setSelectionBit(kLooseKin);

  if (testMediumId(part)) part.setSelectionBit(kMediumId);
  if (testMediumIso(part)) part.setSelectionBit(kMediumIso);
  if (testMediumKin(part)) part.setSelectionBit(kMediumKin);

  if (testTightId(part)) part.setSelectionBit(kTightId);
  if (testTightIso(part)) part.setSelectionBit(kTightIso);
  if (testTightKin(part)) part.setSelectionBit(kTightKin);

  if (testSoftId(part)) part.setSelectionBit(kSoftId);
  if (testSoftIso(part)) part.setSelectionBit(kSoftIso);
  if (testSoftKin(part)) part.setSelectionBit(kSoftKin);

  if (testPreselectionVeto(part)) part.setSelectionBit(kPreselectionVeto);
  if (testPreselectionLoose(part)) part.setSelectionBit(kPreselectionLoose);
  if (testPreselectionAccept(part)) part.setSelectionBit(kPreselectionAccept);
}
