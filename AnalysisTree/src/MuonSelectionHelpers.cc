#include <cassert>

#include <DataFormats/MuonReco/interface/Muon.h>

#include "MuonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


// These are functions hidden from the user
namespace MuonSelectionHelpers{
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

  bool testPreselectionVeto(MuonObject const& part);
  bool testPreselectionLoose(MuonObject const& part);
  bool testPreselectionTight(MuonObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;
using namespace reco;


float MuonSelectionHelpers::absMiniIso(MuonObject const& part){ return part.extras.miniIso_comb_nofsr; }
float MuonSelectionHelpers::relMiniIso(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absMiniIso(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_DR0p3(MuonObject const& part){ return part.extras.pfIso03_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p3(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p3(part)/pt : 0.f); }

float MuonSelectionHelpers::absPFIso_DR0p4(MuonObject const& part){ return part.extras.pfIso04_comb_nofsr; }
float MuonSelectionHelpers::relPFIso_DR0p4(MuonObject const& part){ float pt = part.pt(); return (pt>0. ? absPFIso_DR0p4(part)/pt : 0.f); }

float MuonSelectionHelpers::getIsolationDRmax(MuonObject const& part){
  if (isoType_preselection == kPFIsoDR0p3) return 0.3;
  else if (isoType_preselection == kPFIsoDR0p4) return 0.4;
  else if (isoType_preselection == kMiniIso) return (10. / std::min(std::max(part.pt()/part.currentSystScale, 50.), 200.));
  else{
    MELAerr << "MuonSelectionHelpers::getIsolationDRmax: Isolation type " << isoType_preselection << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

float MuonSelectionHelpers::computeIso(MuonObject const& part){
  if (isoType_preselection == kPFIsoDR0p3) return relPFIso_DR0p3(part);
  else if (isoType_preselection == kPFIsoDR0p4) return relPFIso_DR0p4(part);
  else if (isoType_preselection == kMiniIso) return relMiniIso(part);
  else MELAerr << "MuonSelectionHelpers::computeIso: Isolation " << isoType_preselection << " with id " << idType_preselection << " is not implemented." << endl;
  return 999.f;
}

bool MuonSelectionHelpers::testMuonSystemTime(MuonObject const& part){
  // Test precomputed timing flag first
  if (part.extras.pass_muon_timing) return true;
  // Cut suggestions from Piotr for out-of-time muons from https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf
  // reco::Muon::InTimeMuon selector bit flag also stores the same info
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
  bool rpcok = (rpcndof>=2 && rpcerr==0.);
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
    part.testSelectionBit(bit_preselectionVeto_id)
    &&
    part.testSelectionBit(bit_preselectionVeto_iso)
    &&
    part.testSelectionBit(bit_preselectionVeto_kin)
    );
}
bool MuonSelectionHelpers::testPreselectionLoose(MuonObject const& part){
  return (
    part.testSelectionBit(bit_preselectionLoose_id)
    &&
    part.testSelectionBit(bit_preselectionLoose_iso)
    &&
    part.testSelectionBit(bit_preselectionLoose_kin)
    &&
    (bit_preselection_time != kValidMuonSystemTime || part.testSelectionBit(bit_preselection_time))
    );
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

  // The functions below test the bits set in the steps above.
  part.setSelectionBit(kPreselectionVeto, testPreselectionVeto(part));
  part.setSelectionBit(kPreselectionLoose, testPreselectionLoose(part));
  part.setSelectionBit(kPreselectionTight, testPreselectionTight(part));
}
