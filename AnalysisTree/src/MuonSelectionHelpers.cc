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

#define ID_PASS_VETO CutBasedIdLoose
#define ID_PASS_LOOSE CutBasedIdLoose
#define ID_PASS_MEDIUM CutBasedIdMedium
#define ID_PASS_TIGHT CutBasedIdTight
#define ISO_FCN MuonSelectionHelpers::relPFIso_DR0p3

bool MuonSelectionHelpers::testVetoId(MuonObject const& part){
  return ((part.extras.POG_selector_bits & Muon::ID_PASS_VETO) == Muon::ID_PASS_VETO);
}
bool MuonSelectionHelpers::testLooseId(MuonObject const& part){
  return ((part.extras.POG_selector_bits & Muon::ID_PASS_LOOSE) == Muon::ID_PASS_LOOSE);
}
bool MuonSelectionHelpers::testMediumId(MuonObject const& part){
  return ((part.extras.POG_selector_bits & Muon::ID_PASS_MEDIUM) == Muon::ID_PASS_MEDIUM);
}
bool MuonSelectionHelpers::testTightId(MuonObject const& part){
  return ((part.extras.POG_selector_bits & Muon::ID_PASS_TIGHT) == Muon::ID_PASS_TIGHT);
}

bool MuonSelectionHelpers::testVetoIso(MuonObject const& part){
  return (ISO_FCN(part)<isoThr_veto);
}
bool MuonSelectionHelpers::testLooseIso(MuonObject const& part){
  return (ISO_FCN(part)<isoThr_loose);
}
bool MuonSelectionHelpers::testMediumIso(MuonObject const& part){
  return (ISO_FCN(part)<isoThr_medium);
}
bool MuonSelectionHelpers::testTightIso(MuonObject const& part){
  return (ISO_FCN(part)<isoThr_tight);
}

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

#ifdef ID_PASS_VETO
#undef ID_PASS_VETO
#endif
#ifdef ID_PASS_LOOSE
#undef ID_PASS_LOOSE
#endif
#ifdef ID_PASS_MEDIUM
#undef ID_PASS_MEDIUM
#endif
#ifdef ID_PASS_TIGHT
#undef ID_PASS_TIGHT
#endif
#ifdef ISO_FCN
#undef ISO_FCN
#endif

bool MuonSelectionHelpers::testPtEtaGen(MuonObject const& part){
  return (part.pt()>=ptThr_gen && fabs(part.eta())<etaThr_gen);
}
bool MuonSelectionHelpers::testPreselection(MuonObject const& part){
  return (
    (
    (bit_preselection_iso == kVetoIso && testVetoIso(part))
      ||
      (bit_preselection_iso == kLooseIso && testLooseIso(part))
      ||
      (bit_preselection_iso == kMediumIso && testMediumIso(part))
      ||
      (bit_preselection_iso == kTightIso && testTightIso(part))
      )
    &&
    (
    (bit_preselection_id == kVetoId && testVetoId(part))
      ||
      (bit_preselection_id == kLooseId && testLooseId(part))
      ||
      (bit_preselection_id == kMediumId && testMediumId(part))
      ||
      (bit_preselection_id == kTightId && testTightId(part))
      )
    &&
    (
    (bit_preselection_kin == kVetoKin && testVetoKin(part))
      ||
      (bit_preselection_kin == kLooseKin && testLooseKin(part))
      ||
      (bit_preselection_kin == kMediumKin && testMediumKin(part))
      ||
      (bit_preselection_kin == kTightKin && testTightKin(part))
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

  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
