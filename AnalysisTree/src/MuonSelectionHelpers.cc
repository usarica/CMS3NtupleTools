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
  bool rpcok =(rpcndof && rpcerr==0.);
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

#define ISO_FCN MuonSelectionHelpers::relPFIso_DR0p3

bool MuonSelectionHelpers::testVetoIdIso(MuonObject const& part){
  return (testVetoId(part) && ISO_FCN(part)<isoThr_veto && testMuonSystemTime(part));
}
bool MuonSelectionHelpers::testLooseIdIso(MuonObject const& part){
  return (testLooseId(part) && ISO_FCN(part)<isoThr_loose && testMuonSystemTime(part));
}
bool MuonSelectionHelpers::testMediumIdIso(MuonObject const& part){
  return (testMediumId(part) && ISO_FCN(part)<isoThr_medium && testMuonSystemTime(part));
}
bool MuonSelectionHelpers::testTightIdIso(MuonObject const& part){
  return (testTightId(part) && ISO_FCN(part)<isoThr_tight && testMuonSystemTime(part));
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
bool MuonSelectionHelpers::testPtEtaSkim(MuonObject const& part){
  // pT and eta skim cut
  if (testTightIdIso(part)) return (part.pt()>=ptThr_skim_tight && fabs(part.eta())<etaThr_skim_tight);
  else if (testMediumIdIso(part)) return (part.pt()>=ptThr_skim_medium && fabs(part.eta())<etaThr_skim_medium);
  else if (testLooseIdIso(part)) return (part.pt()>=ptThr_skim_loose && fabs(part.eta())<etaThr_skim_loose);
  else if (testVetoIdIso(part)) return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
  else return false;
}
bool MuonSelectionHelpers::testPreselection(MuonObject const& part){
  return (
    (
    (bit_preselection_idisoreco == kVetoIDIso && testVetoIdIso(part))
      ||
      (bit_preselection_idisoreco == kLooseIDIso && testLooseIdIso(part))
      ||
      (bit_preselection_idisoreco == kMediumIDIso && testMediumIdIso(part))
      ||
      (bit_preselection_idisoreco == kTightIDIso && testTightIdIso(part))
      ) && testPtEtaSkim(part)
    );
}
void MuonSelectionHelpers::setSelectionBits(MuonObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testPtEtaGen(part)) part.setSelectionBit(kGenPtEta);

  if (testMuonSystemTime(part)) part.setSelectionBit(kValidMuonSystemTime);

  if (testVetoId(part)) part.setSelectionBit(kVetoID);
  if (testVetoIdIso(part)) part.setSelectionBit(kVetoIDIso);
  if (testLooseId(part)) part.setSelectionBit(kLooseID);
  if (testLooseIdIso(part)) part.setSelectionBit(kLooseIDIso);
  if (testMediumId(part)) part.setSelectionBit(kMediumID);
  if (testMediumIdIso(part)) part.setSelectionBit(kMediumIDIso);
  if (testTightId(part)) part.setSelectionBit(kTightID);
  if (testTightIdIso(part)) part.setSelectionBit(kTightIDIso);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
  if (testPreselection(part)) part.setSelectionBit(kPreselection);
}
