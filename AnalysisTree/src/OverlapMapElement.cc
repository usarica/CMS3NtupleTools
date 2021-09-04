#include <cassert>
#include "OverlapMapElement.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK8JetObject)
#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_common() const{ \
  ParticleObject::LorentzVector_t res; \
  if (this->isValid()){ \
    res = ParticleObject::PolarLorentzVector_t(commonPFCandidates_sump4_pt, commonPFCandidates_sump4_eta, commonPFCandidates_sump4_phi, commonPFCandidates_sump4_mass); \
  } \
  else{ IVYerr << "OverlapMapElement<" << #T1 << ", " << #T2 << ">::p4_common: Map element is invalid!" << endl; assert(0); } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS


#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK4JetObject)
#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_commonMuCands_goodMET() const{ \
  ParticleObject::LorentzVector_t res; \
  if (this->isValid()){ \
    res = ParticleObject::LorentzVector_t(commonGoodMETPFMuons_sump4_px, commonGoodMETPFMuons_sump4_py, 0, 0); \
  } \
  else{ IVYerr << "OverlapMapElement<" << #T1 << ", " << #T2 << ">::p4_commonMuCands_goodMET: Map element is invalid!" << endl; assert(0); } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS

#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK4JetObject)
#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_commonMuCands_goodMET() const{ \
  ParticleObject::LorentzVector_t res; \
  if (linkedElementPair.first){ \
    if (linkedElementPair.first->extras.has_pfMuon_goodMET) res = ParticleObject::PolarLorentzVector_t(commonPFCandidates_sump4_pt, commonPFCandidates_sump4_eta, commonPFCandidates_sump4_phi, commonPFCandidates_sump4_mass); \
  } \
  else{ IVYerr << "OverlapMapElement<" << #T1 << ", " << #T2 << ">::p4_commonMuCands_goodMET: The first of the linked element pair is not linked!" << endl; assert(0); } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS
