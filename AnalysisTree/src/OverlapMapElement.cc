#include "OverlapMapElement.h"

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
    ParticleObject::LorentzVector_t res; \
    res = ParticleObject::PolarLorentzVector_t(commonPFCandidates_sump4_pt, commonPFCandidates_sump4_eta, commonPFCandidates_sump4_phi, commonPFCandidates_sump4_mass); \
  } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS


#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(ElectronObject, AK8JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(PhotonObject, AK8JetObject)
#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_commonMuCands_goodMET() const{ \
  ParticleObject::LorentzVector_t res; \
  if (this->isValid()){ \
    res = ParticleObject::LorentzVector_t(commonGoodMETPFMuons_sump4_px, commonGoodMETPFMuons_sump4_py, 0, 0); \
  } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS

#define OVERLAPMAP_SPECIALIZATIONS \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK4JetObject) \
OVERLAPMAP_SPECIALIZATION(MuonObject, AK8JetObject)
#define OVERLAPMAP_SPECIALIZATION(T1, T2) \
template<> ParticleObject::LorentzVector_t OverlapMapElement<T1, T2>::p4_commonMuCands_goodMET() const{ \
  ParticleObject::LorentzVector_t res; \
  if (linkedElementPair.first && linkedElementPair.first->extras.has_pfMuon_goodMET){ \
    res = linkedElementPair.first->uncorrected_p4(); \
  } \
  return res; \
}
OVERLAPMAP_SPECIALIZATIONS;
#undef OVERLAPMAP_SPECIALIZATION
#undef OVERLAPMAP_SPECIALIZATIONS
