#include "ParticleSelectionHelpers.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "HelperFunctions.h"


void ParticleSelectionHelpers::setUseProbeLeptonsInLooseSelection(bool flag){
  MuonSelectionHelpers::setAllowProbeIdInLooseSelection(flag);
  ElectronSelectionHelpers::setAllowProbeIdInLooseSelection(flag);
}
void ParticleSelectionHelpers::setUseFakeableLeptonsInLooseSelection(bool flag){
  MuonSelectionHelpers::setAllowFakeableInLooseSelection(flag);
  ElectronSelectionHelpers::setAllowFakeableInLooseSelection(flag);
}


// Veto, loose and tight particle ids
#define SELECTION_TYPES \
SELECTION_TYPE(Veto) \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Tight)
#define SELECTION_TYPE(TYPE) \
template<> bool ParticleSelectionHelpers::is##TYPE##Particle(MuonObject const* part){ \
  return part->testSelectionBit(MuonSelectionHelpers::kPreselection##TYPE); \
} \
template<> bool ParticleSelectionHelpers::is##TYPE##Particle(ElectronObject const* part){ \
  return part->testSelectionBit(ElectronSelectionHelpers::kPreselection##TYPE); \
} \
template<> bool ParticleSelectionHelpers::is##TYPE##Particle(PhotonObject const* part){ \
  return part->testSelectionBit(PhotonSelectionHelpers::kPreselection##TYPE); \
}
SELECTION_TYPES;
#undef SELECTION_TYPE
#undef SELECTION_TYPES


// Implementation of generic veto-tight ids
#define SELECTION_TYPES \
SELECTION_TYPE(Veto) \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Tight)
#define SELECTION_TYPE(TYPE) \
template<> bool ParticleSelectionHelpers::is##TYPE##Particle(ParticleObject const* part){ \
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part); \
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part); \
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part); \
  if (muon) return is##TYPE##Particle(muon); \
  else if (electron) return is##TYPE##Particle(electron); \
  else if (photon) return is##TYPE##Particle(photon); \
  else return false; \
}
SELECTION_TYPES;
#undef SELECTION_TYPE
#undef SELECTION_TYPES


// FSR suitability of muons and electrons requires checking individual bits
template<> bool ParticleSelectionHelpers::isFSRSuitable<ParticleObject>(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return isFSRSuitable(muon);
  else if (electron) return isFSRSuitable(electron);
  else return false;
}
template<> bool ParticleSelectionHelpers::isFSRSuitable<MuonObject>(MuonObject const* part){ return part->testSelectionBit(MuonSelectionHelpers::kPreselectionLoose_NoIso); }
template<> bool ParticleSelectionHelpers::isFSRSuitable<ElectronObject>(ElectronObject const* part){ return part->testSelectionBit(ElectronSelectionHelpers::kPreselectionLoose_NoIso); }


// Functions for jets
#define SELECTION_TYPES \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Tight)

#define SELECTION_TYPE(TYPE) \
template<> bool ParticleSelectionHelpers::is##TYPE##Jet(AK4JetObject const* part){ \
  return part->testSelectionBit(AK4JetSelectionHelpers::kPreselection##TYPE); \
} \
template<> bool ParticleSelectionHelpers::is##TYPE##Jet(AK8JetObject const* part){ \
  return part->testSelectionBit(AK8JetSelectionHelpers::kPreselection##TYPE); \
} \
template<> bool ParticleSelectionHelpers::is##TYPE##Jet(ParticleObject const* part){ \
  AK4JetObject const* ak4jet = dynamic_cast<AK4JetObject const*>(part); \
  AK8JetObject const* ak8jet = dynamic_cast<AK8JetObject const*>(part); \
  if (ak4jet) return is##TYPE##Jet(ak4jet); \
  else if (ak8jet) return is##TYPE##Jet(ak8jet); \
  else return false; \
}

SELECTION_TYPES;

#undef SELECTION_TYPE
#undef SELECTION_TYPES


template<> bool ParticleSelectionHelpers::isJetForHEMVeto<AK4JetObject>(AK4JetObject const* jet){
  // No PU jet id requirement, so do not use the isTightJet check.
  // No eta requirement either...
  return jet->testSelectionBit(AK4JetSelectionHelpers::bit_preselectionTight_id) && jet->pt()>=AK4JetSelectionHelpers::ptThr_skim_tight;
}
template<> bool ParticleSelectionHelpers::isJetForHEMVeto<AK8JetObject>(AK8JetObject const* jet){
  // No eta requirement.
  return jet->testSelectionBit(AK8JetSelectionHelpers::bit_preselectionTight_id) && jet->pt()>=AK8JetSelectionHelpers::ptThr_skim_HEMcheck;
}

template<> bool ParticleSelectionHelpers::isJetForPUJetIdSF<AK4JetObject>(AK4JetObject const* jet){
  return (
    jet->testSelectionBit(AK4JetSelectionHelpers::bit_preselectionTight_id)
    &&
    (!AK4JetSelectionHelpers::getApplyTightLeptonVetoIdToJetsFlag() || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
    &&
    //jet->testSelectionBit(AK4JetSelectionHelpers::kLooseKin)
    jet->testSelectionBit(AK4JetSelectionHelpers::kTightKin)
    &&
    jet->pt()<AK4JetSelectionHelpers::ptThr_PUId
    );
}

template<> bool ParticleSelectionHelpers::isJetForBtagSF<AK4JetObject>(AK4JetObject const* jet){
  return jet->testSelectionBit(AK4JetSelectionHelpers::kBtaggable_NoPUJetId);
}
