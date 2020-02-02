#include "ParticleSelectionHelpers.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "HelperFunctions.h"


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
