#ifndef PARTICLESELECTIONHELPERS_H
#define PARTICLESELECTIONHELPERS_H

#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"


namespace ParticleSelectionHelpers{
  // Functions to select "preselection" id types to particles
#define SELECTION_TYPES \
SELECTION_TYPE(Veto) \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Tight)

#define SELECTION_TYPE(TYPE) \
  template<typename T> bool is##TYPE##Particle(T const* part); \
  template<> bool is##TYPE##Particle(ParticleObject const* part); \
  template<> bool is##TYPE##Particle(MuonObject const* part); \
  template<> bool is##TYPE##Particle(ElectronObject const* part); \
  template<> bool is##TYPE##Particle(PhotonObject const* part);

  SELECTION_TYPES;

#undef SELECTION_TYPE
#undef SELECTION_TYPES

  template<typename T> bool isParticleForJetCleaning(T const* part);
  template<typename T> bool isParticleForTriggerChecking(T const* part);


  // Functions to select "preselection" id types to jets
#define SELECTION_TYPES \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Tight)

#define SELECTION_TYPE(TYPE) \
  template<typename T> bool is##TYPE##Jet(T const* part); \
  template<> bool is##TYPE##Jet(AK4JetObject const* part); \
  template<> bool is##TYPE##Jet(AK8JetObject const* part);

  SELECTION_TYPES;

#undef SELECTION_TYPE
#undef SELECTION_TYPES

  template<typename T> bool isJetForTriggerChecking(T const* part);

}

template<typename T> bool ParticleSelectionHelpers::isParticleForJetCleaning(T const* part){ return isLooseParticle(part); }
template bool ParticleSelectionHelpers::isParticleForJetCleaning<MuonObject>(MuonObject const*);
template bool ParticleSelectionHelpers::isParticleForJetCleaning<ElectronObject>(ElectronObject const*);
template bool ParticleSelectionHelpers::isParticleForJetCleaning<PhotonObject>(PhotonObject const*);

template<typename T> bool ParticleSelectionHelpers::isParticleForTriggerChecking(T const* part){ return isLooseParticle(part); }
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<MuonObject>(MuonObject const*);
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<ElectronObject>(ElectronObject const*);
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<PhotonObject>(PhotonObject const*);

template<typename T> bool ParticleSelectionHelpers::isJetForTriggerChecking(T const* part){ return isTightJet(part); }
template bool ParticleSelectionHelpers::isJetForTriggerChecking<AK4JetObject>(AK4JetObject const*);
template bool ParticleSelectionHelpers::isJetForTriggerChecking<AK8JetObject>(AK8JetObject const*);


#endif
