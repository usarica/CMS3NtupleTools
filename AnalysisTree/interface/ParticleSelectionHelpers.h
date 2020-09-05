#ifndef PARTICLESELECTIONHELPERS_H
#define PARTICLESELECTIONHELPERS_H

#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"


namespace ParticleSelectionHelpers{
  void setUseProbeLeptonsInLooseSelection(bool flag); // Allows probe leptons to be usable in jet cleaning or trigger checking
  void setUseFakeableLeptonsInLooseSelection(bool flag); // Allows fakeable leptons to be usable in jet cleaning or trigger checking

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
  template<typename T> bool isParticleForIsotrackCleaning(T const* part);
  template<typename T> bool isParticleForTriggerChecking(T const* part);

  template<typename T> bool isFSRSuitable(T const* part);
  template<> bool isFSRSuitable<ParticleObject>(ParticleObject const* part);
  template<> bool isFSRSuitable<MuonObject>(MuonObject const* part);
  template<> bool isFSRSuitable<ElectronObject>(ElectronObject const* part);


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

  template<typename T> bool isJetForTriggerChecking(T const*);
  template<typename T> bool isJetForHEMVeto(T const*);
  template<> bool isJetForHEMVeto<AK4JetObject>(AK4JetObject const*);
  template<> bool isJetForHEMVeto<AK8JetObject>(AK8JetObject const*);
  template<typename T> bool isJetForBtagSF(T const*);
  template<> bool isJetForBtagSF<AK4JetObject>(AK4JetObject const*);

}

template<typename T> bool ParticleSelectionHelpers::isParticleForJetCleaning(T const* part){ return isLooseParticle(part); }
template bool ParticleSelectionHelpers::isParticleForJetCleaning<MuonObject>(MuonObject const*);
template bool ParticleSelectionHelpers::isParticleForJetCleaning<ElectronObject>(ElectronObject const*);
template bool ParticleSelectionHelpers::isParticleForJetCleaning<PhotonObject>(PhotonObject const*);

template<typename T> bool ParticleSelectionHelpers::isParticleForIsotrackCleaning(T const* part){ return isLooseParticle(part); }
template bool ParticleSelectionHelpers::isParticleForIsotrackCleaning<MuonObject>(MuonObject const*);
template bool ParticleSelectionHelpers::isParticleForIsotrackCleaning<ElectronObject>(ElectronObject const*);
template bool ParticleSelectionHelpers::isParticleForIsotrackCleaning<PhotonObject>(PhotonObject const*);

template<typename T> bool ParticleSelectionHelpers::isParticleForTriggerChecking(T const* part){ return isLooseParticle(part); }
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<MuonObject>(MuonObject const*);
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<ElectronObject>(ElectronObject const*);
template bool ParticleSelectionHelpers::isParticleForTriggerChecking<PhotonObject>(PhotonObject const*);

template<typename T> bool ParticleSelectionHelpers::isJetForTriggerChecking(T const* jet){ return isTightJet(jet); }
template bool ParticleSelectionHelpers::isJetForTriggerChecking<AK4JetObject>(AK4JetObject const*);
template bool ParticleSelectionHelpers::isJetForTriggerChecking<AK8JetObject>(AK8JetObject const*);

#endif
