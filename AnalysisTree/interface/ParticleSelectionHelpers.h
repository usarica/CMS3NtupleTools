#ifndef PARTICLESELECTIONHELPERS_H
#define PARTICLESELECTIONHELPERS_H

#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"


namespace ParticleSelectionHelpers{
#define SELECTION_TYPES \
SELECTION_TYPE(Veto) \
SELECTION_TYPE(Loose) \
SELECTION_TYPE(Medium) \
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

}


#endif
