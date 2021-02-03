#ifndef GENJETOBJECT_H
#define GENJETOBJECT_H

#include "ParticleObject.h"


#define GENJET_VARIABLES \
GENJET_VARIABLE(float, pt, 0) \
GENJET_VARIABLE(float, eta, 0) \
GENJET_VARIABLE(float, phi, 0) \
GENJET_VARIABLE(float, mass, 0)


typedef ParticleObject GenJetObject;


#endif
