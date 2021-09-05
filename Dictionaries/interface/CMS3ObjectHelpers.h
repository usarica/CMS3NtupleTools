#ifndef CMS3OBJECTHELPERS_H
#define CMS3OBJECTHELPERS_H

#include <DataFormats/Common/interface/View.h>


namespace ParticleObjectHelpers{
  enum ObjectType{
    kMuon = 0,
    kElectron,
    kPhoton,
    kAK4Jet,
    kAK8Jet,

    kGenPrunedParticle,
    kGenPackedParticle,

    kGenAK4Jet,
    kGenAK8Jet,

    nObjectTypes
  };

  // edm::View<obj_type> holds const obj_type references, whereas std::vector<const obj_type*> holds pointers...
  template<typename T> void getObjectPointer(typename edm::View<T>::const_iterator const& it, T const*& ptr){ ptr = &(*it); }

}


#include "ParticleObjectHelpers.h"


#endif
