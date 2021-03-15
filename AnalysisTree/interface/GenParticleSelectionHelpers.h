#ifndef GENPARTICLESELECTIONHELPERS_H
#define GENPARTICLESELECTIONHELPERS_H

#include "GenParticleObject.h"


namespace GenParticleSelectionHelpers{
  enum SelectionBits{
    // Prompt final state lepton or photon that do not come from soft radiation.
    // Excludes taus, neutrinos, quarks since they are not precisely 'visible'.
    kHardPromptFinalVisibleParticle,

    nSelectionBits
  };

  // Kinematic pT thresholds
  constexpr ParticleObject::LorentzVector_t::Scalar ptThr_hardparticle_lepton = 0.5;
  constexpr ParticleObject::LorentzVector_t::Scalar ptThr_hardparticle_photon = 5.;

  void setSelectionBits(GenParticleObject& part);

}


#endif
