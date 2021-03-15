#include <cassert>
#include <cmath>

#include "GenParticleSelectionHelpers.h"
#include "PDGHelpers.h"
#include "MELAStreamHelpers.hh"


namespace GenParticleSelectionHelpers{
  bool testHardPromptFinalVisibleParticle(GenParticleObject const& part);
}


using namespace std;
using namespace MELAStreamHelpers;


bool GenParticleSelectionHelpers::testHardPromptFinalVisibleParticle(GenParticleObject const& part){
  auto const& extras = part.extras;
  cms3_id_t const& part_id = part.pdgId();
  ParticleObject::LorentzVector_t::Scalar const part_pt = part.pt();
  return (
    extras.isPromptFinalState
    && (
      (PDGHelpers::isALepton(part_id) && part_pt>=ptThr_hardparticle_lepton)
      ||
      (PDGHelpers::isAPhoton(part_id) && part_pt>=ptThr_hardparticle_photon)
      )
    );
}
void GenParticleSelectionHelpers::setSelectionBits(GenParticleObject& part){
  static_assert(std::numeric_limits<ParticleObject::SelectionBitsType_t>::digits >= nSelectionBits);

  part.setSelectionBit(kHardPromptFinalVisibleParticle, testHardPromptFinalVisibleParticle(part));
}
