#ifndef CMS3_MCUTILITIES_H
#define CMS3_MCUTILITIES_H

#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"


namespace MCUtilities{

  const reco::GenParticle* motherID(const reco::GenParticle&); // Find the GenParticle mother of a GenParticle
  const reco::GenParticle* motherIDPacked(const pat::PackedGenParticle&); // Find the GenParticle mother of a PackedGenParticle. Return "0" if nothing is found.

  void getAllMothers(const reco::GenParticle*, std::vector<const reco::GenParticle*>&);
  void getAllMothers(const pat::PackedGenParticle*, std::vector<const reco::GenParticle*>&);

}


#endif
