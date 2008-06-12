#include "CMS2/NtupleMaker/interface/MCUtilities.h"

MCUtilities::MCUtilities() {
}

MCUtilities::~MCUtilities() {
}

const reco::GenParticle* MCUtilities::motherID(const reco::GenParticle& gp) {
  //extract motherid from a GenParticle by walking on the left side of the mother graph
  //inlining it here to avoid cyclic depce
    
  //yuck; this is all because mother access in the candidate is not virtual
  const reco::GenParticle* mom = &gp;
  while(mom->numberOfMothers()>0) {
    for(uint j=0; j<mom->numberOfMothers(); ++j) {
      mom = dynamic_cast<const reco::GenParticle*>(mom->mother(j));
      if(mom->pdgId()!=gp.pdgId()){
	return mom;
      }
    }
  }
  return mom;
}
