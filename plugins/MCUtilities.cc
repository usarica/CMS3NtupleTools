#include <algorithm>
#include <vector>
#include <iterator>

#include "CMS3/NtupleMaker/interface/plugins/MCUtilities.h"


using namespace std;
using namespace reco;


const reco::GenParticle* MCUtilities::motherID(const reco::GenParticle& gp){
  //extract motherid from a GenParticle by walking on the left side of the mother graph
  //inlining it here to avoid cyclic depce
  //yuck; this is all because mother access in the candidate is not virtual
  const reco::GenParticle* mom = &gp;
  while (mom->numberOfMothers() > 0) {
    for (uint j = 0; j < mom->numberOfMothers(); ++j) {
      mom = dynamic_cast<const reco::GenParticle*>(mom->mother(j));
      if (mom->pdgId()!=gp.pdgId()) return mom;
    }
  }

  return mom;
}
const reco::GenParticle* MCUtilities::motherIDPacked(const pat::PackedGenParticle& gp) {
  const pat::PackedGenParticle* momPGP = &gp;
  const reco::GenParticle* firstMomGP = (const reco::GenParticle*) momPGP->mother(0); // link to the prunedGenParticles collection
  if (firstMomGP){
    if (firstMomGP->pdgId() != momPGP->pdgId()) return firstMomGP;
    const reco::GenParticle* momGP = MCUtilities::motherID(*firstMomGP); // then call the usual function
    return momGP;
  }
  else return nullptr;
}

void MCUtilities::getAllMothers(const reco::GenParticle* part, std::vector<const reco::GenParticle*>& res){
  if (!part) return;
  for (size_t j=0; j<part->numberOfMothers(); j++){
    const reco::GenParticle* mom = dynamic_cast<const reco::GenParticle*>(part->mother(j));
    if (!mom) continue;
    // Walk back the tree to get a mother that is not Pythia junk; otherwise add the mother to the collection
    if (mom->pdgId()==part->pdgId()) getAllMothers(mom, res);
    else{
      // Check if the mother is already there
      auto it = std::find(std::begin(res), std::end(res), mom);
      if (it!=std::end(res)) res.push_back(mom);
    }
  }
}
void MCUtilities::getAllMothers(const pat::PackedGenParticle* part, std::vector<const reco::GenParticle*>& res){
  if (!part) return;
  for (size_t j=0; j<part->numberOfMothers(); j++){
    const reco::GenParticle* mom = dynamic_cast<const reco::GenParticle*>(part->mother(j));
    if (!mom) continue;
    // Walk back the tree to get a mother that is not Pythia junk; otherwise add the mother to the collection
    if (mom->pdgId()==part->pdgId()) getAllMothers(mom, res);
    else{
      // Check if the mother is already there
      auto it = std::find(std::begin(res), std::end(res), mom);
      if (it!=std::end(res)) res.push_back(mom);
    }
  }
}
