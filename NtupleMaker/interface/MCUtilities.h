#ifndef CMS3_MCUTILITIES_H
#define CMS3_MCUTILITIES_H

#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "IvyFramework/IvyDataTools/interface/HelperFunctions.h"


namespace MCUtilities{

  reco::GenParticle const* motherID(reco::GenParticle const&); // Find the GenParticle mother of a GenParticle
  reco::GenParticle const* motherIDPacked(pat::PackedGenParticle const&); // Find the GenParticle mother of a PackedGenParticle. Return "0" if nothing is found.

  template<typename T> void getAllMothers(T const*, std::vector<reco::GenParticle const*>&, bool ignoreFSR=true);
  template<typename T> void getAllDaughters(T const*, std::vector<reco::GenParticle const*>&, bool followFSR=true);

}

template<typename T> void MCUtilities::getAllMothers(T const* part, std::vector<reco::GenParticle const*>& res, bool ignoreFSR){
  if (!part) return;
  for (size_t j=0; j<part->numberOfMothers(); j++){
    reco::GenParticle const* mom = dynamic_cast<reco::GenParticle const*>(part->mother(j));
    if (!mom) continue;
    // Walk back the tree to get a mother that is not Pythia junk; otherwise add the mother to the collection
    bool doSkip = false;
    if (mom->pdgId()==part->pdgId()){
      if (ignoreFSR || !mom->isLastCopyBeforeFSR()){ getAllMothers(mom, res); doSkip = true; }
    }
    // Check if the mother is already there
    if (!doSkip && !HelperFunctions::checkListVariable(res, mom)) res.push_back(mom);
  }
}
template<typename T> void MCUtilities::getAllDaughters(T const* part, std::vector<reco::GenParticle const*>& res, bool followFSR){
  if (!part) return;
  for (size_t j=0; j<part->numberOfDaughters(); j++){
    reco::GenParticle const* dau = dynamic_cast<reco::GenParticle const*>(part->daughter(j));
    if (!dau) continue;
    // Walk the tree to get a daughter that is not Pythia junk; otherwise add the daughter to the collection
    bool doSkip = false;
    if (dau->pdgId()==part->pdgId()){
      if (followFSR){ getAllDaughters(dau, res); doSkip = true; }
    }
    // Check if the daughter is already there
    if (!doSkip && HelperFunctions::checkListVariable(res, dau)) doSkip = true;
    if (doSkip) continue;
    res.push_back(dau);
  }
}


#endif
