// -*- C++ -*-
//
// Package:    MatchUtilities
// Class:      MatchUtilities
// 
/**\class MatchUtilities MatchUtilities.h CMS2/NtupleMaker/interface/MatchUtilities.h

Description: utilities to match objects

*/
//
// Original Author:  Oliver Gutsche
//         Wed Jun 11 17:20:33 CDT 2008
// $Id: MatchUtilities.h,v 1.1 2008/06/12 00:12:20 gutsche Exp $
//
//
#ifndef CMS2_MATCHUTILITIES_H
#define CMS2_MATCHUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"


class MatchUtilities {
public:
  MatchUtilities();
  ~MatchUtilities();

  static const reco::GenParticle* matchCandToGen(const reco::Candidate&, const std::vector<reco::GenParticle>* getParticles);
  static const reco::GenJet* matchCandToGenJet(const reco::Candidate& jet,  const std::vector<reco::GenJet>* getJet); 
  
  static const reco::Candidate* matchGenToCand(const reco::GenParticle&, std::vector<const reco::Candidate*> cand); 
  static const reco::Candidate* matchGenToCand(const reco::GenJet&, std::vector<const reco::Candidate*> cand);

private:

};

#endif
