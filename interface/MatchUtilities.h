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
// Wed Jun 11 17:20:33 CDT 2008
// $Id: MatchUtilities.h,v 1.5 2008/12/09 00:23:00 kalavase Exp $
//
//
#ifndef CMS2_MATCHUTILITIES_H
#define CMS2_MATCHUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <Math/VectorUtil.h>



class MatchUtilities {
public:
  MatchUtilities();
  ~MatchUtilities();

  static const reco::GenParticle* matchCandToGen(const reco::Candidate&, const std::vector<reco::GenParticle>* genParticles);
  static const reco::GenParticle* matchCandToGen(const reco::Candidate&, const std::vector<reco::GenParticle>* genParticles,
						 int& genidx, int status);
  static const reco::GenParticle* matchCandToGen(const reco::Track&, const std::vector<reco::GenParticle>* genParticles, 
						 int& genidx, int status);
  static const reco::GenParticle* matchCandToGen(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& candp4, 
						 const std::vector<reco::GenParticle>* genParticles, int& genidx, int status);
  static const reco::GenJet* matchCandToGenJet(const reco::Candidate& jet,  const std::vector<reco::GenJet>* genJets);
  static const reco::GenJet* matchCandToGenJet(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& genJetp4, 
					       const std::vector<reco::GenJet>* genJets);
  
  static const reco::Candidate* matchGenToCand(const reco::GenParticle&, std::vector<const reco::Candidate*> cand);
  static const reco::Candidate* matchGenToCand(const reco::GenJet&, std::vector<const reco::Candidate*> cand);
  
  static const bool isStableGenPart(reco::GenParticle);
private:

};

#endif
