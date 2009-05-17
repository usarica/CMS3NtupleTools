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
// $Id: MatchUtilities.h,v 1.6 2009/05/17 19:37:41 kalavase Exp $
//
//
#ifndef CMS2_MATCHUTILITIES_H
#define CMS2_MATCHUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <Math/VectorUtil.h>
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"



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

  static const void alignRecoPatJetCollections(const std::vector<reco::CaloJet>&,
					       std::vector<pat::Jet>&);
  static const void alignRecoPatElectronCollections(const std::vector<reco::GsfElectron>&,
					     std::vector<pat::Electron>&);
  static const void alignRecoPatMuonCollections(const std::vector<reco::Muon>&,
					     std::vector<pat::Muon>&);
  
  template <class T1, class T2> static const void alignCollections(const std::vector<T1>& v_ref,
							    std::vector<T2>& v_toAllign) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
    
    if(v_ref.size() != v_toAllign.size())
      throw cms::Exception("MatchUtilities") 
	<< "The two collections you're trying to allign do not have the same number of entries!!!" 
	<< " Exiting. It's probably the PAT's fault. Trust me.";
  
    std::vector<T2> v_temp = v_toAllign;
    v_toAllign.clear();
  

    //loop over the Reference Collection
    for(unsigned int i = 0; i < v_ref.size(); i++) {

      double dR = 0.05;
      unsigned int matchIdx = 9999;
      LorentzVector ref_p4 = v_ref.at(i).p4();

      //now loop over the collection to allign
      for(unsigned int j = 0; j < v_temp.size(); j++) {
	LorentzVector temp_p4 = v_temp.at(j).p4();
	double newdR = ROOT::Math::VectorUtil::DeltaR(ref_p4, temp_p4);
      
	if(newdR < dR) {
	  dR = newdR;
	  matchIdx = j;
	}
      }//vector to allign loop
      if(matchIdx == 9999) {
	std::cout << "Object not found in collection!!!!" << std::endl;
      }
      v_toAllign.push_back(v_temp.at(matchIdx) );
    }//loop over the reference collection
    
  }



};

#endif
