// -*- C++ -*-
//
// Package:    HypTrilepMaker
// Class:      HypTrilepMaker
// 
/**\class HypTrilepMaker HypTrilepMaker.h CMS2/NtupleMaker/interface/HypTrilepMaker.h

Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypTrilepMaker.h,v 1.1 2008/06/24 00:39:34 gutsche Exp $
//
//
#ifndef CMS2_HYPTRILEPMAKER_H
#define CMS2_HYPTRILEPMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"
typedef math::XYZTLorentzVector LorentzVector;

//
// class decleration
//

class HypTrilepMaker : public edm::EDProducer {
public:
  explicit HypTrilepMaker (const edm::ParameterSet&);
  ~HypTrilepMaker();

  bool goodElectronWithoutIsolation(const int els_tightId,
				    const int els_closestMuon,
				    const float els_d0);
  bool goodMuonWithoutIsolation(const float mus_gfit_chi2,
				const float mus_gfit_ndof,
				const float mus_d0,
				const int mus_validHits);
  bool passElectronIsolation(const float els_tkIso,
			     const LorentzVector &els_p4);
  bool passMuonIsolation(const float mus_iso03_sumPt,
			 const float mus_iso03_emEt,
			 const float mus_iso03_hadEt,
			 const LorentzVector &mus_p4);
  bool goodMuonIsolated(const float mus_gfit_chi2,
			const float mus_gfit_ndof,
			const float mus_d0,
			const int mus_validHits,
			const float mus_iso03_sumPt,
			const float mus_iso03_emEt,
			const float mus_iso03_hadEt,
			const LorentzVector &mus_p4);
  bool goodElectronIsolated(const int els_tightId,
			    const int els_closestMuon,
			    const float els_d0,
			    const float els_tkIso,
			    const LorentzVector &els_p4);
  unsigned int encodeTrileptonCandidate(unsigned int first,
					unsigned int second,
					unsigned int third);

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag jetsInputTag;
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;
  edm::InputTag elToMuAssInputTag;
  edm::InputTag metInputTag;

};


#endif
