// -*- C++ -*-
//
// Package:    ElToPFElAssMaker
// Class:      ElToPFElAssMaker
// 
/**\class ElToPFElAssMaker ElToPFElAssMaker.cc CMS2/NtupleMaker/src/ElToPFElAssMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToPFElAssMaker.cc,v 1.1 2010/06/13 16:00:49 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/ElToPFElAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

ElToPFElAssMaker::ElToPFElAssMaker(const edm::ParameterSet& iConfig) {

     aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
     std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     produces<vector<int>   >(branchprefix+"pfelsidx").setBranchAlias(aliasprefix_+"_pfelsidx");
     
     elsInputTag_ = iConfig.getParameter<edm::InputTag>("elsInputTag");
     pfelsInputTag_ = iConfig.getParameter<edm::InputTag>("pfelsInputTag");
     m_minDR_     = iConfig.getParameter<double>       ("minDR");
}

void ElToPFElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>    > vector_els_pfelsidx(new vector<int>   );
     
     // get reco electron p4's
     edm::Handle<vector<LorentzVector> > els_p4_h;
     iEvent.getByLabel(elsInputTag_.label(), "elsp4", els_p4_h);

     // get particle flow electron p4's
     edm::Handle<vector<LorentzVector> > pfels_p4_h;
     iEvent.getByLabel(pfelsInputTag_.label(), "pfelsp4", pfels_p4_h);
     
     //loop over reco electrons and find the closest particle flow electron
     for (vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(); els_it != els_p4_h->end(); els_it++)
     {       
	  double el_eta = els_it->Eta();
	  double el_phi = els_it->Phi();
       
	  double minDR = 9999.;
	  unsigned int i = 0;
	  int index = -1; 

	  for (vector<LorentzVector>::const_iterator pfels_it = pfels_p4_h->begin(); pfels_it != pfels_p4_h->end(); pfels_it++, i++) {
	 
	       double pfel_eta = pfels_it->Eta();
	       double pfel_phi = pfels_it->Phi();

	       double dR = deltaR(el_eta, el_phi, pfel_eta, pfel_phi);

	       if(dR < minDR) {
		    minDR = dR;
		    index = i;
	       }
	  }

	  if(minDR > m_minDR_) {
	       minDR = -9999.;
	       index = -1;
	  }

	  // fill vector
	  vector_els_pfelsidx->push_back(index);
     }

     // store vectors
     std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     iEvent.put(vector_els_pfelsidx, branchprefix+"pfelsidx");
}

// ------------ method called once each job just before starting event loop  ------------
void ElToPFElAssMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ElToPFElAssMaker::endJob() {}


//define this as a plug-in
DEFINE_FWK_MODULE(ElToPFElAssMaker);
