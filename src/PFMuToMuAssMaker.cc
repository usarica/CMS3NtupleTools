// -*- C++ -*-
//
// Package:    PFMuToMuAssMaker
// Class:      PFMuToMuAssMaker
// 
/**\class PFMuToMuAssMaker PFMuToMuAssMaker.cc CMS2/NtupleMaker/src/PFMuToMuAssMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMuToMuAssMaker.cc,v 1.1 2010/06/13 16:00:49 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/PFMuToMuAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

PFMuToMuAssMaker::PFMuToMuAssMaker(const edm::ParameterSet& iConfig) {

     aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
     std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     produces<vector<int>   >(branchprefix+"musidx").setBranchAlias(aliasprefix_+"_musidx");
     
     musInputTag_ = iConfig.getParameter<edm::InputTag>("musInputTag");
     pfmusInputTag_ = iConfig.getParameter<edm::InputTag>("pfmusInputTag");
     m_minDR_     = iConfig.getParameter<double>       ("minDR");
}

void PFMuToMuAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<int>    > vector_pfmus_musidx(new vector<int>   );
     
     // get reco electron p4's
     edm::Handle<vector<LorentzVector> > mus_p4_h;
     iEvent.getByLabel(musInputTag_.label(), "musp4", mus_p4_h);

     // get particle flow electron p4's
     edm::Handle<vector<LorentzVector> > pfmus_p4_h;
     iEvent.getByLabel(pfmusInputTag_.label(), "pfmusp4", pfmus_p4_h);
     
     //loop over reco electrons and find the closest particle flow electron
     for (vector<LorentzVector>::const_iterator pfmus_it = pfmus_p4_h->begin(); pfmus_it != pfmus_p4_h->end(); pfmus_it++)
     {       
	  double pfmu_eta = pfmus_it->Eta();
	  double pfmu_phi = pfmus_it->Phi();
       
	  double minDR = 9999.;
	  unsigned int i = 0;
	  int index = -1; 

	  for (vector<LorentzVector>::const_iterator mus_it = mus_p4_h->begin(); mus_it != mus_p4_h->end(); mus_it++, i++) {
	 
	       double mu_eta = mus_it->Eta();
	       double mu_phi = mus_it->Phi();

	       double dR = deltaR(pfmu_eta, pfmu_phi, mu_eta, mu_phi);

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
	  vector_pfmus_musidx->push_back(index);
     }

     // store vectors
     std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     iEvent.put(vector_pfmus_musidx, branchprefix+"musidx");
}

// ------------ method called once each job just before starting event loop  ------------
void PFMuToMuAssMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFMuToMuAssMaker::endJob() {}


//define this as a plug-in
DEFINE_FWK_MODULE(PFMuToMuAssMaker);
