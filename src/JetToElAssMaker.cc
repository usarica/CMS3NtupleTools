// -*- C++ -*-
//
// Package:    JetToElAssMaker
// Class:      JetToElAssMaker
// 
/**\class JetToElAssMaker JetToElAssMaker.cc CMS2/NtupleMaker/src/JetToElAssMaker.cc

 Description: make associations between jets and electrons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JetToElAssMaker.cc,v 1.5 2009/09/10 10:51:43 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/JetToElAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;

JetToElAssMaker::JetToElAssMaker(const edm::ParameterSet& iConfig)
{
     produces<std::vector<int>    >("jetsclosestElectron"  ).setBranchAlias("jets_closestElectron"   );
     produces<std::vector<double> >("jetsclosestElectronDR").setBranchAlias("jets_closestElectron_DR");

     m_minDR     = iConfig.getParameter<double>("minDR");
     jetInputTag = iConfig.getParameter<edm::InputTag>("jetInputTag_");
     elInputTag  = iConfig.getParameter<edm::InputTag>("elInputTag_" );
}

void JetToElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     std::auto_ptr<std::vector<int>    > vector_jets_closestElectron   (new std::vector<int>   );        
     std::auto_ptr<std::vector<double> > vector_jets_closestElectron_DR(new std::vector<double>);        

     edm::Handle<std::vector<LorentzVector> > els_p4_h;
     edm::Handle<std::vector<LorentzVector> > jets_p4_h;

     iEvent.getByLabel(elInputTag,  els_p4_h );  
     iEvent.getByLabel(jetInputTag, jets_p4_h);

     for(std::vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin(); jets_it != jets_p4_h->end(); jets_it++) {
       
       double jet_eta = jets_it->Eta();
       double jet_phi = jets_it->Phi();
       
       double minDR   = m_minDR;
       unsigned int i = 0;
       int index      = -1; 

       for(std::vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(); els_it != els_p4_h->end(); els_it++, i++) {
	 
	 double el_eta = els_it->Eta();
	 double el_phi = els_it->Phi();
	 double dR     = deltaR(jet_eta, jet_phi, el_eta, el_phi);

	 if(dR < minDR) {
	   minDR = dR;
	   index = i;
	 }
       }

       vector_jets_closestElectron->push_back(index);
       vector_jets_closestElectron_DR->push_back(minDR);
     }


     iEvent.put(vector_jets_closestElectron   , "jetsclosestElectron"  );
     iEvent.put(vector_jets_closestElectron_DR, "jetsclosestElectronDR");
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetToElAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetToElAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetToElAssMaker);
