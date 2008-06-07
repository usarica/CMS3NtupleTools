// -*- C++ -*-
//
// Package:    MuonMaker
// Class:      MuonMaker
// 
/**\class MuonMaker MuonMaker.cc CMS2/MuonMaker/src/MuonMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonMaker.cc,v 1.1 2008/06/07 10:04:02 jmuelmen Exp $
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
#include "CMS2/NtupleMaker/interface/MuonMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

typedef math::XYZTLorentzVector LorentzVector;
using std::vector;
//
// class decleration
//

//
// constructors and destructor
//
MuonMaker::MuonMaker(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
     produces<vector<LorentzVector> >	("musp4"	).setBranchAlias("mus_p4"	);
     produces<vector<float> >		("muschi2"	).setBranchAlias("mus_chi2"	);

   //now do what ever other initialization is needed

}


MuonMaker::~MuonMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     std::auto_ptr<vector<LorentzVector> > 	vector_muon_p4		(new vector<LorentzVector>);
     std::auto_ptr<vector<float> > 	 	vector_muon_chi2	(new vector<float>);
     Handle<edm::View<reco::Muon> > muon_h;
     iEvent.getByLabel("allLayer1TopMuons", muon_h);
     edm::View<reco::Muon>::const_iterator muons_end = muon_h->end();
     for (edm::View<reco::Muon>::const_iterator i = muon_h->begin(); 
	  i != muons_end; ++i) {
	  printf("%f\n", i->p4().pt());
	  vector_muon_p4	->push_back(i->p4());
	  vector_muon_chi2	->push_back(i->combinedMuon().isNonnull() ? i->combinedMuon()->chi2() : -999);
     }
     iEvent.put(vector_muon_p4		, "musp4"	);
     iEvent.put(vector_muon_chi2	, "muschi2"	);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMaker);
