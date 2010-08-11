//-*- C++ -*-
//
// Package:    TrkJetMaker
// Class:      TrkJetMaker.cc
//
/**\class TrkJetMaker TrkJetMaker.cc CMS2/NtupleMaker/src/TrkJetMaker.cc

   Description: Dumps the TrkJet contents into the ntuple
   Implementation:
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: TrkJetMaker.cc,v 1.9 2010/08/11 18:41:40 fgolf Exp $
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/TrkJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

typedef math::XYZTLorentzVectorF LorentzVector;

bool sortTrkJetsByPt(LorentzVector jet1, LorentzVector jet2) {
	 return jet1.pt() > jet2.pt();
}

//
// constructors and destructor
//
TrkJetMaker::TrkJetMaker(const edm::ParameterSet& iConfig)
{
	 // product of this EDProducer
	 produces<unsigned int>						("evtntrkjets"			).setBranchAlias("evt_ntrkjets"			);
	 produces<std::vector<LorentzVector> >		("trkjetsp4"			).setBranchAlias("trkjets_p4"			);
	 produces<std::vector<float> >				("trkjetscor"			).setBranchAlias("trkjets_cor"			);             
	 produces<std::vector<int> >				("trkjetsntrks"			).setBranchAlias("trkjets_ntrks"		);	 
	 produces<std::vector<int> >				("trkjetsvtxsidx"		).setBranchAlias("trkjets_vtxs_idx"		);

	 // parameters from configuration
	 trkJetsInputTag       = iConfig.getParameter<edm::InputTag>("trkJetsInputTag");
	 trkJetCorrectionL2L3_ = iConfig.getParameter<std::string>  ("trkJetCorrectionL2L3");

}


TrkJetMaker::~TrkJetMaker()
{

}

void
TrkJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	 // get TrkJet collection
	 edm::Handle<edm::View<reco::TrackJet> > trkJets;
	 iEvent.getByLabel(trkJetsInputTag, trkJets);
	 
	 edm::Handle<std::vector<LorentzVector> > vtxs_h;
	 iEvent.getByLabel("vtxsposition", vtxs_h);

	 // create containers
	 std::auto_ptr<unsigned int>                evt_ntrkjets			(new unsigned int(trkJets->size())		);
	 std::auto_ptr<std::vector<LorentzVector> > vector_trkjets_p4		(new std::vector<LorentzVector>			);
	 std::auto_ptr<std::vector<float> >         vector_trkjets_cor		(new std::vector<float>					);
	 std::auto_ptr<std::vector<int> >			vector_trkjets_ntrks	(new std::vector<int>					);
	 std::auto_ptr<std::vector<int> >			vector_trkjets_vtxs_idx	(new std::vector<int>					);	 

	 const JetCorrector* correctorL2L3 = JetCorrector::getJetCorrector (trkJetCorrectionL2L3_, iSetup);
  
	 // loop over jets and fill containers
	 edm::View<reco::TrackJet>::const_iterator jetsEnd = trkJets->end();
	 for ( edm::View<reco::TrackJet>::const_iterator jet = trkJets->begin(); jet != jetsEnd; ++jet) {

		  double cor = correctorL2L3	->correction(jet->p4()					);
		  vector_trkjets_cor			->push_back( cor                        );
		  vector_trkjets_p4				->push_back( LorentzVector( jet->p4() ) );
		  vector_trkjets_ntrks			->push_back( jet->numberOfTracks()		);

		  reco::VertexRef refvtx = jet->primaryVertex();
		  float trkjets_vtx_x = (*refvtx).x();
		  float trkjets_vtx_y = (*refvtx).y();

		  int vtx_idx = -1;
		  int nVtxs = vtxs_h->size();
		  for (int i_vtx = 0; i_vtx < nVtxs; i_vtx++)
		  {
			   float vtxs_x = vtxs_h->at(i_vtx).x();
			   float vtxs_y = vtxs_h->at(i_vtx).y();

			   if (fabs(trkjets_vtx_x - vtxs_x) < 0.001 && fabs(trkjets_vtx_y - vtxs_y) < 0.001)
					vtx_idx = i_vtx;
		  }

		  vector_trkjets_vtxs_idx->push_back(vtx_idx);
	 }

//	 std::sort( vector_trkjets_p4->begin(), vector_trkjets_p4->end(), sortTrkJetsByPt );

	 // put containers into event
	 iEvent.put(evt_ntrkjets			, "evtntrkjets"			);
	 iEvent.put(vector_trkjets_p4		,  "trkjetsp4"			);
	 iEvent.put(vector_trkjets_cor		, "trkjetscor"			);
	 iEvent.put(vector_trkjets_ntrks	, "trkjetsntrks"		);
	 iEvent.put(vector_trkjets_vtxs_idx	, "trkjetsvtxsidx"		);

}


void
TrkJetMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrkJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkJetMaker);
