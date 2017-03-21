// -*- C++ -*-
//
// Package:    NtupleMaker
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
// $Id: MuonMaker.h,v 1.19 2012/07/20 01:24:31 dbarge Exp $
//
//
#ifndef CMS2_MUONMAKER_H
#define CMS2_MUONMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


//
// class declaration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class MuonMaker : public edm::EDProducer {
public:
     explicit MuonMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  double muonIsoValuePF(const reco::Muon& mu, const reco::Vertex& vtx, float coner, float minptn, float dzcut, int filterId);

  void muIsoCustomCone( edm::View<pat::Muon>::const_iterator& mu, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso);
  void muMiniIso( edm::View<pat::Muon>::const_iterator& mu, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float & dbiso);
  
      // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon> > muonsToken;
  edm::InputTag beamSpotInputTag;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  std::string tevMuonsName;

  std::string aliasprefix_;
  std::string branchprefix_;

  edm::Handle<reco::PFCandidateCollection> pfCand_h;
  edm::Handle<pat::PackedCandidateCollection> packPfCand_h;
  edm::Handle<reco::VertexCollection> vertexHandle;
  const pat::PackedCandidateCollection *pfCandidates;

  // Cosmics Compatibility
  edm::InputTag src_;
  //  
  edm::InputTag showerTag_;

};


#endif
