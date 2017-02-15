// -*- C++ -*-    
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.h,v 1.27 2012/08/16 00:00:26 slava77 Exp $
//
//
#ifndef NTUPLEMAKER_ELECTRONMAKER_H
#define NTUPLEMAKER_ELECTRONMAKER_H

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class decleration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class ElectronMaker : public edm::EDProducer {
public:
    explicit ElectronMaker (const edm::ParameterSet&);
    ~ElectronMaker();

private:
//  virtual void beginJob() ;
    virtual void beginJob() ;
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    double electronIsoValuePF(const reco::GsfElectron& el, const reco::Vertex& vtx, float coner, float minptn, float dzcut,
                              float footprintdr, float gammastripveto, float elestripveto, int filterId);
  
    int classify(const edm::RefToBase<pat::Electron> &);
    template<typename T> const edm::ValueMap<T>& getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag);
 
    // for 2012 pf isolation
    void PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection, 
                         const int vertexIndex, const float &R, float &pfiso_ch, float &pfiso_em, float &pfiso_nh);
 
  void elIsoCustomCone(edm::View<pat::Electron>::const_iterator& el, float dr, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float &dbiso);
  void elMiniIso(edm::View<pat::Electron>::const_iterator& el, bool useVetoCones, float ptthresh, float &chiso, float &nhiso, float &emiso, float &dbiso);



    // ----------member data ---------------------------
    edm::InputTag beamSpotInputTag_;
    edm::InputTag trksInputTag_;
    edm::InputTag gsftracksInputTag_;
    edm::InputTag eidLHTag_;
    edm::InputTag cms2scsseeddetidInputTag_;

    edm::EDGetTokenT<reco::VertexCollection> vtxToken;
    edm::EDGetTokenT<edm::View<pat::Electron>  > electronsToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsToken;
    edm::EDGetTokenT<float> bFieldToken;
    edm::EDGetTokenT<LorentzVector> beamSpotToken;
    edm::EDGetTokenT<reco::ConversionCollection> recoConversionToken;


    edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronVIDNonTrigMvaWP80IdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronVIDNonTrigMvaWP90IdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronVIDTrigMvaWP80IdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronVIDTrigMvaWP90IdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > electronVIDNonTrigMvaValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > electronVIDTrigMvaValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int>  > electronVIDNonTrigMvaCatMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int>  > electronVIDTrigMvaCatMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > electronVIDSpring16GPMvaValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int>  > electronVIDSpring16GPMvaCatMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > electronVIDSpring16HZZMvaValueMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int>  > electronVIDSpring16HZZMvaCatMapToken_;

  edm::InputTag pfIsoCharged03InputTag;
  edm::InputTag pfIsoGamma03InputTag;
  edm::InputTag pfIsoNeutral03InputTag;
  edm::InputTag pfIsoCharged04InputTag;
  edm::InputTag pfIsoGamma04InputTag;
  edm::InputTag pfIsoNeutral04InputTag;

    EcalClusterLazyTools* clusterTools_;
    MultiTrajectoryStateTransform *mtsTransform_;

    double minAbsDist_;
    double minAbsDcot_;
    double minSharedFractionOfHits_;
    std::string aliasprefix_;
  bool useVID_;

    std::vector<Int_t> passVetoId_;
    std::vector<Int_t> passLooseId_;
    std::vector<Int_t> passMediumId_;
    std::vector<Int_t> passTightId_;
    std::vector<Int_t> passHEEPId_;
    std::vector<Int_t> passVIDNonTrigMvaWP80Id_;
    std::vector<Int_t> passVIDNonTrigMvaWP90Id_;
    std::vector<Int_t> passVIDTrigMvaWP80Id_;
    std::vector<Int_t> passVIDTrigMvaWP90Id_;
    std::vector<Float_t> VIDNonTrigMvaValue_;
    std::vector<Float_t> VIDTrigMvaValue_;
    std::vector<Float_t> VIDSpring16GPMvaValue_;
    std::vector<Float_t> VIDSpring16HZZMvaValue_;
    std::vector<Int_t> VIDNonTrigMvaCat_;
    std::vector<Int_t> VIDTrigMvaCat_;
    std::vector<Int_t> VIDSpring16GPMvaCat_;
    std::vector<Int_t> VIDSpring16HZZMvaCat_;

    edm::Handle<reco::PFCandidateCollection> pfCand_h;
    edm::Handle<pat::PackedCandidateCollection> packPfCand_h;
    const pat::PackedCandidateCollection *pfCandidates;
    edm::Handle<reco::VertexCollection> vertexHandle;
    edm::Handle<reco::ConversionCollection> convs_h;

    edm::InputTag rhoInputTag_;
    edm::InputTag beamSpot_tag_;

  edm::InputTag ebReducedRecHitCollectionTag;
  edm::InputTag eeReducedRecHitCollectionTag;
  edm::InputTag esReducedRecHitCollectionTag;
  
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection;

    PFPileUpAlgo *pfPileUpAlgo_;
};

#endif

