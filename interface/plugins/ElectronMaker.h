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
#include "FWCore/Framework/interface/stream/EDProducer.h"
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class decleration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class ElectronMaker : public edm::stream::EDProducer<> {
public:
  explicit ElectronMaker(const edm::ParameterSet&);
  ~ElectronMaker();

private:
  //  virtual void beginJob() ;
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  double electronIsoValuePF(
    reco::GsfElectron const&, reco::Vertex const&,
    float, float, float, float, float, float,
    int
  );

  int classify(edm::RefToBase<pat::Electron> const&);
  template<typename T> const edm::ValueMap<T>& getValueMap(edm::Event const&, edm::InputTag const&);

  void elIsoCustomCone(edm::View<pat::Electron>::const_iterator const&, float, bool, float, float&, float&, float&, float&);
  void elMiniIso(edm::View<pat::Electron>::const_iterator const&, bool, float, float&, float&, float&, float&);

  void setMVAIdUserVariables(edm::View<pat::Electron>::const_iterator const&, pat::Electron&, std::string const&, std::string const&) const;
  void setCutBasedIdUserVariables(edm::View<pat::Electron>::const_iterator const&, pat::Electron&, std::string const&, std::string const&) const;


  // ----------member data ---------------------------
  std::string aliasprefix_;
  int year_;

  //edm::InputTag beamSpot_tag_;

  edm::InputTag trksInputTag_;
  edm::InputTag gsftracksInputTag_;

  edm::InputTag eidLHTag_;
  //edm::InputTag cms2scsseeddetidInputTag_;

  edm::InputTag ebReducedRecHitCollectionTag;
  edm::InputTag eeReducedRecHitCollectionTag;
  edm::InputTag esReducedRecHitCollectionTag;

  edm::InputTag rhoInputTag_;


  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<edm::View<pat::Electron>  > electronsToken;
  //edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsToken;
  edm::EDGetTokenT<float> bFieldToken;
  edm::EDGetTokenT<LorentzVector> beamSpotToken;
  edm::EDGetTokenT<reco::ConversionCollection> recoConversionToken;

  edm::EDGetTokenT<edm::ValueMap<float> > miniIsoChgValueMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > miniIsoAllValueMapToken_;

  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetsToken;
  /*
  edm::InputTag pfIsoCharged03InputTag;
  edm::InputTag pfIsoGamma03InputTag;
  edm::InputTag pfIsoNeutral03InputTag;
  edm::InputTag pfIsoCharged04InputTag;
  edm::InputTag pfIsoGamma04InputTag;
  edm::InputTag pfIsoNeutral04InputTag;
  */

  double minAbsDist_;
  double minAbsDcot_;
  double minSharedFractionOfHits_;

  edm::Handle<reco::PFCandidateCollection> pfCand_h;
  edm::Handle<pat::PackedCandidateCollection> packPfCand_h;
  const pat::PackedCandidateCollection* pfCandidates;
  edm::Handle<reco::VertexCollection> vertexHandle;
  edm::Handle<reco::ConversionCollection> convs_h;

  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection;

};

#endif
