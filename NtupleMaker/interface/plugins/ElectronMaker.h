#ifndef NTUPLEMAKER_ELECTRONMAKER_H
#define NTUPLEMAKER_ELECTRONMAKER_H

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <CommonTools/Utils/interface/StringCutObjectSelector.h>

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

  int classify(edm::RefToBase<pat::Electron> const&);

  void setupMVACuts();
  void setMVAIdUserVariables(edm::View<pat::Electron>::const_iterator const&, pat::Electron&, std::string const&, std::string const&) const;
  void setCutBasedIdUserVariables(edm::View<pat::Electron>::const_iterator const&, pat::Electron&, std::string const&, std::string const&) const;
  void applyTriggerEmulationCuts(double const&, edm::View<pat::Electron>::const_iterator const&, pat::Electron&) const;

protected:
  std::string aliasprefix_;
  int year_;

  edm::VParameterSet MVACuts_;

  edm::InputTag trksInputTag_;
  edm::InputTag gsftracksInputTag_;

  edm::InputTag ebReducedRecHitCollectionTag;
  edm::InputTag eeReducedRecHitCollectionTag;
  edm::InputTag esReducedRecHitCollectionTag;

  edm::InputTag rhoInputTag_;
  edm::InputTag rhoCaloInputTag_;


  std::unordered_map< std::string, std::vector< StringCutObjectSelector<pat::Electron, true> > > MVACutObjects;

  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<edm::View<pat::Electron>  > electronsToken;
  edm::EDGetTokenT<LorentzVector> beamSpotToken;
  edm::EDGetTokenT<reco::ConversionCollection> recoConversionToken;

  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection;

  edm::EDGetTokenT< double > rhoToken;
  edm::EDGetTokenT< double > rhoCaloToken;

};


#endif
