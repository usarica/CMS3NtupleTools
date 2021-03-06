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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"

#include "Math/VectorUtil.h"


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
  void setEGammaPFElectronIdSelectionBits(edm::View<pat::Electron>::const_iterator const&, pat::PackedCandidate const*, pat::Electron&) const;
  void applyTriggerEmulationCuts(double const&, edm::View<pat::Electron>::const_iterator const&, pat::Electron&, unsigned int const&) const;

  static float getRecHitEnergyTime(DetId const&, EcalRecHitCollection const*, EcalRecHitCollection const*, unsigned short, unsigned short, float* outtime=nullptr);

protected:
  std::string aliasprefix_;
  int year_;

  edm::VParameterSet MVACuts_;

  std::unordered_map< std::string, std::vector< StringCutObjectSelector<pat::Electron, true> > > MVACutObjects;

  edm::EDGetTokenT<edm::View<pat::Electron>  > electronsToken;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<reco::ConversionCollection> recoConversionToken;

  edm::EDGetTokenT<EcalRecHitCollection> ebhitsToken;
  edm::EDGetTokenT<EcalRecHitCollection> eehitsToken;
  //edm::EDGetTokenT<EcalRecHitCollection> eshitsToken;

  edm::EDGetTokenT< double > rhoToken;
  edm::EDGetTokenT< double > rhoCaloToken;

};


#endif
