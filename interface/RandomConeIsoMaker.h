// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      RandomConeIsoMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/RandomConeIsoMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: RandomConeIsoMaker.h,v 1.3 2010/03/02 19:24:12 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TString.h"
//
// class decleration
//

class RandomConeIsoMaker : public edm::EDProducer {
public:
  explicit RandomConeIsoMaker (const edm::ParameterSet&);
  ~RandomConeIsoMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  reco::isodeposit::IsoDepositExtractor* muIsoExtractorCalo_;
  reco::isodeposit::IsoDepositExtractor* muIsoExtractorTrack_;
  reco::isodeposit::IsoDepositExtractor* muIsoExtractorJet_;
  edm::InputTag ecalBarrelRecHitProducer_;
  edm::InputTag ecalBarrelRecHitCollection_;  
  edm::InputTag ecalEndcapRecHitProducer_;
  edm::InputTag ecalEndcapRecHitCollection_;
  edm::InputTag hbheRecHit_ ;
  edm::InputTag towerProducer_;
  edm::InputTag trackProducer_;
  edm::InputTag beamspotProducer_;
  edm::InputTag srProducerEE_;
  edm::InputTag srProducerEB_;

  double egIsoPtMinBarrel_; //minimum Et noise cut
  double egIsoEMinBarrel_;  //minimum E noise cut
  double egIsoPtMinEndcap_; //minimum Et noise cut
  double egIsoEMinEndcap_;  //minimum E noise cut
  double egIsoConeSizeOut_; //outer cone size
  double egIsoConeSizeInBarrel_; //inner cone size
  double egIsoConeSizeInEndcap_; //inner cone size
  double egIsoJurassicWidth_ ; // exclusion strip width for jurassic veto
  bool useIsolEt_; //switch for isolEt rather than isolE
  bool tryBoth_ ; // use rechits from barrel + endcap 
  bool subtract_ ; // subtract SC energy (allows veto cone of zero size)
  bool useNumCrystals_ ; // veto on number of crystals
  bool vetoClustered_ ;  // veto all clusterd rechits
  double egHcalIsoPtMin_;
  double egHcalIsoConeSizeOut_;
  double egHcalIsoConeSizeIn_;
  signed int egHcalDepth_; 

  double ptMin_;
  double intRadius_;
  double extRadius_;
  double maxVtxDist_;
  double drb_;

  // for track association / propagation to calo
  TrackAssociatorParameters parameters_;
  TrackDetectorAssociator trackAssociator_;


  CLHEP::HepJamesRandom 	*jamesRandom_;
  edm::InputTag primaryVertexInputTag_;	

  double signedRnd(double number, const double &min, const double &max);
  std::pair<int, double> getTrackIso(const reco::TrackBase::Point &beamspot,
				     const reco::TrackCollection *trackCollection, 
				     const math::XYZVector &tmpElectronMomentumAtVtx) const;
  
  // std::string datasetName_;
//   std::string CMS2tag_;
};


#endif
