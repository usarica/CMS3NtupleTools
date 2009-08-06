#ifndef CMS2_HYPISOPRODUCER_H
#define CMS2_HYPISOPRODUCER_H

// -*- C++ -*-
//

#include <vector>
#include <functional>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;

class HypIsoMaker : public edm::EDProducer {
 public:
  explicit HypIsoMaker(const edm::ParameterSet&);
  ~HypIsoMaker();


  virtual void produce(edm::Event&, const edm::EventSetup&);

  //double getHypSum( LorentzVector isop4, vector<LorentzVector> exclp4, bool returnEt, edm::ESHandle<CaloGeometry> theCaloGeom, CaloRecHitMetaCollectionV* caloHits ) const;
  
  //double getHypSum( int objid, int objidx, vector<int> excid, vector<int> excidx, bool returnEt, edm::Handle< edm::View<reco::Candidate> > emObjectHandle, edm::Handle< edm::View<reco::Muon> > muonHandle, edm::ESHandle<CaloGeometry> theCaloGeom, const CaloSubdetectorGeometry* subdet[], CaloRecHitMetaCollectionV* caloHits, edm::Event& iEvent, const edm::EventSetup& iSetup) const;

  double getHypSum( int objid, int objidx, vector<int> excid, vector<int> excidx, bool returnEt, edm::Handle< edm::View<reco::Candidate> > emObjectHandle, edm::Handle< edm::View<reco::Muon> > muonHandle, edm::Event& iEvent, const edm::EventSetup& iSetup) const;
  
  double getSum(const reco::Candidate * emobject, const std::vector<const reco::Candidate*> exclObjects,
				 bool returnEt, edm::ESHandle<CaloGeometry> theCaloGeom, CaloRecHitMetaCollectionV* caloHits ) const;

  //double track_iso(LorentzVector *prip4, float *pri_d0, float *pri_z0, int *priid, vector<LorentzVector> *excp4, int *excid,
  //				   vector<LorentzVector> *trks_trk_p4, vector<float> *trks_d0, vector<float> *trks_z0) ;
  //double track_iso(LorentzVector prip4, float pri_d0, float pri_z0, int priid, vector<LorentzVector> excp4, int excid,
  //vector<LorentzVector> trks_trk_p4, vector<float> trks_d0, vector<float> trks_z0) ;
  double track_iso(LorentzVector prip4, float pri_d0, float pri_z0, int priid, vector<LorentzVector> excp4, vector<int> excid,
				   const vector<LorentzVector> *trksp4, const vector<float> *trks_d0, const vector<float> *trks_z0) ;


 private:
  // ----------member data ---------------------------

  edm::InputTag electronsInputTag_;
  edm::InputTag muonsInputTag_;
  edm::InputTag cms2elsInputTag_;
  edm::InputTag cms2musInputTag_;
  edm::InputTag trackInputTag_;
  edm::InputTag hypInputTag_;
  edm::InputTag emObjectProducer_;
  edm::InputTag ecalBarrelRecHitProducer_;
  edm::InputTag ecalBarrelRecHitCollection_;  
  edm::InputTag ecalEndcapRecHitProducer_;
  edm::InputTag ecalEndcapRecHitCollection_;
  edm::InputTag caloTowersInputTag_;
  
  bool recomputeEcalIso_;
  bool recomputeTrckIso_;
  
  double elsEcalVetoRadBarrel_; //els inner cone size--Barrel
  double elsEcalVetoRadEndcap_; //els inner cone size--Endcap
  double elsEcalExtCone_      ; //els outer cone size
  double musEcalVetoRadBarrel_; //mus inner cone size--Barrel
  double musEcalVetoRadEndcap_; //mus inner cone size--Endcap
  double musEcalExtCone_      ; //mus outer cone size    
  double IsoJurassicWidth_; 	// exclusion strip width for jurassic veto
  
  double elsEtMinBarrel_; //minimum Et noise cut
  double elsEMinBarrel_;  //minimum E noise cut
  double elsEtMinEndcap_; //minimum Et noise cut
  double elsEMinEndcap_;  //minimum E noise cut

  double musEtMinBarrel_; //minimum Et noise cut
  double musEMinBarrel_;  //minimum E noise cut
  double musEtMinEndcap_; //minimum Et noise cut
  double musEMinEndcap_;  //minimum E noise cut

  double trackIsoExtRadius_;   
  double trackIsoElsInRadius_;
  double trackIsoMusInRadius_;
  double trackIsoMinPt_;
  double trackIsoMind0_;
  double trackIsoMinz0_;

  bool useIsolEt_; //switch for isolEt rather than isolE

  TrackAssociatorParameters muonparameters_;

};

#endif
