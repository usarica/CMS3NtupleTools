// -*- C++ -*-
// $Id: ElCaloIsoMaker.cc,v 1.5 2010/03/18 02:11:56 kalavase Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "CMS2/NtupleMaker/interface/ElCaloIsoMaker.h"

ElCaloIsoMaker::ElCaloIsoMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

   produces<std::vector<float> >  (branchprefix+"ecalJuraIso").setBranchAlias(aliasprefix_+"_ecalJuraIso");
   produces<std::vector<float> >  (branchprefix+"ecalJuraTowerIso").setBranchAlias(aliasprefix_+"_ecalJuraTowerIso");
   produces<std::vector<float> >  (branchprefix+"hcalConeIso").setBranchAlias(aliasprefix_+"_hcalConeIso");
   m_electronsInputTag =    iConfig.getParameter<edm::InputTag>("electronsInputTag");
   m_basicClusterInputTag = iConfig.getParameter<edm::InputTag>("basicClusterInputTag");
   m_caloTowersInputTag   = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
   m_maxDR =                iConfig.getParameter<double>("maxDR");
   m_minDR =                iConfig.getParameter<double>("minDR");
   m_minDEta =              iConfig.getParameter<double>("minDEta");
}


ElCaloIsoMaker::~ElCaloIsoMaker() {}

void ElCaloIsoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   produceEcalIso(iEvent,iSetup);
   produceEcalTowerIso(iEvent,iSetup);
   produceHcalIso(iEvent,iSetup);
}


void ElCaloIsoMaker::produceEcalIso(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   std::auto_ptr<std::vector<float> >  els_juraIso( new std::vector<float> ) ;
   
   edm::Handle<edm::View<reco::GsfElectron> > electron_h;
   iEvent.getByLabel(m_electronsInputTag, electron_h);
   
   edm::Handle<reco::BasicClusterCollection> basicClusterHandle;
   iEvent.getByLabel(m_basicClusterInputTag, basicClusterHandle);
   
 std::string branchprefix = aliasprefix_;
   if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
   

   for(edm::View<reco::GsfElectron>::const_iterator electron = electron_h->begin(); 
       electron != electron_h->end(); ++electron){
      math::XYZPoint positionAtEcal = electron->caloPosition();
      double juraIso(0);
      // loop over basic clusters
      for(reco::BasicClusterCollection::const_iterator cluster = basicClusterHandle->begin(); 
	  cluster != basicClusterHandle->end(); ++cluster)
	{
	   double dR = deltaR( positionAtEcal.eta(), positionAtEcal.phi(), cluster->eta(), cluster->phi() );
	   if ( dR > m_maxDR ) continue;
	   if ( dR < m_minDR ) continue;
	   if ( fabs( positionAtEcal.eta() - cluster->eta() ) < m_minDEta ) continue;
	   
	   juraIso += cluster->energy()*sin(cluster->position().theta());
	}
      els_juraIso->push_back( juraIso );
   }
   
   iEvent.put(els_juraIso, branchprefix+"ecalJuraIso");
}

void ElCaloIsoMaker::produceHcalIso(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   std::auto_ptr<std::vector<float> >  els_hcalIso( new std::vector<float> ) ;
   
   edm::Handle<edm::View<reco::GsfElectron> > electron_h;
   iEvent.getByLabel(m_electronsInputTag, electron_h);
   
   edm::Handle<CaloTowerCollection> caloTowers;
   iEvent.getByLabel(m_caloTowersInputTag, caloTowers);
   
   std::string branchprefix = aliasprefix_;
   if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
   
   for(edm::View<reco::GsfElectron>::const_iterator electron = electron_h->begin(); 
       electron != electron_h->end(); ++electron){
      math::XYZPoint positionAtEcal = electron->caloPosition();
      double hcalIso(0);
      // loop over towers
      for(CaloTowerCollection::const_iterator tower = caloTowers->begin();
	  tower != caloTowers->end(); ++tower)
	{
	   double dR = deltaR( positionAtEcal.eta(), positionAtEcal.phi(), tower->eta(), tower->phi() );
	   if ( dR > m_maxDR ) continue;
	   
	   hcalIso += tower->hadEnergy() + tower->outerEnergy();
	}
      els_hcalIso->push_back( hcalIso );
   }
   iEvent.put(els_hcalIso, branchprefix+"hcalConeIso");
}

void ElCaloIsoMaker::produceEcalTowerIso(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   std::auto_ptr<std::vector<float> >  els_juraIso( new std::vector<float> ) ;
   
   edm::Handle<edm::View<reco::GsfElectron> > electron_h;
   iEvent.getByLabel(m_electronsInputTag, electron_h);
   
   edm::Handle<CaloTowerCollection> caloTowers;
   iEvent.getByLabel(m_caloTowersInputTag, caloTowers);
   
   std::string branchprefix = aliasprefix_;
   if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

   for(edm::View<reco::GsfElectron>::const_iterator electron = electron_h->begin(); 
       electron != electron_h->end(); ++electron){
      math::XYZPoint positionAtEcal = electron->caloPosition();
      double juraIso(0);
      // loop over towers
      for(CaloTowerCollection::const_iterator tower = caloTowers->begin();
	  tower != caloTowers->end(); ++tower)
	{
	   double dR = deltaR( positionAtEcal.eta(), positionAtEcal.phi(), tower->eta(), tower->phi() );
	   if ( dR > m_maxDR ) continue;
	   if ( dR < m_minDR ) continue;
	   if ( fabs( positionAtEcal.eta() - tower->eta() ) < m_minDEta ) continue;
	   
	   juraIso += tower->hadEnergy() + tower->outerEnergy();
	}
      els_juraIso->push_back( juraIso );
   }
   iEvent.put(els_juraIso, branchprefix+"ecalJuraTowerIso");
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElCaloIsoMaker);
