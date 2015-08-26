// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      METMaker
// 
/**\class METMaker.cc CMS3/NtupleMaker/src/METMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: METMaker.h,v 1.12 2010/05/24 14:21:23 fgolf Exp $
//
//

#ifndef NTUPLEMAKER_METMAKER_H
#define NTUPLEMAKER_METMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/CaloMET.h"

//
// class decleration
//

class METMaker : public edm::EDProducer {
public:
     explicit METMaker (const edm::ParameterSet&);
     ~METMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     // ----------member data ---------------------------  
     edm::EDGetTokenT<edm::View<reco::CaloMET> > met_token;       
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metHO_token;     
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metNoHF_token;   
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metNoHFHO_token; 
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metOpt_token;       
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metOptHO_token;     
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metOptNoHF_token;   
     edm::EDGetTokenT<edm::View<reco::CaloMET> > metOptNoHFHO_token; 
     edm::EDGetTokenT<edm::View<reco::CaloMET> > corMetGlobalMuons_token;
     edm::EDGetTokenT<edm::View<reco::CaloMET> > MuonJEScorMET_token;
     edm::EDGetTokenT<edm::View<reco::CaloMET> > muon_vm_token;
     edm::EDGetTokenT<edm::View<reco::CaloMET> > muon_token;
     edm::EDGetTokenT<bool> hbheNoiseFilterToken;
     edm::InputTag caloTowerInputTag;

     double towerEtThreshold;
     bool make_eta_rings;
     std::string aliasprefix_;
};

#endif
