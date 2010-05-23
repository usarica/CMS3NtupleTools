// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      METMaker
// 
/**\class METMaker.cc CMS2/NtupleMaker/src/METMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: METMaker.h,v 1.11 2010/05/23 17:59:24 fgolf Exp $
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
     edm::InputTag met_tag;       
     edm::InputTag met36x_tag;
     edm::InputTag metHO_tag;     
     edm::InputTag metNoHF_tag;   
     edm::InputTag metNoHFHO_tag; 
     edm::InputTag metOpt_tag;       
     edm::InputTag metOptHO_tag;     
     edm::InputTag metOptNoHF_tag;   
     edm::InputTag metOptNoHFHO_tag; 
     edm::InputTag corMetGlobalMuons_tag;
     edm::InputTag MuonJEScorMET_tag;
     edm::InputTag muon_vm_tag;
     edm::InputTag muon_tag;
     edm::InputTag caloTowerInputTag;
     edm::InputTag hbheNoiseFilterInputTag;
     double towerEtThreshold;
     bool make_eta_rings;
     std::string aliasprefix_;
};

#endif
