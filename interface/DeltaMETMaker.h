// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      DeltaMETMaker
// 
/**\class DeltaMETMaker.cc CMS2/NtupleMaker/src/DeltaMETMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: DeltaMETMaker.h,v 1.1 2010/05/19 21:57:56 fgolf Exp $
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

class DeltaMETMaker : public edm::EDProducer {
public:
     explicit DeltaMETMaker (const edm::ParameterSet&);
     ~DeltaMETMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     // ----------member data ---------------------------  
     edm::InputTag cms2_metInputTag;       
     edm::InputTag cms2_metphiInputTag;     
     edm::InputTag cms2_sumetInputTag;   
     edm::InputTag metInputTag; 
     std::string aliasprefix_;
};

#endif
