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
// $Id: METMaker.h,v 1.1 2008/06/19 20:04:29 kalavase Exp $
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
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  edm::InputTag genParticlesInputTag;
      // ----------member data ---------------------------
  
};


#endif
