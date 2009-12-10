// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GenMaker
// 
/**\class GenMaker GenMaker.cc CMS2/NtupleMaker/src/GenMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GenMaker.h,v 1.11 2009/12/10 02:15:32 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_GENMAKER_H
#define NTUPLEMAKER_GENMAKER_H

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

class GenMaker : public edm::EDProducer {
public:
     explicit GenMaker (const edm::ParameterSet&);
     ~GenMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     virtual void beginRun(edm::Run&, const edm::EventSetup&);

     // ----------member data ---------------------------
     edm::InputTag genParticlesInputTag;
     edm::InputTag genEventScaleInputTag;

     bool ntupleOnlyStatus3;
     bool ntupleDaughters;
 
     std::vector<int> vmetPIDs;

     double inclusiveCrossSectionValue;
     double exclusiveCrossSectionValue;
     double kfactorValue;

};

#endif

