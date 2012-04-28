// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElectronIsolationMaker
// 
/**\class ElectronIsolationMaker ElectronIsolationMaker.cc CMS2/NtupleMaker/src/ElectronIsolationMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronIsolationMaker.h,v 1.1 2012/04/28 07:55:53 fgolf Exp $
//
//
#ifndef CMS2_ELECTRONISOLATIONMAKER_H
#define CMS2_ELECTRONISOLATIONMAKER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class ElectronIsolationMaker : public edm::EDProducer {
public:
    explicit ElectronIsolationMaker (const edm::ParameterSet&);

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
  
    // ----------member data ---------------------------
    edm::InputTag gsfElectronInputTag;
    edm::InputTag cms2electronInputTag;
    edm::InputTag pfNoPileUpInputTag;

    std::string aliasprefix_;
    std::string branchprefix_;
};


#endif
