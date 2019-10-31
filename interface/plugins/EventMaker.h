// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/EventMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventMaker.h,v 1.16 2010/04/15 11:28:59 jribnik Exp $
//
//
#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
// #include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Scalers/interface/ScalersRaw.h"


#include "TString.h"
//
// class decleration
//

class EventMaker : public edm::stream::EDProducer<> {
public:
    explicit EventMaker (const edm::ParameterSet&);
    ~EventMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    virtual void beginRun (const edm::Run& iRun, const edm::EventSetup& iSetup);
    // virtual void respondToOpenInputFile(const edm::FileBlock& iFileBlock);

    std::string datasetName_;
    std::string CMS3tag_;

    edm::EDGetTokenT<DcsStatusCollection> dcsTag_;
    edm::EDGetTokenT<LumiScalersCollection> scalersTag_;
    std::string aliasprefix_;
    bool isData_;
};


#endif
