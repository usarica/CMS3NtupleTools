// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      LumiFilter
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/LumiFilter.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: LumiFilter.h,v 1.16 2010/04/15 11:28:59 jribnik Exp $
//
//
#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
// #include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Scalers/interface/ScalersRaw.h"

#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/EventRange.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h" 
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"

#include "FWCore/Utilities/interface/EDMException.h"


#include "TString.h"

class LumiFilter : public edm::stream::EDFilter<> {
public:
    explicit LumiFilter (const edm::ParameterSet&);
    ~LumiFilter();

private:
    virtual void beginJob() ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    virtual void beginRun (const edm::Run& iRun, const edm::EventSetup& iSetup);

    std::vector<edm::LuminosityBlockRange> lumisToProcess_;
};


#endif
