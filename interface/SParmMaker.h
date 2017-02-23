// -*- C++ -*-
//
// Package:    SParmMaker
// Class:      SParmMaker
// 
/**\class SParmMaker SParmMaker.h CMS3/NtupleMaker/interface/SParmMaker.h

   Description: copy SUSY mSUGRA parameters into the EDM event tree

   Implementation:
   - extract and fill variables
 
*/
//
// Original BenHooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
//
//
#ifndef CMS2_SPARMMAKER_H
#define CMS2_SPARMMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/BranchType.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"

//
// class declaration
//

class SParmMaker : public edm::EDProducer {
public:
    explicit SParmMaker (const edm::ParameterSet&);
    ~SParmMaker();

private:
    virtual void beginJob() ;
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, const edm::EventSetup&) ;
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, const edm::EventSetup&) ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
      
    // ----------member data ---------------------------
    edm::EDGetTokenT<LHEEventProduct> sparmToken;
    edm::EDGetTokenT<GenLumiInfoHeader> configToken;
    std::string aliasprefix_;
    std::vector<std::string> vsparms_;
    float filtEff;
};


#endif
