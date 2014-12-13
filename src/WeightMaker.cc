//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      WeightMaker
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/WeightsInfo.h"

#include "CMS3/NtupleMaker/interface/WeightMaker.h"

#include <vector>
#include "TString.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

WeightMaker::WeightMaker(const edm::ParameterSet& iConfig) {

  LHEEventInputTag_ = iConfig.getParameter<std::string>("LHEEventInputTag");
  
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<std::vector<TString>>     (branchprefix+"id" ).setBranchAlias(aliasprefix_+"_id" );
  produces<std::vector<float>>       (branchprefix+"wgt").setBranchAlias(aliasprefix_+"_wgt");
  produces<float>                    (branchprefix+"orig").setBranchAlias(aliasprefix_+"_orig");

}

WeightMaker::~WeightMaker() {}

void  WeightMaker::beginJob() {}

void WeightMaker::endJob() {}

void WeightMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  using namespace edm;
  
  auto_ptr<std::vector<TString>> weight_id   ( new std::vector<TString> );
  auto_ptr<std::vector<float>>   weight_wgt  ( new std::vector<float>   );
  auto_ptr<float>                weight_orig ( new float                );

  // get the LHEEventProduct
  edm::Handle<LHEEventProduct> hLHEEvent;
  iEvent.getByLabel(LHEEventInputTag_, hLHEEvent);

  // try to get using the LHEEvent
  if(hLHEEvent.isValid()) {
    
    // get the original weight
    *weight_orig = ((float) hLHEEvent->originalXWGTUP());

    // get the reweighting factors
    const std::vector<gen::WeightsInfo> weights                   = hLHEEvent->weights();
    for (std::vector<gen::WeightsInfo>::const_iterator weights_it = weights.begin(); weights_it != weights.end(); weights_it++) {      
      weight_id ->push_back(weights_it->id);
      weight_wgt->push_back((float) weights_it->wgt);
    }    
  }

  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"), 1, "");

  iEvent.put(weight_id    ,branchprefix+"id"    );
  iEvent.put(weight_wgt   ,branchprefix+"wgt"   );
  iEvent.put(weight_orig  ,branchprefix+"orig"   );

}

// register this maker
DEFINE_FWK_MODULE(WeightMaker);
