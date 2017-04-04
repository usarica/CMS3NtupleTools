//-*- C++ -*-
// $Id: EnergyDensityMaker.cc,v 1.3 2012/05/17 00:25:14 macneill Exp $
//

#include <memory>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class EnergyDensityMaker : public edm::stream::EDProducer<> {
public:
  explicit EnergyDensityMaker (const edm::ParameterSet&);
  ~EnergyDensityMaker(){}

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  edm::EDGetTokenT<double> inputToken;
  std::string m_alias;
  std::string m_branch;
};

EnergyDensityMaker::EnergyDensityMaker(const edm::ParameterSet& iConfig) {
  m_alias = iConfig.getUntrackedParameter<std::string>("alias");
  m_branch = m_alias;
  size_t pos = 0;
  while( (pos = m_branch.find("_")) != std::string::npos)
    m_branch.replace(pos,1,"");
  produces<float>(m_branch).setBranchAlias(m_alias);
  inputToken = consumes<double>(iConfig.getParameter<edm::InputTag>("input"));
}

void EnergyDensityMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::unique_ptr<float> rho(new float(-9999.));
  edm::Handle<double> rhoH;
  iEvent.getByToken( inputToken , rhoH);
  if(rhoH.isValid()){
	*rho = *rhoH; 
  }
  else {
    throw cms::Exception("EnergyDensityMaker: error getting rho collection!");
    return;
  }
  iEvent.put(std::move(rho), m_branch);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EnergyDensityMaker);





  
