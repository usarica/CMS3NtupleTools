//-*- C++ -*-
// $Id: EnergyDensityMaker.cc,v 1.1 2011/06/13 12:29:59 dmytro Exp $
//

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class EnergyDensityMaker : public edm::EDProducer {
public:
  explicit EnergyDensityMaker (const edm::ParameterSet&);
  ~EnergyDensityMaker(){}

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  edm::InputTag m_input;
  std::string m_alias;
  std::string m_branch;
};

EnergyDensityMaker::EnergyDensityMaker(const edm::ParameterSet& iConfig) {
  m_alias = iConfig.getUntrackedParameter<std::string>("alias");
  m_branch = m_alias;
  if(m_branch.find("_") != std::string::npos)
    m_branch.replace(m_branch.find("_"),1,"");
  produces<float>(m_branch).setBranchAlias(m_alias);
  m_input = iConfig.getParameter<edm::InputTag>("input");
}

void EnergyDensityMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<float> rho(new float);
  edm::Handle<double> rhoH;
  iEvent.getByLabel( m_input , rhoH);
  *rho = *rhoH; 
  iEvent.put(rho, m_branch);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EnergyDensityMaker);





  
