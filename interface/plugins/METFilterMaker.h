#ifndef NTUPLEMAKER_METFILTERMAKER_H
#define NTUPLEMAKER_METFILTERMAKER_H

#include <string>
#include <vector>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"


class METFilterMaker : public edm::stream::EDProducer<>{
public:

  explicit METFilterMaker(const edm::ParameterSet&);
  ~METFilterMaker();

private:
  virtual void beginJob();
  virtual void endJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);

  std::string aliasprefix_;
  bool doEcalFilterUpdate;

  edm::InputTag filtersInputTag_;

  edm::EDGetTokenT<edm::TriggerResults> filtersTokenRECO;
  edm::EDGetTokenT<edm::TriggerResults> filtersTokenPAT;
  edm::EDGetTokenT<bool> ecalBadCalibFilterUpdate_token;

};


#endif
