#ifndef NTUPLEMAKER_ECALANOMALOUSEVENTMAKER_H
#define NTUPLEMAKER_ECALANOMALOUSEVENTMAKER_H


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class EcalAnomalousEventMaker : public edm::EDProducer {

  public:
     explicit EcalAnomalousEventMaker (const edm::ParameterSet&);
     ~EcalAnomalousEventMaker();

  private:

    virtual void produce( edm::Event&, const edm::EventSetup& );
    virtual void beginJob();
    virtual void endJob();
    std::string aliasprefix_;
    std::string branchprefix_;

};

#endif
