#ifndef NTUPLEMAKER_LUMINOSITYMAKER_H
#define NTUPLEMAKER_LUMINOSITYMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// MET Filters are stored as triggers in PAT, so need some trigger headers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"



//
class MetFilterMaker : public edm::stream::EDProducer<> {

public:
     
    explicit MetFilterMaker (const edm::ParameterSet&);
    ~MetFilterMaker();

private:
     
    virtual void beginJob ();
    virtual void endJob   ();
    virtual void produce  ( edm::Event&, const edm::EventSetup& );

    std::string   aliasprefix_;
    std::string   branchprefix_;
 
    // std::string   processName_;
    edm::InputTag filtersInputTag_;
    edm::EDGetTokenT<edm::TriggerResults> filtersTokenRECO;
    edm::EDGetTokenT<edm::TriggerResults> filtersTokenPAT;
    edm::EDGetTokenT< bool >ecalBadCalibFilterUpdate_token ;
    bool doEcalFilterUpdate;

    edm::Handle<edm::TriggerResults> metFilterResultsH_;

};

#endif
