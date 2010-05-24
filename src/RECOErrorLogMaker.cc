//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/RECOErrorLogMaker.h" 
#include "FWCore/MessageLogger/interface/ErrorSummaryEntry.h"
#include "FWCore/MessageLogger/interface/ELseverityLevel.h"
#include "TString.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

RECOErrorLogMaker::RECOErrorLogMaker(const edm::ParameterSet& iConfig) {

  
  produces<vector<TString> >    ("evterrSeverity"        ).setBranchAlias("evt_errSeverity"      );
  produces<vector<TString> >    ("evterrCategory"        ).setBranchAlias("evt_errCategory"      );
  produces<vector<TString> >    ("evterrModule"          ).setBranchAlias("evt_errModule"        );

  
  errorSummaryCollInputTag_ = iConfig.getParameter<InputTag>("errorSummaryCollInputTag");
  minSeverity_              = iConfig.getParameter<string>  ("minSeverity"); 

  if(minSeverity_ != "error" && minSeverity_ != "warning")
    throw cms::Exception("RECOErrorLogMaker::RECOErrorLogMaker(). Unsupported minSeverityLevel. \"minSeverity\" config parameter hould be either \"warning\" or \"error\"");
}

RECOErrorLogMaker::~RECOErrorLogMaker()
{
}

void  RECOErrorLogMaker::beginJob()
{
}

void RECOErrorLogMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void RECOErrorLogMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<TString> >    evt_errSeverity       (new vector<TString>   );
  auto_ptr<vector<TString> >    evt_errCategory       (new vector<TString>   );
  auto_ptr<vector<TString> >    evt_errModule         (new vector<TString>   );
  
  edm::Handle<vector<ErrorSummaryEntry> > errorSummary_h;
  iEvent.getByLabel(errorSummaryCollInputTag_, errorSummary_h);
  const vector<ErrorSummaryEntry> *v_errors = errorSummary_h.product();



  for(vector<ErrorSummaryEntry>::const_iterator it = v_errors->begin();
      it != v_errors->end(); it++) {    
    
    if(minSeverity_ == "warning" && it->severity.getLevel() >= ELseverityLevel::ELsev_warning) {
      evt_errSeverity->push_back(it->severity.getName());
      evt_errCategory->push_back(it->category);
      evt_errModule->push_back(it->module);
    }
    if(minSeverity_ == "error" && it->severity.getLevel() >= ELseverityLevel::ELsev_error) {
      evt_errSeverity->push_back(it->severity.getName());
      evt_errCategory->push_back(it->category);
      evt_errModule->push_back(it->module);
    }
  }

  iEvent.put(evt_errSeverity,   "evterrSeverity"	);
  iEvent.put(evt_errCategory,   "evterrCategory"	);
  iEvent.put(evt_errModule,	"evterrModule"		);
  

}

//define this as a plug-in
DEFINE_FWK_MODULE(RECOErrorLogMaker);
