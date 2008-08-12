//-*- C++ -*-
//
// Package:    CSA07ProcessFilter
// Class:      CSA07ProcessFilter
// 
/**\class CSA07ProcessFilter CSA07ProcessFilter.cc CMS2/CSA07ProcessFilter/src/CSA07ProcessFilter.cc

Description: generic filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Frank Wuerthwein
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: CSA07ProcessFilter.cc,v 1.1 2008/08/12 05:30:35 jmuelmen Exp $
//
//


// system include files
#include <memory>
#include <assert.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/CSA07ProcessFilter.h"

using namespace std;

//
// constructors and destructor
//

CSA07ProcessFilter::CSA07ProcessFilter (const edm::ParameterSet& iConfig) 
     : processIds_(iConfig.getParameter<std::vector<int> >("processIds")),
       csa07process_(iConfig.getParameter<edm::InputTag>("csa07process"))
{
     
}

// ------------ method called to produce the data  ------------
bool CSA07ProcessFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
     edm::Handle<int> evt_CSA07Process_h;
     iEvent.getByLabel(csa07process_, evt_CSA07Process_h);
     assert(processIds_.size() == 2);
     return *evt_CSA07Process_h >= processIds_[0] && *evt_CSA07Process_h <= processIds_[1];
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSA07ProcessFilter);





  
