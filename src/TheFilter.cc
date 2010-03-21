//-*- C++ -*-
//
// Package:    TheFilter
// Class:      TheFilter
// 
/**\class TheFilter TheFilter.cc CMS2/TheFilter/src/TheFilter.cc

Description: generic filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Frank Wuerthwein
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TheFilter.cc,v 1.3 2010/03/21 02:55:37 warren Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/TheFilter.h"

using namespace std;

//
// constructors and destructor
//

TheFilter::TheFilter(const edm::ParameterSet& iConfig) {
  filterExpressions = iConfig.getParameter<std::vector<edm::InputTag> >("filterExpressions");
}


TheFilter::~TheFilter() {}

void  TheFilter::beginJob() {
}

void TheFilter::endJob() {
}


// ------------ method called to produce the data  ------------
bool TheFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //loop over strings in input
  for ( vector<edm::InputTag>::const_iterator i=filterExpressions.begin(), end=filterExpressions.end();
	i!=end; ++i) {

	try{
	  edm::Handle<vector<int> > theHandle;
	  iEvent.getByLabel(*i, theHandle);
	  if ( theHandle->size() != 0 ) return true ;
	}
	catch(...){
	  //cout << "not an int" << endl; //just skip if bad
	}
	try{
	  edm::Handle<vector<float> > theHandle;
	  iEvent.getByLabel(*i, theHandle);
	  if ( theHandle->size() != 0 ) return true ;
	}
	catch(...){
	  //cout << "not a float" << endl; //just skip if bad
	}
  }

  //get handle for given input string


  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TheFilter);





  
