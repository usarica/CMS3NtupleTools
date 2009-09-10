//-*- C++ -*-
//
// Package:    ASkimFilter
// Class:      ASkimFilter
// 
/**\class ASkimFilter ASkimFilter.cc CMS2/src/ASkimFilter.cc

Description: filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Ingo Bloch
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ASkimFilter.cc,v 1.2 2009/09/10 10:51:42 fgolf Exp $
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

#include "CMS2/NtupleMaker/interface/ASkimFilter.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;

//
// constructors and destructor
//

ASkimFilter::ASkimFilter(const edm::ParameterSet& iConfig) {
  filterExpressions = iConfig.getParameter<std::vector<edm::InputTag> >("filterExpressions");
  filterPtCut       = iConfig.getParameter<double> ("filterPtCut");
}


ASkimFilter::~ASkimFilter() {}

void  ASkimFilter::beginJob(const edm::EventSetup&) {
}

void ASkimFilter::endJob() {
}


// ------------ method called to produce the data  ------------
bool ASkimFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //  std::cout<<"ASkimFilter called. Pt cut is: "<<filterPtCut<<std::endl;

  //loop over strings in input
  for ( vector<edm::InputTag>::const_iterator i=filterExpressions.begin(), end=filterExpressions.end();
	i!=end; ++i) {
    
    // somewhat dirty - handle needs to be TLorentzVector and correspond to momentum.
    edm::Handle< vector<LorentzVector> > theHandle;
    iEvent.getByLabel(*i, theHandle);
    if ( theHandle->size() != 0 ) {
      
      for ( vector<LorentzVector>::const_iterator lepp4_it = theHandle->begin();
            lepp4_it != theHandle->end(); ++lepp4_it) {
        if( (*lepp4_it).pt() > filterPtCut ) {
          return true;
        }
      }
    }
  }

  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ASkimFilter);





  
