// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TheFilter
// 
/**\class TheFilter TheFilter.h CMS2/NtupleMaker/interface/TheFilter.h

Description: generic filter for cms2

Implementation:
- get list of names of vectors as input
- event passes if any of these vectors have non-zero size

*/
//
// Original Author:  Frank Wuerthwein
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: TheFilter.h,v 1.3 2010/06/15 10:08:36 fgolf Exp $
//
//
#ifndef CMS2_THEFILTER_H
#define CMS2_THEFILTER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class decleration
//

class TheFilter : public edm::EDFilter {
public:
  
    

  explicit TheFilter (const edm::ParameterSet&);
  ~TheFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
   
  // ----------member data ---------------------------
  std::vector<edm::InputTag> filterExpressions;
  
};


#endif
