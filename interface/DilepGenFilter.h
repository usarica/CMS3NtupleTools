// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      DilepGenFilter
// 
/**\class DilepGenFilter DilepGenFilter.h CMS2/NtupleMaker/interface/DilepGenFilter.h

Description: generic filter for cms2

Implementation:
- get list of names of momentum vectors as input
- event passes if any of these vectors have pt larger than configured cut

*/
//
// Original Author:  Ingo Bloch
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: DilepGenFilter.h,v 1.1 2009/06/03 00:13:45 kalavase Exp $
//
//
#ifndef CMS2_ASKIMFILTER_H
#define CMS2_ASKIMFILTER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"
//#include "Math/LorentzVector.h"
//
// class decleration
//

class DilepGenFilter : public edm::EDFilter {
public:
  
    

  explicit DilepGenFilter (const edm::ParameterSet&);
  ~DilepGenFilter();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
   
      
};


#endif
