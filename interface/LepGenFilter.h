// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      LepGenFilter
// 
/**\class LepGenFilter LepGenFilter.h CMS2/NtupleMaker/interface/LepGenFilter.h

Description: generic filter for cms2

Implementation:
- get list of names of momentum vectors as input
- event passes if any of these vectors have pt larger than configured cut

*/
//
// Original Author:  Ingo Bloch
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: LepGenFilter.h,v 1.2 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef CMS2_LEPGENFILTER_H
#define CMS2_LEPGENFILTER_H

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

//
// class decleration
//

class LepGenFilter : public edm::EDFilter {
public:
  
    

  explicit LepGenFilter (const edm::ParameterSet&);
  ~LepGenFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  int nGenLepsRequired_;

      
};


#endif
