// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CSA07ProcessFilter
// 
/**\class CSA07ProcessFilter CSA07ProcessFilter.h CMS2/NtupleMaker/interface/CSA07ProcessFilter.h

Filter on CSA07 process ID (which has been written into a cms2 ntuple)

*/
//
// Original Author:  Frank Wuerthwein
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: CSA07ProcessFilter.h,v 1.1 2008/08/12 05:30:34 jmuelmen Exp $
//
//
#ifndef CMS2_CSA07ProcessFilter_H
#define CMS2_CSA07ProcessFilter_H

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
// class declaration
//

class CSA07ProcessFilter : public edm::EDFilter {
public:
     explicit CSA07ProcessFilter (const edm::ParameterSet&);
  
private:
     virtual bool filter(edm::Event&, const edm::EventSetup&);
     
     // ----------member data ---------------------------
     std::vector<int> 	processIds_;
     edm::InputTag	csa07process_;
  
};


#endif
