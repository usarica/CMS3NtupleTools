// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      SDFilter
// 
/**\class SDFilter SDFilter.h CMS2/NtupleMaker/interface/SDFilter.h

   Description: generic filter for cms2

   Implementation:
   - get list of names of momentum vectors as input
   - event passes if any of these vectors have pt larger than configured cut

*/
//
// Original Author:  Ingo Bloch
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: SDFilter.h,v 1.2 2010/06/15 10:08:36 fgolf Exp $
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
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"
//#include "Math/LorentzVector.h"
//
// class decleration
//

class SDFilter : public edm::EDFilter {
public:
  
    

     explicit SDFilter (const edm::ParameterSet&);
     ~SDFilter();
  
private:
     virtual void beginJob() ;
     virtual bool filter(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
   
     // ----------member data ---------------------------
     edm::InputTag elsInputTag;
     edm::InputTag musInputTag;
     edm::InputTag pfjetsInputTag;
     edm::InputTag photonInputTag;
     edm::InputTag metInputTag;
     edm::InputTag tcmetInputTag;
     edm::InputTag pfmetInputTag;

     double elsPt;
     double musPt;
     double photonPt;
     double pfjetPt;
     double metPt;
     double tcmetPt;
     double pfmetPt;
};


#endif
