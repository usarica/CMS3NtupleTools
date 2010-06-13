// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFMuToMuAssMaker
// 
/**\class PFMuToMuAssMaker PFMuToMuAssMaker.cc CMS2/PFMuToMuAssMaker/src/PFMuToMuAssMaker.cc

   Description: make associations between muons and tracks

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMuToMuAssMaker.h,v 1.1 2010/06/13 16:00:49 fgolf Exp $
//
//
#ifndef CMS2_PFMUTOMUASSMAKER_H
#define CMS2_PFMUTOMUASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class PFMuToMuAssMaker : public edm::EDProducer {
public:
     explicit PFMuToMuAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
     double m_minDR_;
     std::string aliasprefix_;
     edm::InputTag musInputTag_;
     edm::InputTag pfmusInputTag_;
};


#endif
