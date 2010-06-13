// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToPFMuAssMaker
// 
/**\class MuToPFMuAssMaker MuToPFMuAssMaker.cc CMS2/MuToPFMuAssMaker/src/MuToPFMuAssMaker.cc

   Description: make associations between muons and tracks

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToPFMuAssMaker.h,v 1.1 2010/06/13 16:00:49 fgolf Exp $
//
//
#ifndef CMS2_MUTOPFMUASSMAKER_H
#define CMS2_MUTOPFMUASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class MuToPFMuAssMaker : public edm::EDProducer {
public:
     explicit MuToPFMuAssMaker (const edm::ParameterSet&);

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
