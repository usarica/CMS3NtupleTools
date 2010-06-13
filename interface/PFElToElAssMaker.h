// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFElToElAssMaker
// 
/**\class PFElToElAssMaker PFElToElAssMaker.cc CMS2/PFElToElAssMaker/src/PFElToElAssMaker.cc

   Description: make associations between muons and tracks

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFElToElAssMaker.h,v 1.1 2010/06/13 16:00:49 fgolf Exp $
//
//
#ifndef CMS2_PFELTOELASSMAKER_H
#define CMS2_PFELTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class PFElToElAssMaker : public edm::EDProducer {
public:
     explicit PFElToElAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
     double m_minDR_;
     std::string aliasprefix_;
     edm::InputTag elsInputTag_;
     edm::InputTag pfelsInputTag_;
};


#endif
