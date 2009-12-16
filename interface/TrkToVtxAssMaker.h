// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrkToVtxAssMaker
// 
/**\class TrkToVtxAssMaker TrkToVtxAssMaker.cc CMS2/TrkToVtxAssMaker/src/TrkToVtxAssMaker.cc

 Description: calculate the d0 and d0Error of the track wrt the highest-sumpt vertex

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrkToVtxAssMaker.h,v 1.1 2009/12/16 09:52:03 jmuelmen Exp $
//
//
#ifndef CMS2_TRKTOVTXASSMAKER_H
#define CMS2_TRKTOVTXASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//
// class declaration
//

class TrkToVtxAssMaker : public edm::EDProducer {
public:
     explicit TrkToVtxAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
     edm::InputTag m_vtxInputTag;
     edm::InputTag m_trksInputTag;
};

#endif
