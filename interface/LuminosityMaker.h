// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      LuminosityMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/LuminosityMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: LuminosityMaker.h,v 1.1 2011/01/26 01:33:14 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_LUMINOSITYMAKER_H
#define NTUPLEMAKER_LUMINOSITYMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class LuminosityMaker : public edm::EDProducer {
public:
     explicit LuminosityMaker (const edm::ParameterSet&);
     ~LuminosityMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     edm::InputTag lumiSummaryInputTag_;
	 std::string aliasprefix_;
};


#endif
