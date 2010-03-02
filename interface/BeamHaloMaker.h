// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BeamHaloMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/BeamHaloMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: BeamHaloMaker.h,v 1.3 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_BEAMSPOTMAKER_H
#define NTUPLEMAKER_BEAMSPOTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TString.h"

//
// class decleration
//

class BeamHaloMaker : public edm::EDProducer {
public:
     explicit BeamHaloMaker (const edm::ParameterSet&);
      ~BeamHaloMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     std::vector<TString> convertToVectorTString(const std::vector<char>);

   edm::InputTag beamHaloInputTag;
    
};


#endif
