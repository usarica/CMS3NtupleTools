// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      L1Maker
// 
/**\class L1Maker.cc CMS2/NtupleMaker/src/L1Maker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: L1Maker.h,v 1.5 2009/12/05 20:38:11 jribnik Exp $
//
//
#ifndef NTUPLEMAKER_L1DIGIMAKER_H
#define NTUPLEMAKER_L1DIGIMAKER_H

// system include files

// user include files
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TString.h"

class L1Maker : public edm::EDProducer {
    public:
        explicit L1Maker (const edm::ParameterSet&);
        ~L1Maker() {}

    private:
        virtual void beginJob(const edm::EventSetup&) {}
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() {}

        void fillL1Info(
                unsigned int&, unsigned int&, unsigned int&, unsigned int&,
                std::vector<TString>&,
                const L1GtTriggerMenu* menu, const DecisionWord &dWord);
        void fillL1TechnicalInfo(unsigned int&, unsigned int&, const DecisionWord&);

        bool fillL1Particles_;
        std::string l1ParticlesProcessName_;
};

#endif
