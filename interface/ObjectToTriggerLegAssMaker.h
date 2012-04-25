// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ObjectToTriggerLegAssMaker
// 
/**\class ObjectToTriggerLegAssMaker ObjectToTriggerLegAssMaker.cc CMS2/ObjectToTriggerLegAssMaker/src/ObjectToTriggerLegAssMaker.cc

Description: make associations between p4s and trigger objects
of a specified filter...

Implementation:
<Notes on implementation>
 */
//
//
//
#ifndef CMS2_LEPTONTOTRIGGERLEGASSMAKER_H
#define CMS2_LEPTONTOTRIGGERLEGASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include <vector>
#include <string>

//
// class declaration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class ObjectToTriggerLegAssMaker : public edm::EDProducer {
    public:
        explicit ObjectToTriggerLegAssMaker (const edm::ParameterSet&);

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);

        // match an offline p4 to a trigger object
        unsigned int matchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
            const std::string triggerName, const std::string filterName,
            const trigger::TriggerObjectCollection &allObjects,
            const LorentzVector &offlineObject);

        // get version of triggers
        void getTriggerVersions(const std::vector<edm::InputTag> &trigNames, 
            std::vector<unsigned int> &versions);

        // ----------member data ---------------------------

        // electrons and muons
        edm::InputTag objectInputTag_;

        // triggers to match to
        std::vector<edm::InputTag>      triggers_;    
        std::vector<unsigned int>       triggerVersions_;
        std::vector<std::string>        branchNames_;

        const trigger::TriggerEvent*    triggerEvent_;
        HLTConfigProvider               hltConfig_;
        std::string                     processName_;
        double                          cone_;

};

#endif

