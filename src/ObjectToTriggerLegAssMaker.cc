// -*- C++ -*-
//
// Package:    ObjectToTriggerLegAssMaker
// Class:      ObjectToTriggerLegAssMaker
// 
/**\class ObjectToTriggerLegAssMaker ObjectToTriggerLegAssMaker.cc CMS2/NtupleMaker/src/ObjectToTriggerLegAssMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ObjectToTriggerLegAssMaker.cc,v 1.1 2012/04/25 14:06:14 dlevans Exp $
// This code was written by DLE
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CMS2/NtupleMaker/interface/ObjectToTriggerLegAssMaker.h"

#include "TRegexp.h"
#include "TPRegexp.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

ObjectToTriggerLegAssMaker::ObjectToTriggerLegAssMaker(const edm::ParameterSet& iConfig) {

    // get configuration from event
    cone_               = iConfig.getParameter<double>("cone");
    objectInputTag_     = iConfig.getParameter<edm::InputTag>("objectInputTag");
    triggers_           = iConfig.getParameter<std::vector<edm::InputTag> >("triggers");

    // get the branch prefix and remove spurious _s
    std::string aliasprefix = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchPrefix = aliasprefix;
    if (branchPrefix.find("_") != std::string::npos) branchPrefix.replace(branchPrefix.find("_"), 1, "");

    // set up muon trigger branches
    triggerVersions_.reserve(triggers_.size());
    branchNames_.reserve(triggers_.size());

    for (unsigned int i = 0; i < triggers_.size(); ++i) {

        // figure out the branch names, removing spurious _s
        std::string trigName = triggers_[i].process();
        std::string branchName = branchPrefix + trigName;
        if (branchName.find("_") != std::string::npos) branchName.replace(branchName.find("_"), 1, "");
        branchNames_.push_back(branchName);
 
        // set up the branches
        produces<std::vector<int> >     (branchNames_[i]).setBranchAlias            (branchPrefix + "_" + trigName);
        produces<int>                   (branchNames_[i]+"version").setBranchAlias  (branchPrefix + "_" + trigName + "_version");

    }

}

void ObjectToTriggerLegAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //
    // output to store
    //

    std::vector<std::vector<unsigned int> >     prescales;
    std::vector<unsigned int>                   versions;
    prescales.reserve(triggers_.size());
    versions.reserve(triggers_.size());

    //
    // get electrons and muons
    //

    edm::Handle<std::vector<LorentzVector> > obj_p4_h;
    iEvent.getByLabel(objectInputTag_, obj_p4_h);  

    //
    // get trigger information
    //

    // trigger event handle
    edm::Handle<trigger::TriggerEvent> triggerEvent_h_;
    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEvent_h_);
    if (!triggerEvent_h_.isValid()) {
        throw cms::Exception("[ObjectToTriggerLegAssMaker][produce] Error getting TriggerEvent product from Event!");
    }
    triggerEvent_ = triggerEvent_h_.product();

    // set process name if not set
    if (processName_ == "") {
        processName_ = triggerEvent_h_.provenance()->processName();
        bool changed(true);
        hltConfig_.init(iEvent.getRun(), iSetup, processName_, changed);
    }

    //
    //  do the matching
    //

    // online objects for all triggers
    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();

    // loop on muon triggers
    std::vector<LorentzVector>::const_iterator obj_it = obj_p4_h->begin();
    for (unsigned int t = 0; t < triggers_.size(); ++t) {

        // loop on each muon and match to this trigger
        // store zero if no match, otherwise prescale value of trigger
        for (obj_it = obj_p4_h->begin(); obj_it != obj_p4_h->end(); ++obj_it) {
            prescales[t].push_back(matchTriggerObject(iEvent, iSetup,
                triggers_[t].label(), triggers_[t].instance(), allObjects, *obj_it));
        }
    }

    //
    // put the results in the event
    //

    // put electron trigger branches
    for (unsigned int i = 0; i < triggers_.size(); ++i) {

        // get name of this trigger
        std::string trigName = triggers_[i].process();

        // create auto ptrs to hold matches (prescales) of this trigger to each electron
        // and version number of this trigger
        std::auto_ptr<std::vector<unsigned int> > autoptr_prescales(new std::vector<unsigned int> (prescales[i].begin(), prescales[i].end()));
        std::auto_ptr<unsigned int> autoptr_version(new unsigned int (triggerVersions_[i]));

        // put into the event
        iEvent.put(autoptr_prescales,           branchNames_[i]);
        iEvent.put(autoptr_version,             branchNames_[i]);
    }

}

// ------------ method called once each job just before starting event loop  ------------
    void 
ObjectToTriggerLegAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ObjectToTriggerLegAssMaker::endJob() {
}

// ------------ method called when starting to processes a run  ------------
    void
ObjectToTriggerLegAssMaker::beginRun(edm::Run &r, edm::EventSetup const &c)
{

    //
    // configure HLT
    //

    // init HLT config
    if (processName_ != "") {
        bool changed(true);
        if (!hltConfig_.init(r, c, processName_, changed)) {
            throw cms::Exception("[ObjectToTriggerLegAssMaker][beginRun] Config extraction failure with process name " + processName_);
        }
    }

    // find versions for the triggers requested
    getTriggerVersions(triggers_, triggerVersions_);

}

// ------------ method called when ending the processing of a run  ------------
    void
ObjectToTriggerLegAssMaker::endRun(edm::Run&, edm::EventSetup const&)
{
}

unsigned int ObjectToTriggerLegAssMaker::matchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
        const std::string triggerName, const std::string filterName,
        const trigger::TriggerObjectCollection &allObjects,
        const LorentzVector &offlineObject)
{

    unsigned int prescale = 0;

    // loop on triggers
    for (unsigned int i = 0; i < hltConfig_.size(); i++) {
        
        // get name of ith trigger
        TString hltTrigName(hltConfig_.triggerName(i));
        hltTrigName.ToLower();
        
        // pattern to match
        TString pattern(triggerName);
        pattern.ToLower();
        
        // match pattern
        TRegexp reg(Form("%s", pattern.Data()), true);
        
        // if trigger matches
        // then look for the objects corresponding to
        // the specified filter name
        if (hltTrigName.Index(reg) >= 0) {
            
            edm::InputTag filterNameTag(filterName, "", processName_);
            if (filterName == "") {
                const std::vector<std::string> &modules = hltConfig_.saveTagsModules(i);
                filterNameTag = edm::InputTag(modules.back(), "", processName_);
            }
            else {
                filterNameTag = edm::InputTag(filterName, "", processName_);
            }
            
            size_t filterIndex = triggerEvent_->filterIndex(filterNameTag);
            if (filterIndex < triggerEvent_->sizeFilters()) {
                const trigger::Keys &keys = triggerEvent_->filterKeys(filterIndex);
                for (size_t j = 0; j < keys.size(); j++) {
                    trigger::TriggerObject foundObject = allObjects[keys[j]];
                    if (deltaR(foundObject.eta(), foundObject.phi(), offlineObject.eta(), offlineObject.phi()) < cone_) {
                         prescale = hltConfig_.prescaleValue(iEvent, iSetup, hltConfig_.triggerName(i));
                         return prescale;
                    }
                }
            }
        }
    }
    
    assert(prescale == 0);
    return prescale;

}

void ObjectToTriggerLegAssMaker::getTriggerVersions(const std::vector<edm::InputTag> &trigNames, 
        std::vector<unsigned int> &versions)
{

    TPRegexp re("._v(.*)");

    // loop on trigger names
    for (unsigned int t = 0; t < trigNames.size(); ++t) {
        TString trigName = trigNames[t].label();

        // loop on triggers in menu
        for (unsigned int i = 0; i < hltConfig_.size(); i++) {

            // get name of ith trigger
            TString hltTrigName(hltConfig_.triggerName(i));
            hltTrigName.ToLower();

            // test if it matches this trigger name
            // with any version
            TString pattern(trigName);
            pattern.ToLower();
            TRegexp reg(Form("%s", pattern.Data()), true);

            // if trigger matches
            // then extract version number
            if (hltTrigName.Index(reg) >= 0) {

                TObjArray *substrArr = re.MatchS(hltTrigName);
                if (substrArr->GetLast() == 1) {
                    versions[t] = ((TObjString*)substrArr->At(1))->GetString().Atoi();
                } else {
                    versions[t] = 0;
                }

            }
        }

    }

}


//define this as a plug-in
DEFINE_FWK_MODULE(ObjectToTriggerLegAssMaker);

