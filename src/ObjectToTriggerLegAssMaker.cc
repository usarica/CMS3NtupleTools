// -*- C++ -*-
//
// Package:    ObjectToTriggerLegAssMaker
// Class:      ObjectToTriggerLegAssMaker
// 
/**\class ObjectToTriggerLegAssMaker ObjectToTriggerLegAssMaker.cc CMS3/NtupleMaker/src/ObjectToTriggerLegAssMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ObjectToTriggerLegAssMaker.cc,v 1.4 2012/08/24 18:58:55 fgolf Exp $
// This code was written by DLE
//
//

// system include files
#include <memory>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CMS3/NtupleMaker/interface/ObjectToTriggerLegAssMaker.h"

#include "TRegexp.h"
#include "TPRegexp.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

ObjectToTriggerLegAssMaker::ObjectToTriggerLegAssMaker(const edm::ParameterSet& iConfig) {

    // get configuration from event
    cone_               = iConfig.getUntrackedParameter<double>("cone");
    objectToken = consumes<std::vector<LorentzVector> >(iConfig.getUntrackedParameter<edm::InputTag>("objectInputTag"));
    triggers_           = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("triggers");
    processName_        = iConfig.getUntrackedParameter<std::string>("processName");

    triggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults",       "", processName_));

    triggerPrescaleToken= consumes<pat::PackedTriggerPrescales>(iConfig.getUntrackedParameter<std::string>("triggerPrescaleInputTag"));
    triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getUntrackedParameter<std::string>("triggerObjectsName"));

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
        std::string::size_type k = 0;
        while ((k = branchName.find('_', k)) != branchName.npos) branchName.erase(k, 1);
        branchNames_.push_back(branchName);
 
        // set up the branches
        produces<std::vector<unsigned int> >     (branchNames_[i]).setBranchAlias            (branchPrefix + "_" + trigName);
        produces<unsigned int>                   (branchNames_[i]+"version").setBranchAlias  (branchPrefix + "_" + trigName + "_version");

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
    iEvent.getByToken(objectToken, obj_p4_h);  
    if( !obj_p4_h.isValid() ) {
      throw cms::Exception("ObjectToTriggerLegAssMaker::produce: error getting obj_p4_h from Event!");
    }

    //
    // get trigger information
    //

//AOD    // trigger event handle
//AOD    edm::Handle<trigger::TriggerEvent> triggerEvent_h_;
//AOD    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEvent_h_);
//AOD    if (!triggerEvent_h_.isValid()) {
//AOD        throw cms::Exception("[ObjectToTriggerLegAssMaker][produce] Error getting TriggerEvent product from Event!");
//AOD    }
//AOD    triggerEvent_ = triggerEvent_h_.product();
//AOD
//AOD    // set process name if not set
//AOD    if (processName_ == "") {
//AOD        processName_ = triggerEvent_h_.provenance()->processName();
//AOD        bool changed(true);
//AOD        hltConfig_.init(iEvent.getRun(), iSetup, processName_, changed);
//AOD    }

    iEvent.getByToken(triggerResultsToken, triggerResultsH_);
    if (! triggerResultsH_.isValid())
      throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
    triggerNames_ = iEvent.triggerNames(*triggerResultsH_);
    // find versions for the triggers requested
    getTriggerVersions(triggers_, triggerVersions_);

    iEvent.getByToken( triggerPrescaleToken, triggerPrescalesH_);
    if (! triggerPrescalesH_.isValid())
      throw cms::Exception("HLTMaker::produce: error getting PackedTriggerPrescales product from Event!");

    //
    //  do the matching
    //

    // online objects for all triggers
//AOD    const trigger::TriggerObjectCollection &allObjects = triggerEvent_->getObjects();
    
    iEvent.getByToken(triggerObjectsToken, triggerObjectStandAlonesH_);
    if (! triggerObjectStandAlonesH_.isValid())
      throw cms::Exception("HLTMaker::produce: error getting TriggerObjectsStandAlone product from Event!");
    const pat::TriggerObjectStandAloneCollection * allObjects = triggerObjectStandAlonesH_.product();

    /*
    // loop on triggers
    std::vector<LorentzVector>::const_iterator obj_it = obj_p4_h->begin();
    for (unsigned int t = 0; t < triggers_.size(); ++t) {

        // loop on each object and match to this trigger
        // store zero if no match, otherwise prescale value of trigger
        // indexed like the collection of objects
        std::vector<unsigned int> triggerPrescales;
        for (obj_it = obj_p4_h->begin(); obj_it != obj_p4_h->end(); ++obj_it) {
            triggerPrescales.push_back(matchTriggerObject(iEvent, iSetup,
		  triggers_[t].label(), triggers_[t].instance(), t, allObjects, *obj_it));
        }

        // push back the results for this trigger
        prescales.push_back(triggerPrescales);

    }
    */

    for (unsigned int t = 0; t < triggers_.size(); ++t) {
      prescales.push_back(matchTriggerObject(iEvent, iSetup, triggers_[t].label(), triggers_[t].instance(), t, allObjects, obj_p4_h));
    }

    //
    // put the results in the event
    //

    for (unsigned int i = 0; i < triggers_.size(); ++i) {

        // get name of this trigger
        std::string trigName = triggers_[i].process();

        // create auto ptrs to hold matches (prescales) of this trigger to each electron
        // and version number of this trigger
        std::auto_ptr<std::vector<unsigned int> >   autoptr_prescales   (new std::vector<unsigned int> (prescales[i].begin(), prescales[i].end()));
        std::auto_ptr<unsigned int>                 autoptr_version     (new unsigned int (triggerVersions_[i]));

        // put into the event
        iEvent.put(autoptr_prescales,           branchNames_[i]);
        iEvent.put(autoptr_version,             branchNames_[i]+"version");
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
ObjectToTriggerLegAssMaker::beginRun(const edm::Run &r, edm::EventSetup const &c)
{

    //
    // configure HLT
    //

    // init HLT config
//AOD    if (processName_ != "") {
//AOD        bool changed(true);
//AOD        if (!hltConfig_.init(r, c, processName_, changed)) {
//AOD            throw cms::Exception("[ObjectToTriggerLegAssMaker][beginRun] Config extraction failure with process name " + processName_);
//AOD        }
//AOD    }


}

// ------------ method called when ending the processing of a run  ------------
    void
ObjectToTriggerLegAssMaker::endRun(edm::Run&, edm::EventSetup const&)
{
}

std::vector<unsigned int> ObjectToTriggerLegAssMaker::matchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
    const std::string triggerName, const std::string filterName, unsigned int triggerIndex,
    const  pat::TriggerObjectStandAloneCollection* allObjects,
    const edm::Handle<std::vector<LorentzVector> > &offlineObjects)
{
  std::vector<unsigned int> triggerPrescales;

  unsigned int prescale = 0;

  std::unordered_map<int, unsigned int> offlineObjectsPrescales; // map object index to prescale (unordered has O(1) lookup)
  for (std::vector<LorentzVector>::const_iterator obj_it = offlineObjects->begin(); obj_it != offlineObjects->end(); ++obj_it) {
    offlineObjectsPrescales[std::distance(offlineObjects->begin(), obj_it)] = 0;
  }

  // loop over trigger objects
  pat::TriggerObjectStandAlone TO;
  for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
    TO = triggerObjectStandAlonesH_->at(i);
    TO.unpackPathNames( triggerNames_ );
    if ( TO.hasPathName(triggerName, false) ) { 
      if ( (filterName == "" && TO.hasPathName(triggerName, true) ) || TO.hasFilterLabel(filterName) ) { 

        // loop over leptons
        for (std::vector<LorentzVector>::const_iterator obj_it = offlineObjects->begin(); obj_it != offlineObjects->end(); ++obj_it) {
          if(offlineObjectsPrescales[std::distance(offlineObjects->begin(), obj_it)] != 0) continue; // if we already matched this object, skip it in the subsequent loops
          if (deltaR(TO.eta(), TO.phi(), (*obj_it).eta(), (*obj_it).phi()) < cone_) {
            prescale = triggerPrescalesH_.isValid() ? triggerPrescalesH_->getPrescaleForIndex(triggerIndex) : -1;
            offlineObjectsPrescales[std::distance(offlineObjects->begin(), obj_it)] = prescale;
          }
        }

      }
    }
  }

  for (std::vector<LorentzVector>::const_iterator obj_it = offlineObjects->begin(); obj_it != offlineObjects->end(); ++obj_it) {
    triggerPrescales.push_back(offlineObjectsPrescales[std::distance(offlineObjects->begin(), obj_it)]);
  }
  return triggerPrescales;

}

/*
unsigned int ObjectToTriggerLegAssMaker::matchTriggerObject(const edm::Event &iEvent, const edm::EventSetup &iSetup,
        const std::string triggerName, const std::string filterName, unsigned int triggerIndex,
	//const trigger::TriggerObjectCollection &allObjects,
	const  pat::TriggerObjectStandAloneCollection* allObjects,
        const LorentzVector &offlineObject)
{

    unsigned int prescale = 0;

    // For miniAOD, don't loop over triggers. 
    // Loop over triggerObjects, and find one that matches criteria
    for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
      pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
      TO.unpackPathNames( triggerNames_ );
      // TO.hasPathName(triggerName, false) : TO belongs to any of the filters on this path
      // TO.hasPathName(triggerName, true ) : TO belongs to the last EDFilter on this path
      if ( TO.hasPathName(triggerName, false) ) { 
	//std::cout<<"TriggerObject belongs to path "<<triggerName<<std::endl;
	if ( (filterName == "" && TO.hasPathName(triggerName, true) ) || TO.hasFilterLabel(filterName) ) { 
	  //std::cout<<"... and to filter: "<<filterName<<". If not specified, belongs to last EDFilter of this path."<<std::endl;
	  //std::cout<<"Trigger object has eta/phi "<<TO.eta()<<"/"<<TO.phi()<<". Offline object has eta/phi "<<offlineObject.eta()<<"/"<< offlineObject.phi() <<std::endl;
	  if (deltaR(TO.eta(), TO.phi(), offlineObject.eta(), offlineObject.phi()) < cone_) {
	    prescale = triggerPrescalesH_.isValid() ? triggerPrescalesH_->getPrescaleForIndex(triggerIndex) : -1;
	    //std::cout<<"Match!! Prescale is "<<prescale<<std::endl;
	    return prescale;
	  }
	}
      }
    }
    // end of miniAOD


    // loop on triggers
//AOD    for (unsigned int i = 0; i < hltConfig_.size(); i++) {
//AOD        // get name of ith trigger
//AOD        TString hltTrigName(hltConfig_.triggerName(i));
//AOD        hltTrigName.ToLower();
//AOD        
//AOD        // pattern to match
//AOD        TString pattern(triggerName);
//AOD        pattern.ToLower();
//AOD        
//AOD        // match pattern
//AOD        TRegexp reg(Form("%s", pattern.Data()), true);
//AOD        
//AOD        // if trigger matches
//AOD        // then look for the objects corresponding to
//AOD        // the specified filter name
//AOD        if (hltTrigName.Index(reg) >= 0) {
//AOD
//AOD            edm::InputTag filterNameTag(filterName, "", processName_);
//AOD            if (filterName == "") {
//AOD                const std::vector<std::string> &modules = hltConfig_.saveTagsModules(i);
//AOD                filterNameTag = edm::InputTag(modules.back(), "", processName_);
//AOD            }
//AOD            else {
//AOD                filterNameTag = edm::InputTag(filterName, "", processName_);
//AOD            }
//AOD            
//AOD            size_t filterIndex = triggerEvent_->filterIndex(filterNameTag);
//AOD            if (filterIndex < triggerEvent_->sizeFilters()) {
//AOD                const trigger::Keys &keys = triggerEvent_->filterKeys(filterIndex);
//AOD                for (size_t j = 0; j < keys.size(); j++) {
//AOD                    trigger::TriggerObject foundObject = allObjects[keys[j]];
//AOD                    if (deltaR(foundObject.eta(), foundObject.phi(), offlineObject.eta(), offlineObject.phi()) < cone_) {
//AOD                         prescale = hltConfig_.prescaleValue(iEvent, iSetup, hltConfig_.triggerName(i));
//AOD                         return prescale;
//AOD                    }
//AOD                }
//AOD            }
//AOD        }
//AOD    }
   
    assert(prescale == 0);
    return prescale;

}
*/

void ObjectToTriggerLegAssMaker::getTriggerVersions(const std::vector<edm::InputTag> &trigNames, 
        std::vector<unsigned int> &versions)
{

    TPRegexp re("._v(.*)");

    // loop on trigger names
    for (unsigned int t = 0; t < trigNames.size(); ++t) {
        TString trigName = trigNames[t].label();

          // default value of 9999 if we don't find a match
          versions[t] = 9999;
          // loop on triggers in menu
          unsigned int nTriggers = triggerResultsH_->size();
          for(unsigned int i = 0; i < nTriggers; ++i) {

            // get name of ith trigger
            TString hltTrigName(triggerNames_.triggerName(i));
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
                break; // if we found a match, break
              } else {
                versions[t] = 0;
              }
              delete substrArr;
            }
        }

    }

}


//define this as a plug-in
DEFINE_FWK_MODULE(ObjectToTriggerLegAssMaker);

