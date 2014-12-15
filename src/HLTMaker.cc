#include "CMS3/NtupleMaker/interface/HLTMaker.h"
#include <map>
//#include <fstream>
#include "TBits.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace edm;
using namespace reco;
using namespace std;

/*
void PrintTriggerObjectInfo( ofstream& outfile, int id, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4 ){

  outfile.setf( ios::fixed, ios::floatfield );
  outfile << "\t\t\t\t\t"
          <<  id << "\t"
          << "( "
          << setprecision(2) << setw(4) << setfill(' ') << p4.pt() << ", "
          << setprecision(2) << setw(4) << setfill(' ') << p4.eta() << ", "
          << setprecision(2) << setw(4) << setfill(' ') << p4.phi()
          << " )"
          << endl << endl;
  return;
}
*/

HLTMaker::HLTMaker(const edm::ParameterSet& iConfig)
{
  processName_        = iConfig.getUntrackedParameter<string>         ("processName"       );
  fillTriggerObjects_ = iConfig.getUntrackedParameter<bool>           ("fillTriggerObjects");
  prunedTriggerNames_ = iConfig.getUntrackedParameter<vector<string> >("prunedTriggerNames");
  aliasprefix_        = iConfig.getUntrackedParameter<string>         ("aliasPrefix"       );
  processNamePrefix_  = TString(aliasprefix_); //just easier this way....instead of replace processNamePrefix_ everywhere
  triggerObjectsName_ = iConfig.getUntrackedParameter<string>         ("triggerObjectsName");

  produces<TBits>                           (Form("%sbits"      ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits"       ,processNamePrefix_.Data()));
  produces<vector<TString> >                (Form("%strigNames" ,processNamePrefix_.Data())).setBranchAlias(Form("%s_trigNames"  ,processNamePrefix_.Data()));
  produces<vector<unsigned int> >           (Form("%sprescales" ,processNamePrefix_.Data())).setBranchAlias(Form("%s_prescales"  ,processNamePrefix_.Data()));
  produces<vector<vector<int> > >           (Form("%strigObjsid",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_id",processNamePrefix_.Data()));

  produces<vector<vector<LorentzVector> > > (Form("%strigObjsp4",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_p4",processNamePrefix_.Data()));
}

void HLTMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

  // In the case that we are choosing the process name
  // automatically, i.e. the processName_ parameter is
  // an empty string, we can't init  HLTConfigProvider
  // until after we've determined the process name. So
  // don't init here until after we've set processName_
  // in the produce method and init there once and only
  // once. Sounds scary, it is kinda!
  // HLT config _should no longer_ change within runs :)
//AOD  if (processName_ != "") {
//AOD  bool changed(true);
//AOD  if (hltConfig_.init(iRun,iSetup,"*",changed)) {
//AOD  } else 
//AOD    throw cms::Exception("HLTMaker::beginRun: config extraction failure with process name " + processName_);
//AOD  }
}

void HLTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // If the process name is not specified retrieve the  latest
  // TriggerEvent object and the corresponding TriggerResults.
  // We should only have to do this once though, the next time
  // produce is called processName_ should be set.
//AOD  if (processName_ == "") {
//AOD    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEventH_);
//AOD    if (! triggerEventH_.isValid()  )
//AOD      throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
//AOD    // This line is important as it makes sure it is never called
//AOD    // again! A self-terminating code snippet...
//AOD    processName_ = triggerEventH_.provenance()->processName();
//AOD    // This is the once and only once bit described in beginRun
//AOD    bool changed(true);
//AOD    if (hltConfig_.init(iEvent.getRun(),iSetup,processName_,changed)) {
//AOD    } else 
//AOD      throw cms::Exception("HLTMaker::produce: config extraction failure with process name " + processName_);
//AOD  } else {
//AOD    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", processName_), triggerEventH_  );
//AOD    if (! triggerEventH_.isValid()  )
//AOD      throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
//AOD  }

// Now using a single processName_ (set to "HLT" in the configuration file). Is this OK? Do we need the flexibility we had before?
  iEvent.getByLabel(edm::InputTag("TriggerResults",       "", processName_), triggerResultsH_);
  if (! triggerResultsH_.isValid())
    throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
  
  triggerNames_ = iEvent.triggerNames(*triggerResultsH_); // Does this have to be done for every event?

  iEvent.getByLabel(triggerObjectsName_, triggerObjectStandAlonesH_);
  if (! triggerObjectStandAlonesH_.isValid())
    throw cms::Exception("HLTMaker::produce: error getting TriggerObjectsStandAlone product from Event!");

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH_; 
  iEvent.getByLabel( "patTrigger", triggerPrescalesH_);
  if (! triggerPrescalesH_.isValid())
    throw cms::Exception("HLTMaker::produce: error getting PackedTriggerPrescales product from Event!");

////// Printing miniAOD content...
//  for (unsigned int i=0; i<triggerNames_.size(); i++) {
//    std::cout << "triggerNames= " << triggerNames_.triggerName(i) << std::endl;
//  }
//  if ( triggerObjectStandAlonesH_.isValid()) cout<<"Got triggerObjectStandAlonesHandle with size "<<triggerObjectStandAlonesH_->size()<<endl;
//  else cout<<"Couldn't find triggerObjectStandAlonesH"<<endl;
//  for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
//    pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
//    TO.unpackPathNames( triggerNames_ );
//    cout<<"Trigger Object "<<i<<"has pt, eta, phi, id = "<< TO.pt() <<" "<<TO.eta()<<" "<<TO.phi()<<" and hasPathLastFilterAccepted() "<<TO.hasPathLastFilterAccepted()<<endl;
//
//    std::vector< std::string > path_namesPASS = TO.pathNames(true); //make sure they passed!
//    cout<<"Trigger Object "<<i<<" passed "<< path_namesPASS.size() <<" pathNames(true): ";
//    for (uint j  = 0; j < path_namesPASS.size(); j++) cout<<path_namesPASS[j]<<" ";
//    cout<<endl;
//
//   std::vector< std::string > filter_labels = TO.filterLabels();
//   cout<<"Trigger Object "<<i<<" has "<< filter_labels.size() <<" filter_labels: ";
//   for (uint j  = 0; j < filter_labels.size(); j++) cout<<filter_labels[j]<<" ";
//   cout<<endl;
//   
//   cout<<"Trigger Object "<<i<<" has collection() "<<TO.collection()<<endl;
//
//  }
////// END printing miniAOD content

//  // sanity check
//  assert(triggerResultsH_->size()==hltConfig_.size());

  auto_ptr<vector<unsigned int> > prescales (new vector<unsigned int>);

  unsigned int nTriggers = triggerResultsH_->size();
  //if (nTriggers > 768) throw cms::Exception( Form("HLTMaker::produce: number of HLT trigger variables must be increased! ( %d > 768 )", nTriggers) );

  auto_ptr<TBits>                           bits      (new TBits(nTriggers));
  auto_ptr<vector<TString> >                trigNames (new vector<TString>);
  auto_ptr<vector<vector<int> > >           trigObjsid(new vector<vector<int> >);
  auto_ptr<vector<vector<LorentzVector> > > trigObjsp4(new vector<vector<LorentzVector> >);
  trigNames ->reserve(nTriggers);
  trigObjsid->reserve(nTriggers);
  trigObjsp4->reserve(nTriggers);

  for(unsigned int i = 0; i < nTriggers; ++i)
    {
      // Create now because must exist regardless
      // of the accept
      vector<LorentzVector> p4V;
      vector<int> idV;

      // What is your name?
      const string& name = triggerNames_.triggerName(i);
      trigNames->push_back(name);
	
      //What is your prescale?
      prescales->push_back( triggerPrescalesH_.isValid() ? triggerPrescalesH_->getPrescaleForIndex(i) : -1 );
	
	
      // Passed... F+
      if (triggerResultsH_->accept(i)) {
          bits->SetBitNumber(i);

          // Collect desired trigger objects 
          if (fillTriggerObjects_ && doPruneTriggerName(name))
              fillTriggerObjectInfo(i, idV, p4V);
      }

      trigObjsid->push_back(idV);
      trigObjsp4->push_back(p4V);
    }

  // strip upper zeros
  bits->Compact();
  iEvent.put(bits,       Form("%sbits",       processNamePrefix_.Data() ) );

  iEvent.put(prescales,  Form("%sprescales",  processNamePrefix_.Data() ) );
  iEvent.put(trigNames , Form("%strigNames" , processNamePrefix_.Data() ) );
  iEvent.put(trigObjsid, Form("%strigObjsid", processNamePrefix_.Data() ) );
  iEvent.put(trigObjsp4, Form("%strigObjsp4", processNamePrefix_.Data() ) );
}

bool HLTMaker::doPruneTriggerName(const string& name) const
{
  for(unsigned int i = 0; i < prunedTriggerNames_.size(); ++i) {
    // uses wildcard matching like on the command line, not
    // straight up regexp
    TString pattern(prunedTriggerNames_[i]);
    pattern.ToLower();
    TRegexp reg(Form("%s", pattern.Data()), true);
    TString sname(name);
    sname.ToLower();
    if (sname.Index(reg) >= 0)
      return true;
  }
  return false;
}

void HLTMaker::fillTriggerObjectInfo(unsigned int triggerIndex, vector<int>& idV, vector<LorentzVector>& p4V) const {

  // Triggers from miniAOD. 
  // This seems slow, since it does many full loops over the trigger objects. 
  // Would probably be faster to only loop once, and then do multiple loops for the trigger names

  // 1. Get trigger name
  const string& name = triggerNames_.triggerName(triggerIndex);

  // 2. Loop over StandAloneTriggerObjects 
  for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
    pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
    TO.unpackPathNames( triggerNames_ );
    // 3. Check hasPathName( triggerName, true). 
    // This makes sure that the TriggerObjects belongs to the LAST filter in a path. 
    // With respect to CMS2forAOD, we lose the ability to get objects that don't belong to the last filter. 
    // For example, if a path is (Ecut1, Ecut2, Ecut3, Mcut1, Mcut2, Mcut3) we now only get the muon from 
    // Mcut3 and no electron. Instead, in the old CMS2forAOD we would also get the electron from Ecut3.
    if ( TO.hasPathName(name, true) ) {  
      int storeID = 0;
      std::vector<int> IDs = TO.filterIds();
      if (IDs.size() == 1) storeID = IDs[0];
      else if (IDs.size() > 1) {
	// Making some arbitrary choices
	if ( IDs[0]==85 || IDs[1]==85 ) storeID = 85; // when in doubt call it jet (and not the bjet, 86)
	if ( IDs[0]==92 || IDs[1]==92 ) storeID = 92; // when in doubt call it cluster (and not Photon, 81, or Electron, 82)

      }
//      if (IDs.size() > 1) {
//	cout<<"Multiple Types for trigger object related to "<<name<<": ";
//	for (uint k = 0; k < IDs.size(); k++) cout<<IDs[k]<<" ";
//	cout<<endl;
//      }
//      if (IDs.size() == 0) { cout<<"No IDs associated to this trigger object related to "<<name; idV.push_back( 0 ); }
      idV.push_back( storeID );
      p4V.push_back( LorentzVector( TO.p4() ) );
    }
  } // End of loop over trigger objects in miniAOD
  // miniAOD NOTE: Do we need an exception for mixed trigger, as in the CMS2forAOD code? Not sure why this was there in the first place. And impossible to implement for PAT.


  // End of Triggers from miniAOD

//AOD  // trigger name
//AOD  //const string& trigName = hltConfig_.triggerName(triggerIndex);
//AOD
//AOD  // debug 
//AOD  //string fileName = Form("%s.txt", trigName.c_str() );
//AOD  //std::ofstream outfile( fileName.c_str(), ios::out | ios::app );
//AOD
//AOD  //
//AOD  const trigger::TriggerObjectCollection& triggerObjects = triggerEventH_->getObjects();  // trigger objects
//AOD  if (triggerObjects.size() == 0) return;                                                 //
//AOD  const vector<string>& moduleLabels = hltConfig_.moduleLabels(triggerIndex);             // modules on this trigger path
//AOD  const unsigned int    moduleIndex  = triggerResultsH_->index(triggerIndex);             // index (slot position) of module giving the decision of the path
//AOD  unsigned int          nFilters     = triggerEventH_->sizeFilters();                     // number of filters
//AOD
//AOD/////////////////////////////////////
//AOD// Show all stored trigger objects //
//AOD/////////////////////////////////////
//AOD
//AOD  // CMSSW_4_2x
//AOD  // We want to store the filter index corresponding to the last filter for each distinct trigger object type ( the trigger id & p4 can be accessed given the filter index )
//AOD  // For Electrons: take the filter index for the last trigger id of either 82 or 92, we don't consider 82 and 92 distinct trigger objects
//AOD  // This assumes that each filter will have the same trigger ids... this is *NOT* true for overlap triggers, PFTau triggers for instance
//AOD  // If we find that any filter has mixed trigger ids, we store the objects for the last filter only ( this is what was done before and should be right )
//AOD  // Yes this is kludged, and yes a universal solution is preferred
//AOD
//AOD  // map ( trigger id => filter index )
//AOD  map< int, unsigned int > trigObjsToStore;
//AOD
//AOD  // True if different trigger id's are found in the same filter
//AOD  bool mixedTrigIds = false;
//AOD
//AOD  // loop over trigger modules
//AOD  unsigned int lastFilterIndex = nFilters;
//AOD  for(unsigned int j = 0; j <= moduleIndex; ++j) {
//AOD
//AOD    // get module name & filter index
//AOD    const string&      moduleLabel = moduleLabels[j];
//AOD    const unsigned int filterIndex = triggerEventH_->filterIndex(InputTag(moduleLabel, "", processName_));
//AOD
//AOD    // debug
//AOD    //if( j == 0 ){
//AOD    //  outfile << endl << "--------------------------------------------------------------------------------------------------------------" << endl << endl;
//AOD    //  outfile << trigName << endl << endl;
//AOD    //  outfile << "-> Everything:" << endl << endl;
//AOD    //}
//AOD      
//AOD    // these are the filters with trigger objects filled
//AOD    if ( filterIndex < nFilters ){ 
//AOD     
//AOD      //
//AOD      lastFilterIndex = filterIndex;
//AOD   
//AOD      // get trigger objects & ids
//AOD      //unsigned int lastEleFilterIndex = 0;
//AOD      const trigger::Vids& triggerIds  = triggerEventH_->filterIds(filterIndex);
//AOD      const trigger::Keys& triggerKeys = triggerEventH_->filterKeys(filterIndex);
//AOD      assert( triggerIds.size() == triggerKeys.size() );
//AOD
//AOD      // True if a filter has hlt objects ( objects with positive trigger id )
//AOD      bool hlt = false;   
//AOD      
//AOD      // Trigger object id for filters associated with a unique trigger object id
//AOD      int filterId = -1;
//AOD
//AOD      // loop on trigger objects
//AOD      for( unsigned int k = 0; k < triggerKeys.size(); k++){
//AOD
//AOD        // trigger object id
//AOD        int id = triggerIds.at(k);
//AOD
//AOD        // trigger p4, p4id
//AOD        //const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];   
//AOD
//AOD        // hlt objects
//AOD        if( id > 0 ){
//AOD          hlt = true;                                 // True if a filter has hlt objects ( objects with positive trigger id )
//AOD          if( k == 0 ){                               // Assuming all trigger ids are the same for this filter
//AOD            filterId = id;                            // Filter id for all the objects in this filter
//AOD            if( filterId == 82 || filterId == 92 ){   // Remember the index of the last filter with all objects having a trigger id of 82 or 92 
//AOD              //lastEleFilterIndex = filterIndex;
//AOD            }
//AOD          }
//AOD          else {
//AOD            if( id != filterId ){
//AOD              mixedTrigIds = true;                    // True if different trigger id's are found in the same filter
//AOD              //cout << endl << endl << trigName << ", " << moduleLabel << ":\t ERROR in HLTMaker.cc: different trigger object ids found in the same filter... Exiting." << endl << endl;
//AOD              //exit(1);
//AOD            }
//AOD          }
//AOD        } // end if( id > 0 )
//AOD
//AOD        // debug
//AOD        //if( k == 0 ) outfile << "\t" << filterIndex << ": " << moduleLabel << endl << endl;
//AOD        //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );
//AOD
//AOD      } // end loop on trigger objects
//AOD
//AOD      //
//AOD      if( hlt ){                                        // only store hlt objects
//AOD        assert( filterId != -1 );                       // sanity
//AOD        trigObjsToStore[filterId] = filterIndex;        // Store the filter Index ( used to get trigger objects ) for each different trigger type ( filterId )
//AOD        if( filterId == 82 ) trigObjsToStore.erase(92); // If this is an electron trigger ( filterId 82 or 92 ) we only want the last one ( that is either 82 or 92 )
//AOD        if( filterId == 92 ) trigObjsToStore.erase(82); // 
//AOD      }
//AOD
//AOD    } // end if(filterIndex < nFilters)
//AOD  }   // end loop over trigger modules
//AOD 
//AOD//////////////////////////////////////////////
//AOD// Show the trigger objects we want to save //
//AOD//////////////////////////////////////////////
//AOD
//AOD  // debug
//AOD  //outfile << endl << endl << "-> What we want to store";
//AOD  //if(mixedTrigIds) outfile << " ( mixed ids... storing last filter only )";
//AOD  //outfile << ":" << endl << endl;
//AOD
//AOD  // If we find different trigger ids under any filter we only store the objects for the last filter
//AOD  if( mixedTrigIds ){
//AOD
//AOD    // get trigger objects & ids
//AOD    const trigger::Vids& triggerIds  = triggerEventH_->filterIds ( lastFilterIndex );
//AOD    const trigger::Keys& triggerKeys = triggerEventH_->filterKeys( lastFilterIndex );
//AOD    assert( triggerIds.size() == triggerKeys.size() );
//AOD
//AOD    // loop on trigger objects
//AOD    for( unsigned int k = 0; k < triggerKeys.size(); k++){
//AOD
//AOD      // trigger object id
//AOD      int id = triggerIds.at(k);                                                        
//AOD
//AOD      // trigger p4, p4id
//AOD      const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];
//AOD  
//AOD      // store trigger id, trigger p4, & trigger object id
//AOD      p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
//AOD      idV.push_back( id );
//AOD
//AOD      // debug
//AOD      //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );
//AOD    }
//AOD 
//AOD  }
//AOD  // If each filter has trigger objects with the same trigger id, we store the trigger 
//AOD  // objects for the last filter of each distinct trigger object type
//AOD  else {
//AOD
//AOD    map<int,unsigned int>::iterator it;
//AOD    for( it = trigObjsToStore.begin() ; it != trigObjsToStore.end(); it++ ){
//AOD    
//AOD      // get trigger objects & ids
//AOD      const unsigned int filterIndex   = (*it).second;
//AOD      const trigger::Vids& triggerIds  = triggerEventH_->filterIds ( filterIndex );
//AOD      const trigger::Keys& triggerKeys = triggerEventH_->filterKeys( filterIndex );
//AOD      assert( triggerIds.size() == triggerKeys.size() );
//AOD    
//AOD      // loop on trigger objects
//AOD      for( unsigned int k = 0; k < triggerKeys.size(); k++){
//AOD
//AOD        // trigger object id
//AOD        int id = triggerIds.at(k);                                                        
//AOD
//AOD        // trigger p4, p4id
//AOD        const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];
//AOD    
//AOD        // store trigger id, trigger p4, & trigger object id
//AOD        p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
//AOD        idV.push_back( id );
//AOD
//AOD        // debug
//AOD        //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );
//AOD      } 
//AOD
//AOD    }
//AOD  }

/*
/////////////////////////////////
// Show what was stored (2010) //
/////////////////////////////////

  // these are the filters with trigger objects filled
  if (lastFilterIndex < nFilters) {

    // debug
    outfile << endl << endl << "-> What is stored now:" << endl << endl;

    // get trigger objects & ids
    const trigger::Vids& triggerIds = triggerEventH_->filterIds(lastFilterIndex);
    const trigger::Keys& triggerKeys = triggerEventH_->filterKeys(lastFilterIndex);
    assert( triggerIds.size() == triggerKeys.size() );

    // loop on trigger objects
    for(unsigned int j = 0; j < triggerKeys.size(); ++j) {

      // trigger p4, p4id
      const trigger::TriggerObject& triggerObject = triggerObjects[triggerKeys[j]];

      // store trigger id, trigger p4, & trigger object id
      //p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
      //idV.push_back( triggerObject.id() );

      // debug
      PrintTriggerObjectInfo( outfile, 0, triggerObject.particle().p4() );
    }
  }
*/

} // end HLTMaker::fillTriggerObjectInfo()

//define this as a plug-in
DEFINE_FWK_MODULE(HLTMaker);
