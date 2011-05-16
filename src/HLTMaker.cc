#include "CMS2/NtupleMaker/interface/HLTMaker.h"
#include <map>
//#include <fstream>

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

  produces<unsigned int>                    (Form("%sbits1"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits1"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits2"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits2"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits3"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits3"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits4"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits4"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits5"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits5"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits6"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits6"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits7"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits7"      ,processNamePrefix_.Data()));
  produces<unsigned int>                    (Form("%sbits8"     ,processNamePrefix_.Data())).setBranchAlias(Form("%s_bits8"      ,processNamePrefix_.Data()));
  produces<vector<TString> >                (Form("%strigNames" ,processNamePrefix_.Data())).setBranchAlias(Form("%s_trigNames"  ,processNamePrefix_.Data()));
  produces<vector<unsigned int> >           (Form("%sprescales" ,processNamePrefix_.Data())).setBranchAlias(Form("%s_prescales"  ,processNamePrefix_.Data()));
  produces<vector<vector<int> > >           (Form("%strigObjsid",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_id",processNamePrefix_.Data()));

  produces<vector<vector<LorentzVector> > > (Form("%strigObjsp4",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_p4",processNamePrefix_.Data()));
}

void HLTMaker::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{
  // In the case that we are choosing the process name
  // automatically, i.e. the processName_ parameter is
  // an empty string, we can't init  HLTConfigProvider
  // until after we've determined the process name. So
  // don't init here until after we've set processName_
  // in the produce method and init there once and only
  // once. Sounds scary, it is kinda!
  // HLT config _should no longer_ change within runs :)
  if (processName_ != "") {
    bool changed(true);
    if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    } else 
      throw cms::Exception("HLTMaker::beginRun: config extraction failure with process name " + processName_);
  }
}

void HLTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // If the process name is not specified retrieve the  latest
  // TriggerEvent object and the corresponding TriggerResults.
  // We should only have to do this once though, the next time
  // produce is called processName_ should be set.
  if (processName_ == "") {
    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEventH_);
    if (! triggerEventH_.isValid()  )
      throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
    // This line is important as it makes sure it is never called
    // again! A self-terminating code snippet...
    processName_ = triggerEventH_.provenance()->processName();
    // This is the once and only once bit described in beginRun
    bool changed(true);
    if (hltConfig_.init(iEvent.getRun(),iSetup,processName_,changed)) {
    } else 
      throw cms::Exception("HLTMaker::produce: config extraction failure with process name " + processName_);
  } else {
    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", processName_), triggerEventH_  );
    if (! triggerEventH_.isValid()  )
      throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
  }
  iEvent.getByLabel(edm::InputTag("TriggerResults",       "", processName_), triggerResultsH_);
  if (! triggerResultsH_.isValid())
    throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
  // sanity check
  assert(triggerResultsH_->size()==hltConfig_.size());

  auto_ptr<unsigned int> bits1 (new unsigned int);
  auto_ptr<unsigned int> bits2 (new unsigned int);
  auto_ptr<unsigned int> bits3 (new unsigned int);
  auto_ptr<unsigned int> bits4 (new unsigned int);
  auto_ptr<unsigned int> bits5 (new unsigned int);
  auto_ptr<unsigned int> bits6 (new unsigned int);
  auto_ptr<unsigned int> bits7 (new unsigned int);
  auto_ptr<unsigned int> bits8 (new unsigned int);
  auto_ptr<vector<unsigned int> > prescales (new vector<unsigned int>);
  *bits1 = 0;
  *bits2 = 0;
  *bits3 = 0;
  *bits4 = 0;
  *bits5 = 0;
  *bits6 = 0;
  *bits7 = 0;
  *bits8 = 0;

  unsigned int nTriggers = triggerResultsH_->size();
  if (nTriggers > 512) throw cms::Exception( Form("HLTMaker::produce: number of HLT trigger variables must be increased! ( nTriggers = %d )", nTriggers) );

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
      const string& name = hltConfig_.triggerName(i);
      trigNames->push_back(name);
	
      //What is your prescale?
      prescales->push_back(hltConfig_.prescaleValue(iEvent, iSetup, name));
	
	
      // Passed... F+
      if (triggerResultsH_->accept(i)) {
	// Encode trigger bits
	unsigned int bitmask = 1;
	if (i <= 31) {
	  bitmask <<=i;
	  *bits1 |= bitmask;
	}
	if (i >= 32 && i <= 63) {
	  bitmask <<=(i-32);
	  *bits2 |= bitmask;
	}
	if (i >= 64 && i <= 95) {
	  bitmask <<=(i-64);
	  *bits3 |= bitmask;
	}
	if (i >= 96 && i <= 127) {
	  bitmask <<=(i-96);
	  *bits4 |= bitmask;
	}
	if (i >= 128 && i <= 159) {
	  bitmask <<=(i-128);
	  *bits5 |= bitmask;
	}
	if (i >= 160 && i <= 191) {
	  bitmask <<=(i-160);
	  *bits6 |= bitmask;
	}
	if (i >= 192 && i <= 223) {
	  bitmask <<=(i-192);
	  *bits7 |= bitmask;
	}
	if (i >= 224 && i <= 255) {
	  bitmask <<=(i-224);
	  *bits8 |= bitmask;
	}

	// Collect desired trigger objects 
	if (fillTriggerObjects_ && doPruneTriggerName(name))
	  fillTriggerObjectInfo(i, idV, p4V);
      }

      trigObjsid->push_back(idV);
      trigObjsp4->push_back(p4V);
    }

  iEvent.put(bits1,      Form("%sbits1",   processNamePrefix_.Data()));
  iEvent.put(bits2,      Form("%sbits2",   processNamePrefix_.Data()));
  iEvent.put(bits3,      Form("%sbits3",   processNamePrefix_.Data()));
  iEvent.put(bits4,      Form("%sbits4",   processNamePrefix_.Data()));
  iEvent.put(bits5,      Form("%sbits5",   processNamePrefix_.Data()));
  iEvent.put(bits6,      Form("%sbits6",   processNamePrefix_.Data()));
  iEvent.put(bits7,      Form("%sbits7",   processNamePrefix_.Data()));
  iEvent.put(bits8,      Form("%sbits8",   processNamePrefix_.Data()));
  iEvent.put(prescales,  Form("%sprescales",  processNamePrefix_.Data()));
  iEvent.put(trigNames , Form("%strigNames" , processNamePrefix_.Data()));
  iEvent.put(trigObjsid, Form("%strigObjsid", processNamePrefix_.Data()));
  iEvent.put(trigObjsp4, Form("%strigObjsp4", processNamePrefix_.Data()));
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

  // trigger name
  //const string& trigName = hltConfig_.triggerName(triggerIndex);

  // debug 
  //string fileName = Form("%s.txt", trigName.c_str() );
  //std::ofstream outfile( fileName.c_str(), ios::out | ios::app );

  //
  const trigger::TriggerObjectCollection& triggerObjects = triggerEventH_->getObjects();  // trigger objects
  if (triggerObjects.size() == 0) return;                                                 //
  const vector<string>& moduleLabels = hltConfig_.moduleLabels(triggerIndex);             // modules on this trigger path
  const unsigned int    moduleIndex  = triggerResultsH_->index(triggerIndex);             // index (slot position) of module giving the decision of the path
  unsigned int          nFilters     = triggerEventH_->sizeFilters();                     // number of filters

/////////////////////////////////////
// Show all stored trigger objects //
/////////////////////////////////////

  // CMSSW_4_2x
  // We want to store the filter index corresponding to the last filter for each distinct trigger object type ( the trigger id & p4 can be accessed given the filter index )
  // For Electrons: take the filter index for the last trigger id of either 82 or 92, we don't consider 82 and 92 distinct trigger objects
  // This assumes that each filter will have the same trigger ids... this is *NOT* true for overlap triggers, PFTau triggers for instance
  // If we find that any filter has mixed trigger ids, we store the objects for the last filter only ( this is what was done before and should be right )
  // Yes this is kludged, and yes a universal solution is preferred

  // map ( trigger id => filter index )
  map< int, unsigned int > trigObjsToStore;

  // True if different trigger id's are found in the same filter
  bool mixedTrigIds = false;

  // loop over trigger modules
  unsigned int lastFilterIndex = nFilters;
  for(unsigned int j = 0; j <= moduleIndex; ++j) {

    // get module name & filter index
    const string&      moduleLabel = moduleLabels[j];
    const unsigned int filterIndex = triggerEventH_->filterIndex(InputTag(moduleLabel, "", processName_));

    // debug
    //if( j == 0 ){
    //  outfile << endl << "--------------------------------------------------------------------------------------------------------------" << endl << endl;
    //  outfile << trigName << endl << endl;
    //  outfile << "-> Everything:" << endl << endl;
    //}
      
    // these are the filters with trigger objects filled
    if ( filterIndex < nFilters ){ 
     
      //
      lastFilterIndex = filterIndex;
   
      // get trigger objects & ids
      unsigned int lastEleFilterIndex = 0;
      const trigger::Vids& triggerIds  = triggerEventH_->filterIds(filterIndex);
      const trigger::Keys& triggerKeys = triggerEventH_->filterKeys(filterIndex);
      assert( triggerIds.size() == triggerKeys.size() );

      // True if a filter has hlt objects ( objects with positive trigger id )
      bool hlt = false;   
      
      // Trigger object id for filters associated with a unique trigger object id
      int filterId = -1;

      // loop on trigger objects
      for( unsigned int k = 0; k < triggerKeys.size(); k++){

        // trigger object id
        int id = triggerIds.at(k);

        // trigger p4, p4id
        //const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];   

        // hlt objects
        if( id > 0 ){
          hlt = true;                                 // True if a filter has hlt objects ( objects with positive trigger id )
          if( k == 0 ){                               // Assuming all trigger ids are the same for this filter
            filterId = id;                            // Filter id for all the objects in this filter
            if( filterId == 82 || filterId == 92 ){   // Remember the index of the last filter with all objects having a trigger id of 82 or 92 
              lastEleFilterIndex = filterIndex;
            }
          }
          else {
            if( id != filterId ){
              mixedTrigIds = true;                    // True if different trigger id's are found in the same filter
              //cout << endl << endl << trigName << ", " << moduleLabel << ":\t ERROR in HLTMaker.cc: different trigger object ids found in the same filter... Exiting." << endl << endl;
              //exit(1);
            }
          }
        } // end if( id > 0 )

        // debug
        //if( k == 0 ) outfile << "\t" << filterIndex << ": " << moduleLabel << endl << endl;
        //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );

      } // end loop on trigger objects

      //
      if( hlt ){                                        // only store hlt objects
        assert( filterId != -1 );                       // sanity
        trigObjsToStore[filterId] = filterIndex;        // Store the filter Index ( used to get trigger objects ) for each different trigger type ( filterId )
        if( filterId == 82 ) trigObjsToStore.erase(92); // If this is an electron trigger ( filterId 82 or 92 ) we only want the last one ( that is either 82 or 92 )
        if( filterId == 92 ) trigObjsToStore.erase(82); // 
      }

    } // end if(filterIndex < nFilters)
  }   // end loop over trigger modules
 
//////////////////////////////////////////////
// Show the trigger objects we want to save //
//////////////////////////////////////////////

  // debug
  //outfile << endl << endl << "-> What we want to store";
  //if(mixedTrigIds) outfile << " ( mixed ids... storing last filter only )";
  //outfile << ":" << endl << endl;

  // If we find different trigger ids under any filter we only store the objects for the last filter
  if( mixedTrigIds ){

    // get trigger objects & ids
    const trigger::Vids& triggerIds  = triggerEventH_->filterIds ( lastFilterIndex );
    const trigger::Keys& triggerKeys = triggerEventH_->filterKeys( lastFilterIndex );
    assert( triggerIds.size() == triggerKeys.size() );

    // loop on trigger objects
    for( unsigned int k = 0; k < triggerKeys.size(); k++){

      // trigger object id
      int id = triggerIds.at(k);                                                        

      // trigger p4, p4id
      const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];
  
      // store trigger id, trigger p4, & trigger object id
      p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
      idV.push_back( id );

      // debug
      //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );
    }
 
  }
  // If each filter has trigger objects with the same trigger id, we store the trigger 
  // objects for the last filter of each distinct trigger object type
  else {

    map<int,unsigned int>::iterator it;
    for( it = trigObjsToStore.begin() ; it != trigObjsToStore.end(); it++ ){
    
      // get trigger objects & ids
      const unsigned int filterIndex   = (*it).second;
      const trigger::Vids& triggerIds  = triggerEventH_->filterIds ( filterIndex );
      const trigger::Keys& triggerKeys = triggerEventH_->filterKeys( filterIndex );
      assert( triggerIds.size() == triggerKeys.size() );
    
      // loop on trigger objects
      for( unsigned int k = 0; k < triggerKeys.size(); k++){

        // trigger object id
        int id = triggerIds.at(k);                                                        

        // trigger p4, p4id
        const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys[k] ];
    
        // store trigger id, trigger p4, & trigger object id
        p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
        idV.push_back( id );

        // debug
        //PrintTriggerObjectInfo( outfile, id, triggerObject.particle().p4() );
      } 

    }
  }

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
