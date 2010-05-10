#include "CMS2/NtupleMaker/interface/HLTMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace edm;
using namespace reco;
using namespace std;

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
  produces<vector<unsigned int> >           (Form("%sprescales"  ,processNamePrefix_.Data())).setBranchAlias(Form("%s_prescales"   ,processNamePrefix_.Data()));
  produces<vector<vector<int> > >           (Form("%strigObjsid",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_id",processNamePrefix_.Data()));

  produces<vector<vector<LorentzVector> > > (Form("%strigObjsp4",processNamePrefix_.Data())).setBranchAlias(Form("%s_trigObjs_p4",processNamePrefix_.Data()));
    
}

void HLTMaker::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{
  // HLT config _should no longer_ change within runs :)
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
  } else 
    throw cms::Exception("HLTMaker::beginRun: config extraction failure with process name " + processName_);
}

void HLTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iEvent.getByLabel(edm::InputTag("TriggerResults",       "", processName_), triggerResultsH_);
  if (! triggerResultsH_.isValid())
    throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
  if (fillTriggerObjects_) {
    iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", processName_), triggerEventH_  );
    if (! triggerEventH_.isValid()  )
      throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
  }
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
  if (nTriggers > 256)
    throw cms::Exception("HLTMaker::produce: number of HLT trigger variables must be increased!");

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
      prescales->push_back(hltConfig_.prescaleValue(0, name));
	
	
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
    TRegexp reg(Form("%s", prunedTriggerNames_[i].c_str()), true);
    TString sname(name);
    if (sname.Index(reg) >= 0)
      return true;
  }
  return false;
}

void HLTMaker::fillTriggerObjectInfo(unsigned int triggerIndex,
				     vector<int>& idV,
				     vector<LorentzVector>& p4V) const
{
  const trigger::TriggerObjectCollection& triggerObjects = triggerEventH_->getObjects();
  if (triggerObjects.size() == 0) return;

  // modules on this trigger path
  const vector<string>& moduleLabels = hltConfig_.moduleLabels(triggerIndex);
  // index (slot position) of module giving the decision of the path
  const unsigned int moduleIndex = triggerResultsH_->index(triggerIndex);

  unsigned int nFilters = triggerEventH_->sizeFilters();
  // the first and last filter information is stored
  // but we just want the last filter
  unsigned int lastFilterIndex = nFilters;
  for(unsigned int j = 0; j <= moduleIndex; ++j) {
    const string& moduleLabel = moduleLabels[j];
    const unsigned int filterIndex = triggerEventH_->filterIndex(InputTag(moduleLabel, "", processName_));
    if (filterIndex < nFilters)
      lastFilterIndex = filterIndex;
  }
  if (lastFilterIndex < nFilters) {
    const trigger::Vids& triggerIds = triggerEventH_->filterIds(lastFilterIndex);
    const trigger::Keys& triggerKeys = triggerEventH_->filterKeys(lastFilterIndex);
    assert(triggerIds.size()==triggerKeys.size());

    for(unsigned int j = 0; j < triggerKeys.size(); ++j) {
      const trigger::TriggerObject& triggerObject = triggerObjects[triggerKeys[j]];
      p4V.push_back( LorentzVector( triggerObject.particle().p4() ) );
      idV.push_back(triggerObject.id());
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTMaker);
