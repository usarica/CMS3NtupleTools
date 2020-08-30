#include <limits>

#include <FWCore/Utilities/interface/EDMException.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include <CMS3/NtupleMaker/interface/plugins/HLTMaker.h>
#include <CMS3/NtupleMaker/interface/TriggerObjectInfo.h>

#include "MELAStreamHelpers.hh"


typedef math::XYZTLorentzVectorF LorentzVector;

using namespace edm;
using namespace reco;
using namespace std;
using namespace MELAStreamHelpers;


HLTMaker::HLTMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),
  processName_(iConfig.getUntrackedParameter<string>("processName")),

  prunedTriggerNames_(iConfig.getUntrackedParameter< std::vector<std::string> >("prunedTriggerNames")),
  recordFilteredTrigObjects_(iConfig.getParameter<bool>("recordFilteredTrigObjects")),

  hltConfig_(iConfig, consumesCollector(), *this),
  doFillInformation(true)
{
  triggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", processName_));
  triggerPrescaleToken = consumes<pat::PackedTriggerPrescales>(iConfig.getUntrackedParameter<std::string>("triggerPrescaleInputTag"));
  if (recordFilteredTrigObjects_){
    triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getUntrackedParameter<std::string>("triggerObjectsName"));
  }

  produces< std::vector<TriggerInfo> >().setBranchAlias(aliasprefix_);
  if (recordFilteredTrigObjects_){
    produces< std::vector<TriggerObjectInfo> >("filteredTriggerObjectInfos");
  }
}

void HLTMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){
  // In the case that we are choosing the process name
  // automatically, i.e. the processName_ parameter is
  // an empty string, we can't init  HLTConfigProvider
  // until after we've determined the process name. So
  // don't init here until after we've set processName_
  // in the produce method and init there once and only
  // once. Sounds scary, it is kinda!
  // HLT config _should no longer_ change within runs :)
  if (processName_ != ""){
    bool changed(true);
    // if (hltConfig_.init(iRun,iSetup,"*",changed)) {
    if (!hltConfig_.init(iRun, iSetup, processName_, changed)) throw cms::Exception("HLTMaker::beginRun: config extraction failure with process name " + processName_);
  }
}

void HLTMaker::beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&){
  doFillInformation = true;
  cached_triggerinfos.clear();
  cached_allTriggerNames.clear();
}

void HLTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique< std::vector<TriggerInfo> >();
  std::unique_ptr< std::vector<TriggerObjectInfo> > filteredTriggerObjectInfos = (recordFilteredTrigObjects_ ? std::make_unique< std::vector<TriggerObjectInfo> >() : nullptr);

  // if it is data, cache for only a single lumi block, otherwise cache for whole job
  bool isdata = iEvent.isRealData();
  bool make_cache = doFillInformation;

  // If the process name is not specified retrieve the latest
  // TriggerEvent object and the corresponding TriggerResults.
  // We should only have to do this once though, the next time
  // produce is called processName_ should be set.

  // Now using a single processName_ (set to "HLT" in the configuration file). Is this OK? Do we need the flexibility we had before?
  edm::Handle<edm::TriggerResults> triggerResultsH_;
  iEvent.getByToken(triggerResultsToken, triggerResultsH_);
  if (!triggerResultsH_.isValid()) throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
  edm::TriggerNames const& triggerNames_ = iEvent.triggerNames(*triggerResultsH_);
  size_t nTriggers = triggerResultsH_->size();
  result->reserve(nTriggers);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjectStandAlonesH_;
  if (recordFilteredTrigObjects_){
    iEvent.getByToken(triggerObjectsToken, triggerObjectStandAlonesH_);
    if (!triggerObjectStandAlonesH_.isValid()) throw cms::Exception("HLTMaker::produce: error getting TriggerObjectsStandAlone product from Event!");
  }

  if (make_cache){
    cached_triggerNamesPSetId = triggerNames_.parameterSetID();
    cached_allTriggerNames.reserve(nTriggers);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH_;
    iEvent.getByToken(triggerPrescaleToken, triggerPrescalesH_);
    if (!triggerPrescalesH_.isValid()) throw cms::Exception("HLTMaker::produce: error getting PackedTriggerPrescales product from Event!");

    // flip this flag for subsequent events
    doFillInformation = false;
    cached_triggerinfos.reserve(nTriggers);

    for (size_t i = 0; i < nTriggers; ++i){
      // What is your name?
      const std::string& name = triggerNames_.triggerName(i);
      cached_allTriggerNames.emplace_back(name);
      //std::string name = getTrimmedTriggerName(triggerNames_.triggerName(i));
      if (!pruneTriggerByName(name)) continue;

      bool passTrigger = triggerResultsH_->accept(i);
      int HLTprescale = 1;
      int L1prescale = 1;

      // What is your prescale?
      // Buggy way in miniAOD:
      // prescales->push_back( triggerPrescalesH_.isValid() ? triggerPrescalesH_->getPrescaleForIndex(i) : -1 );
      if (isdata){
        // get prescale info from hltConfig_
        std::pair< std::vector< std::pair<std::string, int> >, int> detailedPrescaleInfo = hltConfig_.prescaleValuesInDetail(iEvent, iSetup, name);

        // prescale
        HLTprescale = (triggerPrescalesH_.isValid() ? detailedPrescaleInfo.second : -1);

        // save l1 prescale values in standalone vector
        std::vector<int> l1prescalevals;
        for (auto const& vv:detailedPrescaleInfo.first) l1prescalevals.push_back(vv.second);

        // find and save minimum l1 prescale of any ORed L1 that seeds the HLT
        bool isAllZeros = std::all_of(l1prescalevals.begin(), l1prescalevals.end(), [] (int i) { return i==0; });
        if (isAllZeros) L1prescale = 0;
        else{
          // first remove all values that are 0
          std::vector<int>::iterator new_end = std::remove(l1prescalevals.begin(), l1prescalevals.end(), 0);

          // now find the minimum
          std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals), new_end);
          int minind = std::distance(l1prescalevals.begin(), result);
          int l1prescale = (minind < std::distance(l1prescalevals.begin(), new_end) ? l1prescalevals.at(minind) : -1);
          // Update 5 March 2020: Instead of harmonic mean, calculate prescale from 1/(1-f), where f is the probability to fail all L1 seeds.
          // Until 5 March 2020: Now, IFF (the l1prescale isn't 1 and there's more than 1 nonzero L1 feeding into a HLT path),
          // compute the harmonic mean (or whatever it's called -- think "resistors in parallel") of the L1 prescales
          // and use that as the `l1prescale`.
          // E.g., In 2017, HLT_Mu8_v8 has Mu pt 3, 5, 7 L1 seeds with prescales 16000, 3600, 1500.
          // Taking the minimum gives a prescale of 1500, but the harmonic mean gives 993, so that we are
          // overestimating the rate in Data by ~50%. This showed up when comparing Z peak integrals with MC.
          // Checking on some 2017 data, the condition to use the harmonic mean instead of the minimum
          // only happens for some prescaled single muon and electron triggers used for RA5 fake rate derivation
          // and results in 15-40% differences in Z peak ratios.
          int harmonicmean = 1;
          if (l1prescale != 1 && (new_end-std::begin(l1prescalevals) > 1)){
            double s = 1.;
            for (auto it = l1prescalevals.begin(); it != new_end; it++) s *= (1.-1./double(*it));
            harmonicmean = (int) (1./(1.-s) + 0.5); // +0.5 is to make sure rounding is done to the nearest integer, not just rounding down
          }
          if ((harmonicmean != 1) && (harmonicmean != l1prescale)) l1prescale = harmonicmean;
          L1prescale = l1prescale;
        }
      }
      else{
        HLTprescale = hltConfig_.prescaleValue(iEvent, iSetup, name);
        // L1prescale = 1 is already the case
      }

      // Must pass 'i', the absolute index
      cached_triggerinfos.emplace_back(name, i, passTrigger, HLTprescale, L1prescale);
    }
  }

  for (auto& obj:cached_triggerinfos){
    result->emplace_back(obj);
    // Update trigger accept flag
    auto& tmp_obj = result->back();
    tmp_obj.passTrigger = triggerResultsH_->accept(tmp_obj.index); // Must refer to absolute trigger collection index
    /*
    if (tmp_obj.name!=triggerNames_.triggerName(tmp_obj.index)) throw cms::Exception("InvalidTriggerName")
      << "Cached trigger info name " << tmp_obj.name << " != " << triggerNames_.triggerName(tmp_obj.index)
      << endl;
    */
  }
  //std::copy(cached_triggerinfos.begin(), cached_triggerinfos.end(), result->begin());

  if (recordFilteredTrigObjects_){
    if (cached_triggerinfos.size()>std::numeric_limits<cms3_triggerIndex_t>::max()) throw cms::Exception("NumericPrecision")
      << "The size of cached trigger objects (" << cached_triggerinfos.size()
      << ") exceeds the limit of cms3_triggerIndex_t (" << std::numeric_limits<cms3_triggerIndex_t>::max() << ")";

    size_t nTOs = triggerObjectStandAlonesH_->size();
    filteredTriggerObjectInfos->reserve(nTOs);
    size_t iTO=0;
    for (auto const& triggerObjectStandAlone:(*triggerObjectStandAlonesH_)){
      std::vector<int> filter_ids = triggerObjectStandAlone.filterIds();
      bool isL1Object = false;
      for (auto const& id:filter_ids){ if (id<0){ isL1Object=true; break; } }
      if (isL1Object) continue; // Do not look at L1 trigger objects

      pat::TriggerObjectStandAlone TO = triggerObjectStandAlone; // Copy because we will have to unpack.
      TO.unpackPathNames(triggerNames_);

      std::vector< std::string > path_namesASSOCIATED = TO.pathNames(false); // 'false' refers to making sure that the object associated with the filters for a given trigger
      std::vector< std::string > path_namesPASS = TO.pathNames(true); // 'true' refers to making sure that the object passed all filters for a given trigger
      std::vector<std::pair<cms3_triggerIndex_t, bool>> associatedTriggerIndices; associatedTriggerIndices.reserve(cached_triggerinfos.size());
      {
        unsigned int iCachedTriggerInfo = 0;
        for (auto const& cached_triggerinfo:cached_triggerinfos){
          // Check if associated, then if passed
          if (std::find(path_namesASSOCIATED.begin(), path_namesASSOCIATED.end(), cached_triggerinfo.name)!=path_namesASSOCIATED.end()){
            bool passAllFilters = (std::find(path_namesPASS.begin(), path_namesPASS.end(), cached_triggerinfo.name)!=path_namesPASS.end());
#if TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL == 0
            associatedTriggerIndices.emplace_back(iCachedTriggerInfo, passAllFilters);
#else
            associatedTriggerIndices.emplace_back(cached_triggerinfo.index, passAllFilters);
#endif
          }
          iCachedTriggerInfo++;
        }
      }
      if (!associatedTriggerIndices.empty()){
        // Need to filter the filter id ints for 0s
        unsigned int nFilterIds = filter_ids.size();
        std::vector<cms3_triggertype_t> filter_id_types; filter_id_types.reserve(nFilterIds);
        for (auto const& filter_id:filter_ids){ if (nFilterIds==1 || filter_id!=0) filter_id_types.push_back(filter_id); }
        if (filter_id_types.empty()){
          MELAerr << "Filter id types for trigger object " << iTO << " are all 0!" << endl;
          MELAerr << "\t- path_namesASSOCIATED = " << path_namesASSOCIATED << endl;
          MELAerr << "\t- path_namesPASS = " << path_namesPASS << endl;
          MELAerr << "\t- (pt, eta, phi, mass) = ( " << TO.p4().Pt() << ", " << TO.p4().Eta() << ", " << TO.p4().Phi() << ", " << TO.p4().M() << " )" << endl;
          filter_id_types.emplace_back(0);
        }

        filteredTriggerObjectInfos->emplace_back(iTO, filter_id_types, TO.p4());
        auto& trig_obj = filteredTriggerObjectInfos->back();
        for (auto const& ati_pair:associatedTriggerIndices) trig_obj.addTriggerCollectionIndex(ati_pair.first, ati_pair.second);
      }

      iTO++;
    }

    /*
    // Printing miniAOD content...
    if (make_cache){
      cout << "Trigger names pset id = " << cached_triggerNamesPSetId << endl;
      for (unsigned int i=0; i<triggerNames_.size(); i++) std::cout << "triggerName = " << triggerNames_.triggerName(i) << std::endl;
    }
    else{
      for (unsigned int i=0; i<triggerNames_.size(); i++){
        if (cached_allTriggerNames.at(i)!=triggerNames_.triggerName(i)) std::cerr << "triggerName " << triggerNames_.triggerName(i) << " IS NOT THE SAME AS " << cached_allTriggerNames.at(i) << std::endl;
      }
    }
    for (uint i = 0; i < nTOs; i++){
      pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
      TO.unpackPathNames(triggerNames_);
      TO.unpackFilterLabels(iEvent, *triggerResultsH_);
      std::vector< std::string > path_names = TO.pathNames(false); //TO associated to path
      std::vector< std::string > path_namesPASS = TO.pathNames(true); //make sure they passed!
      std::vector< std::string > filter_labels = TO.filterLabels();
      std::vector< int > filter_ids = TO.filterIds();
      bool doPrint = false;
      for (auto const& cached_triggerinfo:cached_triggerinfos){
        if (cached_triggerinfo.name.find("Photon")!=std::string::npos && std::find(path_namesPASS.begin(), path_namesPASS.end(), cached_triggerinfo.name)!=path_namesPASS.end()){
          doPrint = true;
          break;
        }
      }
      if (doPrint){
        cout<<"Trigger Object "<<i<<" / "<<nTOs<<":"<<endl;
        cout<<"\t- pt, eta, phi, id = "<< TO.pt() <<" "<<TO.eta()<<" "<<TO.phi()<<" "<<TO.pdgId()<<endl;
        cout<<"\t- hasPathLastFilterAccepted() "<<TO.hasPathLastFilterAccepted()<<endl;
        cout<<"\t- Has collection() "<<TO.collection()<<endl;
        cout<<"\t- Associated to "<< path_names.size() <<" pathNames(false): ";
        for (uint j = 0; j < path_names.size(); j++) cout<<path_names[j]<<" ";
        cout<<endl;
        cout<<"\t- Passed "<< path_namesPASS.size() <<" pathNames(true): ";
        for (uint j = 0; j < path_namesPASS.size(); j++) cout<<path_namesPASS[j]<<" ";
        cout<<endl;

        cout<<"\t- Has "<< filter_labels.size() <<" filter_labels: ";
        for (uint j = 0; j < filter_labels.size(); j++) cout<<filter_labels[j]<<" ";
        cout<<endl;

        // Filter ids are listed in DataFormats/HLTReco/interface/TriggerTypeDefs.h
        cout<<"\t- Has "<< filter_ids.size() <<" filter_ids: ";
        for (uint j = 0; j < filter_ids.size(); j++) cout<<filter_ids[j]<<" ";
        cout<<endl;
      }
    }
    // END printing miniAOD content
    */
  }

  iEvent.put(std::move(result));
  if (recordFilteredTrigObjects_) iEvent.put(std::move(filteredTriggerObjectInfos), "filteredTriggerObjectInfos");
}

bool HLTMaker::pruneTriggerByName(const string& name) const{
  if (prunedTriggerNames_.empty()) return true;
  for (TString const& ptrigname:prunedTriggerNames_){
    TString pattern = ptrigname;
    pattern.ToLower();
    TRegexp reg(Form("%s", pattern.Data()), true);
    TString sname(name);
    sname.ToLower();
    if (sname.Index(reg) >= 0) return true;
  }
  return false;
}

std::string HLTMaker::getTrimmedTriggerName(std::string const& name){
  std::vector<std::string> splitname;
  HelperFunctions::splitOptionRecursive(name, splitname, 'v', false);
  std::string newname;
  if (splitname.size()==1) newname = name;
  else{
    for (size_t is=0; is<splitname.size()-1; is++) newname += splitname.at(is)+"v";
  }
  //MELAout << "HLTMaker::getTrimmedTriggerName: Trigger name " << name << " -> " << newname << endl;
  return newname;
}


DEFINE_FWK_MODULE(HLTMaker);
