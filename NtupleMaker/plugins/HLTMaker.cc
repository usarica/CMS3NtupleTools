#include "CMS3/NtupleMaker/interface/plugins/HLTMaker.h"


typedef math::XYZTLorentzVectorF LorentzVector;

using namespace edm;
using namespace reco;
using namespace std;


HLTMaker::HLTMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),

  processName_(iConfig.getUntrackedParameter<string>("processName")),
  processNamePrefix_(aliasprefix_),

  prunedTriggerNames_(iConfig.getUntrackedParameter< std::vector<std::string> >("prunedTriggerNames")),

  hltConfig_(iConfig, consumesCollector(), *this),
  doFillInformation(true)
{
  triggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", processName_));
  triggerPrescaleToken = consumes<pat::PackedTriggerPrescales>(iConfig.getUntrackedParameter<std::string>("triggerPrescaleInputTag"));

  produces< std::vector<TriggerInfo> >().setBranchAlias(aliasprefix_);
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
}

void HLTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique< std::vector<TriggerInfo> >();

  // If the process name is not specified retrieve the latest
  // TriggerEvent object and the corresponding TriggerResults.
  // We should only have to do this once though, the next time
  // produce is called processName_ should be set.

  // Now using a single processName_ (set to "HLT" in the configuration file). Is this OK? Do we need the flexibility we had before?
  edm::Handle<edm::TriggerResults> triggerResultsH_;
  iEvent.getByToken(triggerResultsToken, triggerResultsH_);
  if (!triggerResultsH_.isValid()) throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
  size_t nTriggers = triggerResultsH_->size();
  //if (nTriggers > 768) throw cms::Exception( Form("HLTMaker::produce: number of HLT trigger variables must be increased! ( %d > 768 )", nTriggers) );
  result->reserve(nTriggers);

  // if it's data, cache for only a single lumi block, otherwise cache for whole job
  bool isdata = iEvent.isRealData();
  bool make_cache = doFillInformation;
  if (make_cache){
    edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResultsH_); // Does this have to be done for every event?

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH_;
    iEvent.getByToken(triggerPrescaleToken, triggerPrescalesH_);
    if (!triggerPrescalesH_.isValid()) throw cms::Exception("HLTMaker::produce: error getting PackedTriggerPrescales product from Event!");

    //// Printing miniAOD content...
    //for (unsigned int i=0; i<triggerNames_.size(); i++) std::cout << "triggerNames= " << triggerNames_.triggerName(i) << std::endl;
    //if ( triggerObjectStandAlonesH_.isValid()) cout<<"Got triggerObjectStandAlonesHandle with size "<<triggerObjectStandAlonesH_->size()<<endl;
    //else cout<<"Couldn't find triggerObjectStandAlonesH"<<endl;
    //for ( uint i = 0; i < triggerObjectStandAlonesH_->size(); i++ ) {
    //  pat::TriggerObjectStandAlone TO = triggerObjectStandAlonesH_->at(i);
    //  TO.unpackPathNames( triggerNames_ );
    //  cout<<"Trigger Object "<<i<<"has pt, eta, phi, id = "<< TO.pt() <<" "<<TO.eta()<<" "<<TO.phi()<<" "<<TO.pdgId()<<endl;
    //  cout<<"Trigger Object has hasPathLastFilterAccepted() "<<TO.hasPathLastFilterAccepted()<<endl;
    //  cout<<"Trigger Object "<<i<<" has collection() "<<TO.collection()<<endl;
    //  std::vector< std::string > path_names = TO.pathNames(false); //TO associated to path
    //  cout<<"Trigger Object "<<i<<" associated to "<< path_names.size() <<" pathNames(false): ";
    //  for (uint j  = 0; j < path_names.size(); j++) cout<<path_names[j]<<" ";
    //  cout<<endl;
    //  std::vector< std::string > path_namesPASS = TO.pathNames(true); //make sure they passed!
    //  cout<<"Trigger Object "<<i<<" passed "<< path_namesPASS.size() <<" pathNames(true): ";
    //  for (uint j  = 0; j < path_namesPASS.size(); j++) cout<<path_namesPASS[j]<<" ";
    //  cout<<endl;
    //
    //  std::vector< std::string > filter_labels = TO.filterLabels();
    //  cout<<"Trigger Object "<<i<<" has "<< filter_labels.size() <<" filter_labels: ";
    //  for (uint j  = 0; j < filter_labels.size(); j++) cout<<filter_labels[j]<<" ";
    //  cout<<endl;
    //
    //  std::vector< int > filter_ids = TO.filterIds();
    //  cout<<"Trigger Object "<<i<<" has "<< filter_ids.size() <<" filter_ids: ";
    //  for (uint j  = 0; j < filter_ids.size(); j++) cout<<filter_ids[j]<<" ";
    //  cout<<endl;
    //}
    //// END printing miniAOD content

    // flip this flag for subsequent events
    doFillInformation = false;
    cached_triggerinfos.reserve(nTriggers);

    for (unsigned int i = 0; i < nTriggers; ++i){
      // What is your name?
      const std::string& name = triggerNames_.triggerName(i);
      bool passTrigger = triggerResultsH_->accept(i);
      unsigned int HLTprescale = 1;
      unsigned int L1prescale = 1;

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
          // Now, IFF (the l1prescale isn't 1 and there's more than 1 nonzero L1 feeding into a HLT path),
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
            float s = 0.;
            for (auto it = l1prescalevals.begin(); it != new_end; it++) s += 1./(*it);
            harmonicmean = (int) (1./s);
          }
          if ((harmonicmean != 1) && (harmonicmean != l1prescale)) l1prescale = harmonicmean;
          L1prescale = l1prescale;
        }
      }
      else{
        HLTprescale = hltConfig_.prescaleValue(iEvent, iSetup, name);
        // L1prescale = 1 is already the case
      }

      cached_triggerinfos.emplace_back(name, passTrigger, HLTprescale, L1prescale);
    }
  }

  assert(cached_triggerinfos.size() == nTriggers);
  for (auto const& obj:cached_triggerinfos) result->emplace_back(obj);
  //std::copy(cached_triggerinfos.begin(), cached_triggerinfos.end(), result->begin());

  iEvent.put(std::move(result));
}

bool HLTMaker::doPruneTriggerName(const string& name) const{
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


DEFINE_FWK_MODULE(HLTMaker);
