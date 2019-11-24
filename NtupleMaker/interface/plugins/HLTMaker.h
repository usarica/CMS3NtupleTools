#ifndef NTUPLEMAKER_HLTMAKER_H
#define NTUPLEMAKER_HLTMAKER_H

#include <algorithm>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "CommonTools/TriggerUtils/interface/PrescaleWeightProvider.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//NOT IN miniAOD #include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CMS3/NtupleMaker/interface/TriggerInfo.h"

#include "TRegexp.h"
#include "TString.h"
#include "TBits.h"


typedef math::XYZTLorentzVectorF LorentzVector;


class HLTMaker : public edm::stream::EDProducer<>{
public:
  explicit HLTMaker(const edm::ParameterSet&);
  ~HLTMaker(){}

private:
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;

  std::string processName_;
  TString processNamePrefix_;

  std::vector<std::string> prunedTriggerNames_;

  HLTPrescaleProvider hltConfig_;
  bool doFillInformation;

  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescaleToken;

  std::vector<TriggerInfo> cached_triggerinfos;

  bool doPruneTriggerName(const std::string&) const;

};


#endif
