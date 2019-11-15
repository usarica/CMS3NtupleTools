#ifndef NTUPLEMAKER_SUBJETMAKER_H
#define NTUPLEMAKER_SUBJETMAKER_H

#include <string>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "NNKit/FatJetNN/interface/FatJetNN.h"


class SubJetMaker : public edm::stream::EDProducer<>{
public:
  explicit SubJetMaker(const edm::ParameterSet&);
  ~SubJetMaker();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

protected:
  const std::string jetCollection_;
  const bool isMC;

  deepntuples::FatJetNN* fatjetNN_;

  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetsToken;
  edm::InputTag pfCandidatesTag_;
  double pfJetPtCut_;
  bool keepless_;
  std::string aliasprefix_;
  std::string PFJetCorrectorL2L3_;
  std::string PFJetCorrectorL1L2L3_;
  std::string PFJetCorrectorL1FastL2L3_;
  std::string PFJetCorrectorL1Fast_;
  std::string PFJetCorrectorL1FastL2L3residual_;
};


#endif
