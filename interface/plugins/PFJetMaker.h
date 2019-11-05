#ifndef NTUPLEMAKER_PFJETMAKER_H
#define NTUPLEMAKER_PFJETMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


class PFJetMaker : public edm::stream::EDProducer<>{
public:
  explicit PFJetMaker(const edm::ParameterSet&);
  ~PFJetMaker();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;

  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetsToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

};


#endif
