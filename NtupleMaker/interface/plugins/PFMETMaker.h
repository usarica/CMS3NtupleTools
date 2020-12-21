#ifndef NTUPLEMAKER_PFMETMAKER_H
#define NTUPLEMAKER_PFMETMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"


class PFMETMaker : public edm::stream::EDProducer<>{
public:
  explicit PFMETMaker(const edm::ParameterSet&);
  ~PFMETMaker();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;

  edm::EDGetTokenT< edm::View<pat::MET> > metToken;

  bool applyMETfix;
  edm::EDGetTokenT< edm::View<pat::MET> > metDefaultToken;

};


#endif
