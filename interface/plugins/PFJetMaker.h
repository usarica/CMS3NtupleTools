#ifndef NTUPLEMAKER_PFJETMAKER_H
#define NTUPLEMAKER_PFJETMAKER_H

#include <string>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"


class PFJetMaker : public edm::stream::EDProducer<>{
public:
  explicit PFJetMaker(const edm::ParameterSet&);
  ~PFJetMaker();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  const std::string aliasprefix_;
  const std::string jetCollection_;
  const bool isMC;
  bool isFatJet;
  bool isPuppi;

  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT< reco::VertexCollection > vtxToken;

  edm::EDGetTokenT< edm::View<pat::Jet> > pfJetsToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

  edm::EDGetTokenT< edm::View<reco::GenJet> > genJetsToken;


  std::unordered_map<pat::Jet const*, reco::GenJet const*> get_reco_gen_matchMap(edm::Event const&, edm::Handle< edm::View<pat::Jet> > const&) const;

};


#endif
