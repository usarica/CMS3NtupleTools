#ifndef NTUPLEMAKER_MUONMAKER_H
#define NTUPLEMAKER_MUONMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


typedef math::XYZTLorentzVectorF LorentzVector;


class MuonMaker : public edm::stream::EDProducer<>{
public:
  explicit MuonMaker(const edm::ParameterSet&);

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;
  int year_;

  bool refurbishSelections_;

  edm::EDGetTokenT< edm::View<pat::Muon> > muonsToken;

  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

  edm::EDGetTokenT<reco::VertexCollection> vtxToken;

  edm::EDGetTokenT< double > rhoToken;

  void setCutBasedH4lIdSelectionBits(edm::View<pat::Muon>::const_iterator const&, pat::Muon&) const;

};


#endif
