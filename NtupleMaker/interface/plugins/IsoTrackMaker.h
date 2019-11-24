#ifndef NTUPLEMAKER_ISOTRACKMAKER_H
#define NTUPLEMAKER_ISOTRACKMAKER_H

#include <string>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

#include <CMS3/NtupleMaker/interface/IsoTrackInfo.h>


class IsoTrackMaker : public edm::stream::EDProducer<>{
public:
  explicit IsoTrackMaker(const edm::ParameterSet&);
  ~IsoTrackMaker();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;

  edm::EDGetTokenT<pat::IsolatedTrackCollection> isoTracksToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

};


#endif
