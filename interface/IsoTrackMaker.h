#ifndef NTUPLEMAKER_ISOTRACKMAKER_H
#define NTUPLEMAKER_ISOTRACKMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

class IsoTrackMaker : public edm::stream::EDProducer<> {
public:
    explicit IsoTrackMaker (const edm::ParameterSet&);
    ~IsoTrackMaker();

private:
    //  virtual void beginJob() ;
    virtual void beginJob() ;
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken;
    edm::EDGetTokenT<pat::IsolatedTrackCollection> isoTracksToken;
    double pt_cut_;
    double pt_cut_noIso_;

    const pat::PackedCandidateCollection *pfCandidates;

};

#endif
