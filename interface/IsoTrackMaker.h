#ifndef NTUPLEMAKER_ISOTRACKMAKER_H
#define NTUPLEMAKER_ISOTRACKMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

class IsoTrackMaker : public edm::EDProducer {
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
  edm::InputTag pfCandidatesTag_;
  double isotrack_dz_cut_;
  double isolation_dz_cut_;
  double pflep_pt_cut_;
  double pfhad_pt_cut_;
  double coneR_;

  const pat::PackedCandidateCollection *pfCandidates;

};

#endif

