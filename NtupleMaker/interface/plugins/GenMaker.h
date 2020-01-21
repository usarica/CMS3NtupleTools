#ifndef NTUPLEMAKER_GENMAKER_H
#define NTUPLEMAKER_GENMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "CommonLHETools/LHEHandler/interface/LHEHandler.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include <CMS3/NtupleMaker/interface/GenInfo.h>


class GenMaker : public edm::one::EDProducer<edm::one::WatchRuns, edm::one::SharedResources>{
//class GenMaker : public edm::EDProducer{
public:
  explicit GenMaker(const edm::ParameterSet&);
  ~GenMaker();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&){}

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  std::string aliasprefix_;
  int year;

  float xsecOverride;
  float brOverride;

  edm::InputTag LHEInputTag_;
  edm::InputTag genEvtInfoInputTag_;
  edm::InputTag prunedGenParticlesInputTag_;
  edm::InputTag packedGenParticlesInputTag_;
  edm::InputTag genMETInputTag_;
  bool ntuplePackedGenParticles_;

  float superMH;

  bool doHiggsKinematics;
  MELAEvent::CandidateVVMode candVVmode;
  int decayVVmode;
  std::vector<std::string> lheMElist;

  edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken;

  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;

  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticlesToken;
  edm::EDGetTokenT< edm::View<pat::MET> > genMETToken;

  std::shared_ptr<LHEHandler> lheHandler_default; // LHEHandler for default PDFs
  std::shared_ptr<LHEHandler> lheHandler_NNPDF30_NLO; // LHEHandler for the 2016-like PDFs

  /******************/
  /* ME COMPUTATION */
  /******************/
  CMS3MELAHelpers::GMECBlock lheMEblock;
  void setupMELA();
  void doMELA(MELACandidate*, GenInfo&);
  void cleanMELA();

};


#endif
