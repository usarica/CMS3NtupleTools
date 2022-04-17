#ifndef NTUPLEMAKER_GENMAKER_H
#define NTUPLEMAKER_GENMAKER_H

#include <memory>
#include <utility>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include <CommonLHETools/LHEHandler/interface/LHEHandler.h>
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include <CMS3/NtupleMaker/interface/GenInfo.h>
#include <CMS3/NtupleMaker/interface/KFactorHelpers.h>


class GenMaker : public edm::one::EDProducer<edm::one::WatchRuns, edm::one::SharedResources>{
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
  std::string recoMode;

  float xsecOverride;
  float brOverride;

  std::vector<std::pair<KFactorHelpers::KFactorType, KFactorHelpers::KFactorType>> kfactor_num_denum_list;

  edm::InputTag LHEInputTag_;
  edm::InputTag genEvtInfoInputTag_;
  edm::InputTag prunedGenParticlesInputTag_;
  edm::InputTag packedGenParticlesInputTag_;
  edm::InputTag genJetsInputTag_;
  edm::InputTag genMETInputTag_;
  bool ntuplePackedGenParticles_;

  float superMH;

  bool doHiggsKinematics;
  MELAEvent::CandidateVVMode candVVmode;
  int decayVVmode;
  std::vector<std::string> lheMElist;
  
  std::shared_ptr<KFactorHelpers::KFactorHandler_QCD_ggVV_Sig> KFactor_QCD_ggVV_Sig_handle;
  std::shared_ptr<KFactorHelpers::KFactorHandler_QCD_qqVV_Bkg> KFactor_QCD_qqVV_Bkg_handle;
  std::shared_ptr<KFactorHelpers::KFactorHandler_EW_qqVV_Bkg> KFactor_EW_qqVV_Bkg_handle;

  edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken;

  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;

  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticlesToken;
  edm::EDGetTokenT< edm::View<reco::GenJet> > genJetsToken;
  edm::EDGetTokenT< edm::View<pat::MET> > genMETToken;

  std::shared_ptr<LHEHandler> lheHandler_default; // LHEHandler for default PDFs
  std::shared_ptr<LHEHandler> lheHandler_NNPDF30_NLO; // LHEHandler for the 2016-like PDFs

  LHEHandler::RunMode getLHEHandlerRunMode() const;

  /******************/
  /* ME COMPUTATION */
  /******************/
  IvyMELAHelpers::GMECBlock lheMEblock;
  void setupMELA();
  void doMELA(MELACandidate*, GenInfo&);
  void cleanMELA();

  /************************/
  /* K FACTOR COMPUTATION */
  /************************/
  void setupKFactorHandles(edm::ParameterSet const&);

};


#endif
