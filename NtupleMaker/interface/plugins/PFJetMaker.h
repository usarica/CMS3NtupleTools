#ifndef NTUPLEMAKER_PFJETMAKER_H
#define NTUPLEMAKER_PFJETMAKER_H

#include <string>
#include <vector>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"


class PFJetMaker : public edm::stream::EDProducer<>{
public:
  explicit PFJetMaker(const edm::ParameterSet&);
  ~PFJetMaker();

private:
  enum METShiftType{
    kMETShift_JECNominal,
    kMETShift_JECDn,
    kMETShift_JECUp,
    kMETShift_JECNominal_JERNominal,
    kMETShift_JECNominal_JERDn,
    kMETShift_JECNominal_JERUp,
    kMETShift_JECDn_JERNominal,
    kMETShift_JECUp_JERNominal,

    nMETShiftTypes
  };

  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

protected:
  bool printWarnings;

  const std::string aliasprefix_;
  const std::string jetCollection_; // ==JECpayload in main_pset.py

  const bool isMC;
  bool isFatJet;
  bool isPuppi;
  bool METshift_fixEE2017;
  std::vector<std::string> JEClevels;

  unsigned long long cacheId_rcdJEC;

  std::shared_ptr<FactorizedJetCorrector> jetCorrector;
  std::shared_ptr<JetCorrectionUncertainty> jetUncEstimator;

  edm::EDGetTokenT<double> rhoToken;
  edm::EDGetTokenT< reco::VertexCollection > vtxToken;

  edm::EDGetTokenT< edm::View<pat::Jet> > pfJetsToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

  edm::EDGetTokenT< edm::View<reco::GenJet> > genJetsToken;

  void get_reco_gen_matchMap(
    edm::Event const&,
    edm::Handle< edm::View<pat::Jet> > const&, edm::Handle< edm::View<reco::GenJet> > const&,
    std::unordered_map<pat::Jet const*, reco::GenJet const*>&
  ) const;

  void run_JetCorrector_JEC_L123_L1(
    double const& jet_pt_uncorrected, double const& jet_eta, double const& jet_phi,
    double const& jet_area, double const& rho, int const& npv,
    double& JEC_L123, double& JEC_L1
  );
  void run_JetUncertainty(
    double const& jet_pt_corrected, double const& jet_eta, double const& jet_phi,
    double& relJECUnc
  );

  void compute_METShift(
    reco::Particle::LorentzVector const& p4_jet_uncorrected, reco::Particle::LorentzVector const& p4_mucands,
    double const& JEC_L1L2L3, double const& JEC_L1, double const& JERval,
    char const& iJECshift,
    bool& flag_isGoodMET, reco::Particle::LorentzVector& p4_metShift
  );

};


#endif
