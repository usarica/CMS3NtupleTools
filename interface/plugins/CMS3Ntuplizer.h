#ifndef NTUPLEMAKER_CMS3NTUPLIZER_H
#define NTUPLEMAKER_CMS3NTUPLIZER_H

#include <cassert>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <memory>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/Run.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <FWCore/ServiceRegistry/interface/Service.h>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <DataFormats/JetReco/interface/PFJet.h>

#include <SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include <CMSDataTools/AnalysisTree/interface/SimpleEntry.h>
#include <CMSDataTools/AnalysisTree/interface/BaseTree.h>

#include <CMS3/NtupleMaker/interface/GenInfo.h>
#include <CMS3/NtupleMaker/interface/TriggerInfo.h>
#include <CMS3/NtupleMaker/interface/METFilterInfo.h>
#include <CMS3/NtupleMaker/interface/METInfo.h>


class CMS3Ntuplizer : public edm::EDAnalyzer{
public:
  explicit CMS3Ntuplizer(const edm::ParameterSet&);
  ~CMS3Ntuplizer();

protected:
  const edm::ParameterSet pset;
  BaseTree* outtree;
  SimpleEntry commonEntry;
  bool firstEvent;

  int year;
  TString treename;
  //TString outfilename;
  bool isMC;

  std::string prefiringWeightsTag;
  bool applyPrefiringWeights;

  edm::EDGetTokenT< edm::View<pat::Electron> > electronsToken;
  edm::EDGetTokenT< edm::View<pat::Photon> > photonsToken;
  edm::EDGetTokenT< edm::View<pat::Muon> > muonsToken;
  edm::EDGetTokenT< edm::View<pat::Jet> > ak4jetsToken;
  edm::EDGetTokenT< edm::View<pat::Jet> > ak8jetsToken;

  edm::EDGetTokenT< METInfo > pfmetToken;
  edm::EDGetTokenT< METInfo > puppimetToken;

  edm::EDGetTokenT< reco::VertexCollection > vtxToken;

  edm::EDGetTokenT< double > rhoToken;
  edm::EDGetTokenT< edm::View<TriggerInfo> > triggerInfoToken;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > puInfoToken;
  edm::EDGetTokenT< METFilterInfo > metFilterInfoToken;

  edm::EDGetTokenT< double > prefiringWeightToken;
  edm::EDGetTokenT< double > prefiringWeightToken_Dn;
  edm::EDGetTokenT< double > prefiringWeightToken_Up;

  edm::EDGetTokenT< GenInfo > genInfoToken;
  edm::EDGetTokenT< reco::GenParticleCollection > prunedGenParticlesToken;
  edm::EDGetTokenT< pat::PackedGenParticleCollection > packedGenParticlesToken;
  edm::EDGetTokenT< edm::View<reco::GenJet> > genAK4JetsToken;
  edm::EDGetTokenT< edm::View<reco::GenJet> > genAK8JetsToken;


  void recordGenInfo(const edm::Event&);
  void recordGenParticles(const edm::Event&, std::vector<reco::GenParticle const*>*, std::vector<pat::PackedGenParticle const*>*);
  void recordGenJets(const edm::Event&, bool const&, std::vector<reco::GenJet const*>*);

  size_t fillElectrons(const edm::Event&, std::vector<pat::Electron const*>*);
  size_t fillPhotons(const edm::Event&, std::vector<pat::Photon const*>*);
  size_t fillMuons(const edm::Event&, std::vector<pat::Muon const*>*);
  size_t fillAK4Jets(const edm::Event&, std::vector<pat::Jet const*>*);
  size_t fillAK8Jets(const edm::Event&, std::vector<pat::Jet const*>*);
  size_t fillVertices(const edm::Event&, std::vector<reco::Vertex const*>*);

  bool fillEventVariables(const edm::Event&);
  bool fillTriggerInfo(const edm::Event&);
  bool fillMETFilterVariables(const edm::Event&);
  bool fillMETVariables(const edm::Event&);

  bool fillGenVariables(
    const edm::Event&,
    std::vector<reco::GenParticle const*>*,
    std::vector<pat::PackedGenParticle const*>*,
    std::vector<reco::GenJet const*>*,
    std::vector<reco::GenJet const*>*
  );

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);

  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

};


#endif
