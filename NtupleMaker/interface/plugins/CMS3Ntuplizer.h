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
#include <FWCore/Framework/interface/one/EDAnalyzer.h>
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
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
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
#include <CMS3/NtupleMaker/interface/TriggerObjectInfo.h>
#include <CMS3/NtupleMaker/interface/METFilterInfo.h>
#include <CMS3/NtupleMaker/interface/METInfo.h>
#include <CMS3/NtupleMaker/interface/IsotrackInfo.h>
#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>


class CMS3Ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>{
public:
  explicit CMS3Ntuplizer(const edm::ParameterSet&);
  ~CMS3Ntuplizer();

protected:
  enum ParticleRecordLevel{
    kNone=0,
    kPromptFinalStatePhotons,
    kReducedFinalStates,
    kAllFinalStates,
    kReducedFinalStatesAndHardProcesses,
    kAll,
    nParticleRecordLevels
  };

  static const std::string colName_muons;
  static const std::string colName_electrons;
  static const std::string colName_photons;
  static const std::string colName_fsrcands;
  static const std::string colName_isotracks;
  static const std::string colName_ak4jets;
  static const std::string colName_ak8jets;
  static const std::string colName_vtxs;
  static const std::string colName_pfcands;
  static const std::string colName_triggerinfos;
  static const std::string colName_triggerobjects;

protected:
  const edm::ParameterSet pset;
  std::shared_ptr<BaseTree> outtree;
  SimpleEntry commonEntry;
  bool firstEvent;

  int year;
  TString treename;
  //TString outfilename;
  bool isMC;
  bool is80X;

  bool processTriggerObjectInfos;
  bool keepMuonTimingInfo;

  std::string prefiringWeightsTag;
  bool applyPrefiringWeights;

  ParticleRecordLevel keepGenParticles;
  bool keepGenJets;

  bool const includeLJetsSelection;
  int const minNmuons;
  int const minNelectrons;
  int const minNleptons;
  int const minNphotons;
  int const minNak4jets;
  int const minNak8jets;

  edm::EDGetTokenT< edm::View<pat::Muon> > muonsToken;
  edm::EDGetTokenT< edm::View<pat::Electron> > electronsToken;
  edm::EDGetTokenT< edm::View<pat::Photon> > photonsToken;
  edm::EDGetTokenT< edm::View<pat::Jet> > ak4jetsToken;
  edm::EDGetTokenT< edm::View<pat::Jet> > ak8jetsToken;
  edm::EDGetTokenT< edm::View<IsotrackInfo> > isotracksToken;
  edm::EDGetTokenT< edm::View<pat::PackedCandidate> > pfcandsToken;

  edm::EDGetTokenT< METInfo > pfmetToken;
  edm::EDGetTokenT< reco::Particle::LorentzVector > pfmetshiftToken_JERNominal;
  edm::EDGetTokenT< reco::Particle::LorentzVector > pfmetshiftToken_JERUp;
  edm::EDGetTokenT< reco::Particle::LorentzVector > pfmetshiftToken_JERDn;
  edm::EDGetTokenT< METInfo > puppimetToken;

  edm::EDGetTokenT< reco::VertexCollection > vtxToken;

  edm::EDGetTokenT< double > rhoToken;
  edm::EDGetTokenT< edm::View<TriggerObjectInfo> > triggerObjectInfoToken;
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


  void recordGenInfo(edm::Event const&);
  void recordGenParticles(edm::Event const&, std::vector<reco::GenParticle const*>*, std::vector<pat::PackedGenParticle const*>*);
  void recordGenJets(edm::Event const&, bool const&, std::vector<reco::GenJet const*>*);

  size_t fillMuons(edm::Event const&, std::vector<pat::Muon const*>*);
  size_t fillElectrons(edm::Event const&, std::vector<pat::Electron const*>*);
  size_t fillFSRCandidates(
    edm::Event const&,
    std::vector<pat::Muon const*> const*, std::vector<pat::Electron const*> const*,
    std::vector<FSRCandidateInfo>*
  );
  size_t fillPhotons(edm::Event const&, std::vector<FSRCandidateInfo>*, std::vector<pat::Photon const*>*); // Not std::vector<FSRCandidateInfo> const* because the veto photon lists need to be modified.

  size_t fillAK4Jets(edm::Event const&, std::vector<pat::Jet const*>*);
  size_t fillAK8Jets(edm::Event const&, std::vector<pat::Jet const*>*);
  size_t fillIsotracks(edm::Event const&, std::vector<IsotrackInfo const*>*);

  size_t fillVertices(edm::Event const&, std::vector<reco::Vertex const*>*);

  bool fillEventVariables(edm::Event const&);
  bool fillTriggerInfo(edm::Event const&);
  bool fillMETFilterVariables(edm::Event const&);
  bool fillMETVariables(edm::Event const&);

  bool fillGenVariables(
    edm::Event const&,
    std::vector<reco::GenParticle const*>*,
    std::vector<pat::PackedGenParticle const*>*,
    std::vector<reco::GenJet const*>*,
    std::vector<reco::GenJet const*>*
  );

  static CMS3Ntuplizer::ParticleRecordLevel getParticleRecordLevel(std::string);

  template<typename T> void cleanUnusedCollection(bool const&, TString const&, T&);

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void analyze(edm::Event const&, const edm::EventSetup&);

};

template<typename T> void CMS3Ntuplizer::cleanUnusedCollection(bool const& isSelected, TString const& bname, T& vlist){
  if (!this->isMC || isSelected) return; // No need to clean, will not be recorded anyway

  if (
    bname.BeginsWith(CMS3Ntuplizer::colName_muons.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_electrons.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_photons.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_fsrcands.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_isotracks.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_ak4jets.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_ak8jets.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_vtxs.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_pfcands.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_triggerinfos.data())
    ||
    bname.BeginsWith(CMS3Ntuplizer::colName_triggerobjects.data())
    ){
    //std::cout << "CMS3Ntuplizer::cleanUnusedCollection: Collection " << bname << " can be cleaned because isMC=" << this->isMC << " and isSelected=" << isSelected << "." << std::endl;
    vlist.clear();
  }
}


#endif
