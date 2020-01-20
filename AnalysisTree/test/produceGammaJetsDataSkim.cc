#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void produceGammaJetsDataSkim(TString strSampleSet, TString period){
  SampleHelpers::configure(period, "hadoop:200101");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle,
      OffshellTriggerHelpers::kSinglePho
    }
  );

  // Get handlers
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  ParticleDisambiguator particleDisambiguator;

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  TString const stroutputcore = Form("output/GammaJetsDataSkim/%s", SampleHelpers::theDataPeriod.Data());

  bool isFirstSample = true;
  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (!isData) continue;

    TString cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();

    // Create output
    TString stroutput = stroutputcore + "/" + strSample;
    gSystem->Exec(Form("mkdir -p %s", stroutput.Data()));
    stroutput += "/allevents.root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TDirectory* outdir = foutput->mkdir("cms3ntuple");
    outdir->cd();
    sample_tree.unmuteAllBranches();
    TTree* outtree = sample_tree.getSelectedTree()->CloneTree(0);

    // Loop over the tree
    const int nEntries = sample_tree.getSelectedNEvents();
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    size_t n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      if (ev%10000==0) MELAout << sample_tree.sampleIdentifier << " events: " << n_acc << " / " << ev << " / " << nEntries << endl;

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      if (trigwgt==0.) continue;

      muonHandler.constructMuons(SystematicsHelpers::sNominal);
      electronHandler.constructElectrons(SystematicsHelpers::sNominal);
      photonHandler.constructPhotons(SystematicsHelpers::sNominal);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      size_t n_muons_veto = 0;
      for (auto const& part:muons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_muons_veto++; }
      if (n_muons_veto>0) continue;

      auto const& electrons = electronHandler.getProducts();
      size_t n_electrons_veto = 0;
      for (auto const& part:electrons){ if (ParticleSelectionHelpers::isVetoParticle(part)) n_electrons_veto++; }
      if (n_electrons_veto>0) continue;

      auto const& photons = photonHandler.getProducts();
      size_t n_photons_tight = 0;
      PhotonObject* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }
      if (n_photons_tight!=1) continue;

      jetHandler.constructJetMET(SystematicsHelpers::sNominal, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();

      bool hasBtaggedJet=false;
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          if (jet->getBtagValue()>=btagvalue_thr){
            hasBtaggedJet=true;
            break;
          }
        }
      }
      if (hasBtaggedJet) continue;

      if (!eventFilter.test2018HEMFilter(nullptr, &electrons, &photons, &ak4jets, &ak8jets)) continue;

      outtree->Fill();
      n_acc++;
    }

    outdir->WriteTObject(outtree);
    delete outtree;
    foutput->cd();
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  }
}
