#include <cassert>
#include <limits>
#include "common_includes.h"


using namespace std;


void produceBtaggingEfficiencies(TString strSampleSet, TString period, BtagHelpers::BtagWPType btagwptype, TString strdate=""){
  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "191212");

  constexpr SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;

  BtagHelpers::setBtagWPType(btagwptype);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle
    }
  );

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  VertexHandler vertexHandler;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  EventFilterHandler eventFilter;
  DileptonHandler dileptonHandler;

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampleList);

  TString const coutput_main = "output/" + strdate + "/BtaggingEffs";
  gSystem->mkdir(coutput_main, true);

  MELAout << "List of samples to process: " << sampleList << endl;
  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData) continue;

    TString cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    float sum_wgts=0;
    const int nEntries = sample_tree.getSelectedNEvents();

    simEventHandler.bookBranches(&sample_tree);
    simEventHandler.wrapTree(&sample_tree);

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

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

    sample_tree.silenceUnused();

    MELAout << "Completed getting the rest of the handles..." << endl;
    sample_tree.silenceUnused();

    // Create output tree
    TString cSample = SampleHelpers::getDatasetDirectoryCoreName(strSample.Data()).data();
    HelperFunctions::replaceString(cSample, "/", "_");
    TString stroutput = coutput_main + "/" + cSample + "_btageff_" + SystematicsHelpers::getSystName(theGlobalSyst).data() + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    ExtendedBinning ptbins({ 0, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000 });
    ExtendedBinning etabins(24, -2.4, 2.4);
    etabins.addBinBoundary(-2.5); etabins.addBinBoundary(2.5);
    etabins.addBinBoundary(-5.); etabins.addBinBoundary(5.);
    std::vector<TH2F> h_All{
      TH2F("AllJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("AllJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };
    std::vector<TH2F> h_DeepFlavor_Loose{
      TH2F("DeepFlavor_LooseJets_b", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_LooseJets_c", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning()),
      TH2F("DeepFlavor_LooseJets_udsg", "", ptbins.getNbins(), ptbins.getBinning(), etabins.getNbins(), etabins.getBinning())
    };

    // Loop over the tree
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    unsigned int n_acc=0;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();
      float genwgt = genInfo->getGenWeight(true);

      simEventHandler.constructSimEvent(theGlobalSyst);
      float puwgt = simEventHandler.getPileUpWeight();

      float wgt = genwgt * puwgt;
      sum_wgts += wgt;

      eventFilter.constructFilters();
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float trigwgt = eventFilter.getTriggerWeight(triggerCheckList);
      wgt *= trigwgt;

      if (wgt==0.f) continue;

      muonHandler.constructMuons(theGlobalSyst);
      auto const& muons = muonHandler.getProducts();

      electronHandler.constructElectrons(theGlobalSyst);
      auto const& electrons = electronHandler.getProducts();

      photonHandler.constructPhotons(theGlobalSyst);
      auto const& photons = photonHandler.getProducts();

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();

      if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons, &photons, &ak4jets, &ak8jets)) continue;

      ParticleObject::LorentzVector_t p4_alljets;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto const& jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          p4_alljets = p4_alljets + jet->p4();
        }
      }
      size_t n_ak4jets_tight = ak4jets_tight.size();
      if (n_ak4jets_tight==0) continue;

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }
      if (nTightDilep!=1) continue;

      bool is_ee=false, is_mumu=false, is_emu=false;
      if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -121) is_ee=true;
      else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -143) is_emu=true;
      else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -169) is_mumu=true;
      if (is_emu) continue;

      ParticleObject* leadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(0) : theChosenDilepton->daughter(1));
      ParticleObject* subleadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(1) : theChosenDilepton->daughter(0));
      float pTl1 = leadingLepton->pt();
      float pTl2 = subleadingLepton->pt();
      float pTll = theChosenDilepton->pt();

      float dPhi_pTlljets_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(puppimet->phi()), dPhi_pTlljets_puppimet_pTmiss);

      bool pass_pTl1 = pTl1>=25.;
      bool pass_pTl2 = pTl2>=20.;
      bool pass_pTll = pTll>=35.;
      bool pass_puppimet_thr = puppimet->pt()>=85.;
      bool pass_puppimet_dPhilljets_thr = (std::abs(dPhi_pTlljets_puppimet_pTmiss)>=2.6);
      if (!(pass_pTl1 && pass_pTl2 && pass_pTll && pass_puppimet_thr && pass_puppimet_dPhilljets_thr)) continue;

      for (auto const& jet:ak4jets_tight){
        float btagval = jet->getBtagValue();
        int const& jetFlavor = jet->extras.hadronFlavour;

        unsigned int iflav = 2;
        if (abs(jetFlavor)==5) iflav = 0;
        else if (abs(jetFlavor)==4) iflav = 1;
        h_All.at(iflav).Fill(jet->pt(), jet->eta(), wgt);
        if (btagval>=btagvalue_thr) h_DeepFlavor_Loose.at(iflav).Fill(jet->pt(), jet->eta(), wgt);
      }
    }

    for (unsigned int iflav=0; iflav<3; iflav++){
      h_DeepFlavor_Loose.at(iflav).Divide(&(h_All.at(iflav)));
      foutput->WriteTObject(&h_DeepFlavor_Loose.at(iflav));
      foutput->WriteTObject(&h_All.at(iflav));
    }
    foutput->Close();
  }
}
