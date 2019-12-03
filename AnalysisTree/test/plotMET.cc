#include "common_includes.h"


struct SampleSpecs{
  std::string name;
  std::string label;
  std::string path;
  int mass;

  std::vector<TH1F> hlist_1D;

  SampleSpecs();
  SampleSpecs(std::string n, std::string l, std::string p, int m);
  SampleSpecs(const SampleSpecs& other);
};
SampleSpecs::SampleSpecs() :
  mass(-1)
{}
SampleSpecs::SampleSpecs(std::string n, std::string l, std::string p, int m):
  name(n),
  label(l),
  path(p),
  mass(m)
{}
SampleSpecs::SampleSpecs(const SampleSpecs& other):
  name(other.name),
  label(other.label),
  path(other.path),
  mass(other.mass),
  hlist_1D(other.hlist_1D)
{}


void plotMET(){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;

  TString const cinput_main = "/home/users/usarica/work/Width_AC_Run2/2L2Nu_MiniAOD/CMSSW_10_2_18/src/CMS3/NtupleMaker/test/manualruns";
  TString const coutput_main = "output";
  gSystem->mkdir(coutput_main, true);

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("ggHWW_M125", "gg#rightarrowH(125)#rightarrowWW", "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root", 125);
  sampleList.emplace_back("ggHWW_M500", "gg#rightarrowH(500)#rightarrowWW", "GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root", 500);
  sampleList.emplace_back("ggHWW_M3000", "gg#rightarrowH(3000)#rightarrowWW", "GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root", 3000);
  sampleList.emplace_back("ggHZZ_M200", "gg#rightarrowH(125)#rightarrowWW", "GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_0.root", 200);
  sampleList.emplace_back("ggHZZ_M500", "gg#rightarrowH(500)#rightarrowWW", "GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_0.root", 500);
  sampleList.emplace_back("ggHZZ_M3000", "gg#rightarrowH(3000)#rightarrowWW", "GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2_0.root", 3000);

  std::vector<float> genmetthresholds{ 50, 150, 300 };

  for (auto& sample:sampleList){
    BaseTree sample_tree(cinput_main + "/" + sample.path, "cms3ntuple/Events", "", "");

    TFile* foutput = TFile::Open(Form("%s/%s%s", coutput_main.Data(), sample.name.data(), ".root"), "recreate");

    GenInfoHandler genInfoHandler;
    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    MuonHandler muonHandler;
    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    ElectronHandler electronHandler;
    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    PhotonHandler photonHandler;
    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    JetMETHandler jetHandler;
    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    DileptonHandler dileptonHandler;

    MELAout << "Completed getting the handles..." << endl;

    foutput->cd();
    sample.hlist_1D.emplace_back(
      Form("h_%s_%s", sample.name.data(), "mll"),
      "m_{ll}",
      250, 0., 1000.
    );

    sample.hlist_1D.emplace_back(
      Form("h_%s_%s", sample.name.data(), "gen_pTmiss"),
      "p_{T}^{miss,true}",
      200, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s_%s", sample.name.data(), "reco_pfpuppi_pTmiss"),
      "p_{T}^{miss,PUPPI}",
      200, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s_%s", sample.name.data(), "reco_pfchs_pTmiss"),
      "p_{T}^{miss,CHS}",
      200, 0., 1000.
    );

    for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
      float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
      float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
      TString strgenmetbin;
      if (genmetlow==-1) strgenmetbin = Form("genmet_lt_%.0f", genmethigh);
      else if (genmethigh==-1) strgenmetbin = Form("genmet_gt_%.0f", genmetlow);
      else strgenmetbin = Form("genmet_%.0f_%.0f", genmetlow, genmethigh);

      sample.hlist_1D.emplace_back(
        Form("h_%s_%s_%s", sample.name.data(), "dPhi_pTll_pfpuppi_pTmiss", strgenmetbin.Data()),
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,PUPPI})",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s_%s", sample.name.data(), "resolution_pfpuppi_pTmiss", strgenmetbin.Data()),
        "p_{T}^{miss,PUPPI} resolution",
        40, -120, 120.
      );

      sample.hlist_1D.emplace_back(
        Form("h_%s_%s_%s", sample.name.data(), "dPhi_pTll_pfchs_pTmiss", strgenmetbin.Data()),
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,CHS})",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s_%s", sample.name.data(), "resolution_pfchs_pTmiss", strgenmetbin.Data()),
        "p_{T}^{miss,CHS} resolution",
        40, -120, 120.
      );
    }
    // Configure histograms
    for (TH1F& hh:sample.hlist_1D){
      hh.Sumw2();
    }

    const int nevents = sample_tree.getSelectedNEvents();
    for (int ev=0; ev<nevents; ev++){
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();

      muonHandler.constructMuons(theGlobalSyst);
      auto const& muons = muonHandler.getProducts();

      electronHandler.constructElectrons(theGlobalSyst);
      auto const& electrons = electronHandler.getProducts();

      photonHandler.constructPhotons(theGlobalSyst);
      auto const& photons = photonHandler.getProducts();

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfpuppimet = jetHandler.getPFPUPPIMET();
      auto const& pfchsmet = jetHandler.getPFCHSMET();
      ParticleObject::LorentzVector_t genmet;
      genmet = ParticleObject::PolarLorentzVector_t(genInfo->extras.genmet_met, 0, genInfo->extras.genmet_metPhi, 0);

      //MELAout << "MET values (PFPUPPI, PFCHS) = ( " << pfpuppimet->met() << ", " << pfchsmet->met() << " )" << endl;

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      //MELAout << "Ndileptons: " << dileptons.size() << " | pTs = ";
      //for (auto const& dilepton:dileptons) MELAout << dilepton->pt() << " ";
      //MELAout << endl;

      DileptonObject* theChosenDilepton = nullptr;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          theChosenDilepton = dilepton;
          break;
        }
      }

      if (theChosenDilepton){
        bool is_ee=false, is_mumu=false, is_emu=false;
        if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -121) is_ee=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -143) is_emu=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -169) is_mumu=true;

        float mll = theChosenDilepton->m();

        float dPhi_pTll_pfpuppi_pTmiss = theChosenDilepton->deltaPhi(pfpuppimet->phi());
        float resolution_pfpuppi_pTmiss = pfpuppimet->pt() - genmet.Pt();
        float reco_pfpuppi_pTmiss = pfpuppimet->pt();

        float dPhi_pTll_pfchs_pTmiss = theChosenDilepton->deltaPhi(pfchsmet->phi());
        float resolution_pfchs_pTmiss = pfchsmet->pt() - genmet.Pt();
        float reco_pfchs_pTmiss = pfchsmet->pt();

        float gen_pTmiss = genmet.Pt();

        float wgt = genInfo->extras.genHEPMCweight_default;
        {
          size_t ih=0;
          sample.hlist_1D.at(ih).Fill(mll, wgt);
          ih++;
          sample.hlist_1D.at(ih).Fill(gen_pTmiss, wgt);
          ih++;
          sample.hlist_1D.at(ih).Fill(reco_pfpuppi_pTmiss, wgt);
          ih++;
          sample.hlist_1D.at(ih).Fill(reco_pfchs_pTmiss, wgt);
          ih++;
          for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
            float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
            float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
            if (
              (genmetlow>=0. && gen_pTmiss<genmetlow)
              ||
              (genmethigh>=0. && gen_pTmiss>=genmethigh)
              ) continue;
            sample.hlist_1D.at(ih).Fill(dPhi_pTll_pfpuppi_pTmiss, wgt);
            ih++;
            sample.hlist_1D.at(ih).Fill(resolution_pfpuppi_pTmiss, wgt);
            ih++;
            sample.hlist_1D.at(ih).Fill(dPhi_pTll_pfchs_pTmiss, wgt);
            ih++;
            sample.hlist_1D.at(ih).Fill(resolution_pfchs_pTmiss, wgt);
            ih++;
          }
        }

      }

    } // End loop over events

    for (auto& hh:sample.hlist_1D) foutput->WriteTObject(&hh);
    foutput->Close();
  }
}

