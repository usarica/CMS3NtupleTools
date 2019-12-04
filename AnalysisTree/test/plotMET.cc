#include "common_includes.h"


struct HistogramObject{
  TString name;
  TString title;
  TString xlabel;
  TString ylabel;
  int nbins;
  float xlow;
  float xhigh;

  TH1F hist;

  HistogramObject();
  HistogramObject(
    TString name_,
    TString title_,
    TString xlabel_,
    TString ylabel_,
    int nbins_,
    float xlow_,
    float xhigh_
  );
  HistogramObject(HistogramObject const& other);

};

HistogramObject::HistogramObject() :
  nbins(-1),
  xlow(0),
  xhigh(0)
{}
HistogramObject::HistogramObject(
  TString name_,
  TString title_,
  TString xlabel_,
  TString ylabel_,
  int nbins_,
  float xlow_,
  float xhigh_
) :
  name(name_),
  title(title_),
  xlabel(xlabel_),
  ylabel(ylabel_),
  nbins(nbins_),
  xlow(xlow_),
  xhigh(xhigh_),
  hist(name, title, nbins, xlow, xhigh)
{
  hist.GetXaxis()->SetTitle(xlabel);
  hist.GetYaxis()->SetTitle(ylabel);
}
HistogramObject::HistogramObject(HistogramObject const& other) :
  name(other.name),
  title(other.title),
  xlabel(other.xlabel),
  ylabel(other.ylabel),
  nbins(other.nbins),
  xlow(other.xlow),
  xhigh(other.xhigh),
  hist(other.hist)
{}


struct HistogramProperties{
  int color;
  int dashtype;
  int width;

  HistogramProperties();
  HistogramProperties(int c, int d, int w);
  HistogramProperties(HistogramProperties const& other);
};
HistogramProperties::HistogramProperties() :
  color((int) kBlack),
  dashtype(1),
  width(2)
{}
HistogramProperties::HistogramProperties(int c, int d, int w) :
  color(c),
  dashtype(d),
  width(w)
{}
HistogramProperties::HistogramProperties(HistogramProperties const& other) :
  color(other.color),
  dashtype(other.dashtype),
  width(other.width)
{}


struct SampleSpecs{
  std::string name;
  std::string label;
  std::string path;
  int mass;
  HistogramProperties props;

  std::vector<HistogramObject> hlist_1D;

  SampleSpecs();
  SampleSpecs(std::string n, std::string l, std::string p, int m, HistogramProperties pp = HistogramProperties());
  SampleSpecs(const SampleSpecs& other);

  void setup();
};
SampleSpecs::SampleSpecs() :
  mass(-1)
{}
SampleSpecs::SampleSpecs(std::string n, std::string l, std::string p, int m, HistogramProperties pp):
  name(n),
  label(l),
  path(p),
  props(pp),
  mass(m)
{}
SampleSpecs::SampleSpecs(const SampleSpecs& other):
  name(other.name),
  label(other.label),
  path(other.path),
  mass(other.mass),
  props(other.props),
  hlist_1D(other.hlist_1D)
{}
void SampleSpecs::setup(){
  for (auto& hh:hlist_1D){
    TH1F& hist = hh.hist;

    hist.Sumw2();

    hist.SetLineColor(props.color);
    hist.SetMarkerColor(props.color);
    hist.SetLineWidth(props.width);
    hist.SetLineStyle(props.dashtype);
  }
}



void plotMET(){
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;

  TString const cinput_main = "/home/users/usarica/work/Width_AC_Run2/Samples/101124";
  TString const coutput_main = "output";
  gSystem->mkdir(coutput_main, true);

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("ggHWW_M125", "gg#rightarrowH(125)#rightarrowWW", "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 125, HistogramProperties((int)kViolet, 2, 2));
  sampleList.emplace_back("ggHWW_M500", "gg#rightarrowH(500)#rightarrowWW", "GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 500, HistogramProperties((int) kBlue, 2, 2));
  sampleList.emplace_back("ggHWW_M3000", "gg#rightarrowH(3000)#rightarrowWW", "GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 3000, HistogramProperties((int) kRed, 2, 2));
  sampleList.emplace_back("ggHZZ_M200", "gg#rightarrowH(125)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 200, HistogramProperties((int) kViolet, 7, 2));
  sampleList.emplace_back("ggHZZ_M500", "gg#rightarrowH(500)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", 500, HistogramProperties((int) kBlue, 7, 2));
  sampleList.emplace_back("ggHZZ_M3000", "gg#rightarrowH(3000)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", 3000, HistogramProperties((int) kRed, 7, 2));
  sampleList.emplace_back("DY_M10-50", "DY ll (m_{ll}=10-50 GeV)", "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", -1, HistogramProperties((int) kGreen+2, 1, 2));
  sampleList.emplace_back("DY_M50", "DY ll (m_{ll}>50 GeV)", "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("TT2L2Nu", "t#bar{t} ll", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", -1, HistogramProperties((int) kOrange-3, 1, 2));

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
    sample_tree.silenceUnused();

    foutput->cd();

    sample.hlist_1D.emplace_back(
      Form("h_%s", "gen_pTmiss"), "",
      "p_{T}^{miss,true} (GeV)", "",
      200, 0., 1000.
    );

    sample.hlist_1D.emplace_back(
      Form("h_%s", "mll"), "",
      "m_{ll} (GeV)", "",
      250, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "pTll"), "",
      "p_{T}^{ll} (GeV)", "",
      250, 0., 1000.
    );

    sample.hlist_1D.emplace_back(
      Form("h_%s", "n_ak4jets_tight"), "",
      "N_{j}", "",
      6, 0, 6
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "mjj"), "",
      "m_{j1j2} (GeV)", "",
      250, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "mjets"), "",
      "m_{jets} (GeV)", "",
      250, 0., 1000.
    );

    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfpuppi_pTmiss"), "",
      "p_{T}^{miss,PUPPI} (GeV)", "",
      200, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfpuppi_pTmiss_significance"), "",
      "p_{T}^{miss,PUPPI} significance", "",
      100, 0., 100.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfpuppi_pTmiss_over_pTll"), "",
      "p_{T}^{miss,PUPPI}/p_{T}^{ll}", "",
      100, 0., 20.
    );

    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfchs_pTmiss"), "",
      "p_{T}^{miss,CHS} (GeV)", "",
      200, 0., 1000.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfchs_pTmiss_significance"), "",
      "p_{T}^{miss,CHS} significance", "",
      100, 0., 100.
    );
    sample.hlist_1D.emplace_back(
      Form("h_%s", "reco_pfchs_pTmiss_over_pTll"), "",
      "p_{T}^{miss,CHS}/p_{T}^{ll}", "",
      100, 0., 20.
    );

    for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
      float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
      float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
      TString strgenmetbin;
      if (genmetlow==-1) strgenmetbin = Form("genmet_lt_%.0f", genmethigh);
      else if (genmethigh==-1) strgenmetbin = Form("genmet_gt_%.0f", genmetlow);
      else strgenmetbin = Form("genmet_%.0f_%.0f", genmetlow, genmethigh);
      TString strgenmetbinlabel;
      if (genmetlow==-1) strgenmetbinlabel = Form("p_{T}^{miss,true} < %.0f GeV", genmethigh);
      else if (genmethigh==-1) strgenmetbinlabel = Form("p_{T}^{miss,true} > %.0f GeV", genmetlow);
      else strgenmetbinlabel = Form("p_{T}^{miss,true}: [%.0f, %.0f]", genmetlow, genmethigh);

      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTll_pfpuppi_pTmiss", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,PUPPI})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTll_pfpuppi_pTmiss_Nj_eq_0", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,PUPPI})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTllj_pfpuppi_pTmiss_Nj_ge_1", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+j1}, #vec{p}_{T}^{miss,PUPPI})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTlljj_pfpuppi_pTmiss_Nj_ge_2", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+j1j2}, #vec{p}_{T}^{miss,PUPPI})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTlljets_pfpuppi_pTmiss_Nj_ge_2", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+jets}, #vec{p}_{T}^{miss,PUPPI})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "resolution_pfpuppi_pTmiss", strgenmetbin.Data()), strgenmetbinlabel,
        "p_{T}^{miss,PUPPI} resolution (GeV)", "",
        40, -120, 120.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "reco_pfpuppi_pTmiss_significance", strgenmetbin.Data()), strgenmetbinlabel,
        "p_{T}^{miss,PUPPI} significance", "",
        100, 0., 100.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mTZZ_pfpuppi", strgenmetbin.Data()), strgenmetbinlabel,
        "m_{T}^{ZZ,PUPPI} (GeV)", "",
        300, 0., 3000.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mZZ_plus_pfpuppi", strgenmetbin.Data()), strgenmetbinlabel,
        "m^{ZZ+,PUPPI} (GeV)", "",
        300, 0., 3000.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mZZ_minus_pfpuppi", strgenmetbin.Data()), strgenmetbinlabel,
        "m^{ZZ-,PUPPI} (GeV)", "",
        300, 0., 3000.
      );

      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTll_pfchs_pTmiss", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,CHS})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTll_pfchs_pTmiss_Nj_eq_0", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,CHS})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTllj_pfchs_pTmiss_Nj_ge_1", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+j1}, #vec{p}_{T}^{miss,CHS})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTlljj_pfchs_pTmiss_Nj_ge_2", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+j1j2}, #vec{p}_{T}^{miss,CHS})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "dPhi_pTlljets_pfchs_pTmiss_Nj_ge_2", strgenmetbin.Data()), strgenmetbinlabel,
        "#Delta#phi(#vec{p}_{T}^{ll+jets}, #vec{p}_{T}^{miss,CHS})", "",
        40, -TMath::Pi(), +TMath::Pi()
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "resolution_pfchs_pTmiss", strgenmetbin.Data()), strgenmetbinlabel,
        "p_{T}^{miss,CHS} resolution (GeV)", "",
        40, -120, 120.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "reco_pfchs_pTmiss_significance", strgenmetbin.Data()), strgenmetbinlabel,
        "p_{T}^{miss,CHS} significance", "",
        100, 0., 100.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mTZZ_pfchs", strgenmetbin.Data()), strgenmetbinlabel,
        "m_{T}^{ZZ,CHS} (GeV)", "",
        300, 0., 3000.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mZZ_plus_pfchs", strgenmetbin.Data()), strgenmetbinlabel,
        "m^{ZZ+,CHS} (GeV)", "",
        300, 0., 3000.
      );
      sample.hlist_1D.emplace_back(
        Form("h_%s_%s", "mZZ_minus_pfchs", strgenmetbin.Data()), strgenmetbinlabel,
        "m^{ZZ-,CHS} (GeV)", "",
        300, 0., 3000.
      );
    }

    // Configure histograms
    sample.setup();

    const int nevents = sample_tree.getSelectedNEvents();
    for (int ev=0; ev<nevents; ev++){
      HelperFunctions::progressbar(ev, nevents);
      if (ev>10000) break;
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

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto* jet:ak4jets){ if (ParticleSelectionHelpers::isTightJet(jet)) ak4jets_tight.push_back(jet); }

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
        float pTll = theChosenDilepton->pt();
        float gen_pTmiss = genmet.Pt();

        size_t n_ak4jets_tight = ak4jets_tight.size();
        ParticleObject::LorentzVector_t p4_j;
        ParticleObject::LorentzVector_t p4_jj;
        ParticleObject::LorentzVector_t p4_alljets;
        for (size_t ijet=0; ijet<n_ak4jets_tight; ijet++){
          if (ijet<1) p4_j = p4_j + ak4jets_tight.at(ijet)->p4();
          if (ijet<2) p4_jj = p4_jj + ak4jets_tight.at(ijet)->p4();
          p4_alljets = p4_alljets + ak4jets_tight.at(ijet)->p4();
        }
        float mjj = p4_jj.M();
        float mjets = p4_alljets.M();
        float pZmiss_approx = -(theChosenDilepton->p4()+p4_jj).Z();

        float dPhi_pTll_pfpuppi_pTmiss = theChosenDilepton->deltaPhi(pfpuppimet->phi());
        float dPhi_pTllj_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_j).Phi()), float(pfpuppimet->phi()), dPhi_pTllj_pfpuppi_pTmiss);
        float dPhi_pTlljj_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_jj).Phi()), float(pfpuppimet->phi()), dPhi_pTlljj_pfpuppi_pTmiss);
        float dPhi_pTlljets_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfpuppimet->phi()), dPhi_pTlljets_pfpuppi_pTmiss);
        float resolution_pfpuppi_pTmiss = pfpuppimet->pt() - genmet.Pt();
        float reco_pfpuppi_pTmiss = pfpuppimet->pt();
        float reco_pfpuppi_pTmiss_significance = pfpuppimet->extras.metSignificance;
        // Compute ZZ-style mTs
        float mTZZ_pfpuppi = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(reco_pfpuppi_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfpuppimet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t pfpuppi_p4_ZZplusapprox; pfpuppi_p4_ZZplusapprox = ParticleObject::PolarLorentzVector_t(reco_pfpuppi_pTmiss, -std::asinh(pZmiss_approx/reco_pfpuppi_pTmiss), pfpuppimet->phi(), PDGHelpers::Zmass);
        float mZZ_plus_pfpuppi = (pfpuppi_p4_ZZplusapprox + theChosenDilepton->p4()).M();
        ParticleObject::LorentzVector_t pfpuppi_p4_ZZminusapprox; pfpuppi_p4_ZZminusapprox = ParticleObject::PolarLorentzVector_t(reco_pfpuppi_pTmiss, +std::asinh(pZmiss_approx/reco_pfpuppi_pTmiss), pfpuppimet->phi(), PDGHelpers::Zmass);
        float mZZ_minus_pfpuppi = (pfpuppi_p4_ZZminusapprox + theChosenDilepton->p4()).M();

        float dPhi_pTll_pfchs_pTmiss = theChosenDilepton->deltaPhi(pfchsmet->phi());
        float dPhi_pTllj_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_j).Phi()), float(pfchsmet->phi()), dPhi_pTllj_pfchs_pTmiss);
        float dPhi_pTlljj_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_jj).Phi()), float(pfchsmet->phi()), dPhi_pTlljj_pfchs_pTmiss);
        float dPhi_pTlljets_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfchsmet->phi()), dPhi_pTlljets_pfchs_pTmiss);
        float resolution_pfchs_pTmiss = pfchsmet->pt() - genmet.Pt();
        float reco_pfchs_pTmiss = pfchsmet->pt();
        float reco_pfchs_pTmiss_significance = pfchsmet->extras.metSignificance;
        // Compute ZZ-style mTs
        float mTZZ_pfchs = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(reco_pfchs_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfchsmet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t pfchs_p4_ZZplusapprox; pfchs_p4_ZZplusapprox = ParticleObject::PolarLorentzVector_t(reco_pfchs_pTmiss, -std::asinh(pZmiss_approx/reco_pfchs_pTmiss), pfchsmet->phi(), PDGHelpers::Zmass);
        float mZZ_plus_pfchs = (pfchs_p4_ZZplusapprox + theChosenDilepton->p4()).M();
        ParticleObject::LorentzVector_t pfchs_p4_ZZminusapprox; pfchs_p4_ZZminusapprox = ParticleObject::PolarLorentzVector_t(reco_pfchs_pTmiss, +std::asinh(pZmiss_approx/reco_pfchs_pTmiss), pfchsmet->phi(), PDGHelpers::Zmass);
        float mZZ_minus_pfchs = (pfchs_p4_ZZminusapprox + theChosenDilepton->p4()).M();

        float wgt = genInfo->getGenWeight(true);
        {
          size_t ih=0;
          sample.hlist_1D.at(ih).hist.Fill(gen_pTmiss, wgt); ih++;

          sample.hlist_1D.at(ih).hist.Fill(mll, wgt); ih++;
          sample.hlist_1D.at(ih).hist.Fill(pTll, wgt); ih++;

          sample.hlist_1D.at(ih).hist.Fill(n_ak4jets_tight, wgt); ih++;
          if (n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(mjj, wgt); ih++;
          if (n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(mjets, wgt); ih++;

          sample.hlist_1D.at(ih).hist.Fill(reco_pfpuppi_pTmiss, wgt); ih++;
          sample.hlist_1D.at(ih).hist.Fill(reco_pfpuppi_pTmiss_significance, wgt); ih++;
          sample.hlist_1D.at(ih).hist.Fill(reco_pfpuppi_pTmiss/pTll, wgt); ih++;

          sample.hlist_1D.at(ih).hist.Fill(reco_pfchs_pTmiss, wgt); ih++;
          sample.hlist_1D.at(ih).hist.Fill(reco_pfchs_pTmiss_significance, wgt); ih++;
          sample.hlist_1D.at(ih).hist.Fill(reco_pfchs_pTmiss/pTll, wgt); ih++;
          for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
            float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
            float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
            bool doFill = !(
              (genmetlow>=0. && gen_pTmiss<genmetlow)
              ||
              (genmethigh>=0. && gen_pTmiss>=genmethigh)
              );

            if (doFill) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTll_pfpuppi_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight==0) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTll_pfpuppi_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=1) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTllj_pfpuppi_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTlljj_pfpuppi_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTlljets_pfpuppi_pTmiss, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(resolution_pfpuppi_pTmiss, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(reco_pfpuppi_pTmiss_significance, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mTZZ_pfpuppi, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mZZ_plus_pfpuppi, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mZZ_minus_pfpuppi, wgt); ih++;

            if (doFill) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTll_pfchs_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight==0) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTll_pfchs_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=1) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTllj_pfchs_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTlljj_pfchs_pTmiss, wgt); ih++;
            if (doFill && n_ak4jets_tight>=2) sample.hlist_1D.at(ih).hist.Fill(dPhi_pTlljets_pfchs_pTmiss, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(resolution_pfchs_pTmiss, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(reco_pfchs_pTmiss_significance, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mTZZ_pfchs, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mZZ_plus_pfchs, wgt); ih++;
            if (doFill) sample.hlist_1D.at(ih).hist.Fill(mZZ_minus_pfchs, wgt); ih++;
          }
        }

      }

    } // End loop over events

    for (auto& hh:sample.hlist_1D) foutput->WriteTObject(&(hh.hist));
    foutput->Close();
  }
}

