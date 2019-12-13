#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"


struct HistogramObject{
  TString name;
  TString title;
  TString xlabel;
  TString ylabel;
  int nbins;
  float xlow;
  float xhigh;
  TDirectory* srcDir;

  TH1F hist;

  HistogramObject();
  HistogramObject(
    TString name_,
    TString title_,
    TString xlabel_,
    TString ylabel_,
    int nbins_,
    float xlow_,
    float xhigh_,
    TDirectory* srcDir_=nullptr
  );
  HistogramObject(HistogramObject const& other);

  void write();
};

HistogramObject::HistogramObject() :
  nbins(-1),
  xlow(0),
  xhigh(0),
  srcDir(nullptr)
{}
HistogramObject::HistogramObject(
  TString name_,
  TString title_,
  TString xlabel_,
  TString ylabel_,
  int nbins_,
  float xlow_,
  float xhigh_,
  TDirectory* srcDir_
) :
  name(name_),
  title(title_),
  xlabel(xlabel_),
  ylabel(ylabel_),
  nbins(nbins_),
  xlow(xlow_),
  xhigh(xhigh_),
  srcDir(srcDir_),
  hist(name, title, nbins, xlow, xhigh)
{
  hist.GetXaxis()->SetTitle(xlabel);
  hist.GetYaxis()->SetTitle(ylabel);

  MELAout << "Created histogram " << hist.GetName() << " [" << hist.GetTitle() << "]" << endl;
}
HistogramObject::HistogramObject(HistogramObject const& other) :
  name(other.name),
  title(other.title),
  xlabel(other.xlabel),
  ylabel(other.ylabel),
  nbins(other.nbins),
  xlow(other.xlow),
  xhigh(other.xhigh),
  srcDir(other.srcDir),
  hist(other.hist)
{}
void HistogramObject::write(){
  if (srcDir){
    TDirectory* curdir = gDirectory;
    srcDir->cd();
    srcDir->WriteTObject(&hist);
    curdir->cd();
  }
}



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
  void writeHistograms();
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

    hist.GetXaxis()->SetNdivisions(505);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetLabelOffset(0.007);
    hist.GetXaxis()->SetLabelSize(0.04);
    hist.GetXaxis()->SetTitleSize(0.06);
    hist.GetXaxis()->SetTitleOffset(0.9);
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetNdivisions(505);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelOffset(0.007);
    hist.GetYaxis()->SetLabelSize(0.04);
    hist.GetYaxis()->SetTitleSize(0.06);
    hist.GetYaxis()->SetTitleOffset(1.1);
    hist.GetYaxis()->SetTitleFont(42);
  }
}
void SampleSpecs::writeHistograms(){ for (auto& hh:hlist_1D) hh.write(); }


void getChannelTitleLabel(int ichannel, TString& title, TString& label){
  if (ichannel==0){
    title = "ee";
    label = "ee";
  }
  else if (ichannel==1){
    title = "mumu";
    label = "#mu#mu";
  }
  else if (ichannel==2){
    title = "emu";
    label = "e#mu";
  }
  else{
    title="AllChannels";
    label = "ee+#mu#mu+e#mu";
  }
}

void getHistograms_ZZCuts(int doZZWW, int procsel, TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  constexpr int nchannels = 4; // ichannel=-1, 0, 1, 2 for any, ee, mumu, emu

  TString const cinput_main = "/home/users/usarica/work/Width_AC_Run2/Samples/101124";
  TString const coutput_main = "output/" + strdate + (doZZWW==0 ? "/ZZCuts" : "/WWCuts");

  gSystem->mkdir(coutput_main, true);

  std::vector<SampleSpecs> sampleList;
  if (procsel<0 || procsel==0){
    sampleList.emplace_back("DY_M10-50", "DY ll (m_{ll}=10-50 GeV)", "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", -1, HistogramProperties((int) kGreen+2, 1, 2));
    sampleList.emplace_back("DY_M50", "DY ll (m_{ll}>50 GeV)", "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", -1, HistogramProperties((int) kCyan, 1, 2));
    sampleList.emplace_back("TT2L2Nu", "t#bar{t} ll", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", -1, HistogramProperties((int) kOrange-3, 1, 2));
    sampleList.emplace_back("ZZ2L2Nu", "ZZ#rightarrow2l2#nu", "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/allevents.root", -1, HistogramProperties((int) kYellow-3, 1, 2));
    sampleList.emplace_back("WW2L2Nu", "WW#rightarrow2l2#nu", "WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", -1, HistogramProperties((int) kTeal-1, 1, 2));
  }
  if (procsel<0 || procsel==1){
    sampleList.emplace_back("ggHWW_M125", "ggH(125)#rightarrowWW", "GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 125, HistogramProperties((int) kViolet, 2, 2));
    sampleList.emplace_back("ggHWW_M500", "ggH(500)#rightarrowWW", "GluGluHToWWTo2L2Nu_M500_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 500, HistogramProperties((int) kBlue, 2, 2));
    sampleList.emplace_back("ggHWW_M3000", "ggH(3000)#rightarrowWW", "GluGluHToWWTo2L2Nu_M3000_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 3000, HistogramProperties((int) kRed, 2, 2));
    sampleList.emplace_back("ggHZZ_M200", "ggH(200)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/allevents.root", 200, HistogramProperties((int) kViolet, 7, 2));
    sampleList.emplace_back("ggHZZ_M500", "ggH(500)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", 500, HistogramProperties((int) kBlue, 7, 2));
    sampleList.emplace_back("ggHZZ_M3000", "ggH(3000)#rightarrowZZ", "GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/allevents.root", 3000, HistogramProperties((int) kRed, 7, 2));
  }

  //std::vector<float> genmetthresholds{ 50, 150, 300 };
  std::vector<float> genmetthresholds;

  for (auto& sample:sampleList){
    BaseTree sample_tree(cinput_main + "/" + sample.path, "cms3ntuple/Events", "", "");

    TFile* foutput = TFile::Open(Form("%s/%s%s", coutput_main.Data(), sample.name.data(), ".root"), "recreate");
    MELAout.open(Form("%s/%s%s", coutput_main.Data(), sample.name.data(), ".txt"));

    // Get cross section
    sample_tree.bookBranch<float>("xsec", 0.f);

    // Get handlers
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

    //EventFilterHandler eventFilter;
    //eventFilter.bookBranches(&sample_tree);
    //eventFilter.wrapTree(&sample_tree);

    DileptonHandler dileptonHandler;

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    foutput->cd();

    for (int ichannel=-1; ichannel<nchannels-1; ichannel++){
      TString strChannel, strChannelLabel;
      getChannelTitleLabel(ichannel, strChannel, strChannelLabel);

      TDirectory* channeldir = foutput->mkdir(strChannel);
      channeldir->cd();

      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "gen_pTmiss"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,true} (GeV)", "",
        160, 0., 800.,
        channeldir
      );

      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "mll"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "m_{ll} (GeV)", "",
        250, 0., 1000.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "pTll"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{ll} (GeV)", "",
        250, 0., 1000.,
        channeldir
      );

      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "n_ak4jets_tight"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "N_{j}", "",
        6, 0, 6,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "n_ak4jets_tight_btagged"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "N_{b}", "",
        4, 0, 4,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "mjj"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "m_{j1j2} (GeV)", "",
        125, 0., 1000.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "mjets"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "m_{jets} (GeV)", "",
        125, 0., 1000.,
        channeldir
      );

      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfpuppi_pTmiss"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,PUPPI} (GeV)", "",
        160, 0., 800.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfpuppi_pTmiss_significance"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,PUPPI} significance", "",
        100, 0., 10000.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfpuppi_pTmiss_over_pTll"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,PUPPI}/p_{T}^{ll}", "",
        100, 0., 20.,
        channeldir
      );

      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfchs_pTmiss"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,CHS} (GeV)", "",
        160, 0., 800.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfchs_pTmiss_significance"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,CHS} significance", "",
        75, 0., 150.,
        channeldir
      );
      sample.hlist_1D.emplace_back(
        Form("h1D_%s_%s", strChannel.Data(), "reco_pfchs_pTmiss_over_pTll"), Form("%s|%s", sample.label.data(), strChannelLabel.Data()),
        "p_{T}^{miss,CHS}/p_{T}^{ll}", "",
        100, 0., 20.,
        channeldir
      );

      for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
        float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
        float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
        TString strgenmetbin;
        if (genmetlow==-1 && genmethigh==-1) strgenmetbin = "genmet_inclusive";
        else if (genmetlow==-1) strgenmetbin = Form("genmet_lt_%.0f", genmethigh);
        else if (genmethigh==-1) strgenmetbin = Form("genmet_gt_%.0f", genmetlow);
        else strgenmetbin = Form("genmet_%.0f_%.0f", genmetlow, genmethigh);
        TString strgenmetbinlabel;
        if (genmetlow==-1 && genmethigh==-1) strgenmetbinlabel = "Inclusive p_{T}^{miss,true}";
        else if (genmetlow==-1) strgenmetbinlabel = Form("p_{T}^{miss,true} < %.0f GeV", genmethigh);
        else if (genmethigh==-1) strgenmetbinlabel = Form("p_{T}^{miss,true} > %.0f GeV", genmetlow);
        else strgenmetbinlabel = Form("p_{T}^{miss,true}: [%.0f, %.0f] GeV", genmetlow, genmethigh);

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dEta_j1j2_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2"),
          "|#Delta#eta(#vec{p}_{j1}, #vec{p}_{j2})|", "",
          50, 0., +5.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dR_highest_btagval_jets_Nj_ge_2_Nb_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2, N_{b}>=1"),
          "#DeltaR_{jj} (highest b-tag val.)", "",
          50, 0., +5.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dR_highest_btagval_nonbtagged_jets_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2, N_{b}=0"),
          "#DeltaR_{jj} (highest b-tag val.)", "",
          50, 0., +5.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_pfpuppi_pTmiss_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTll_pfpuppi_pTmiss", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "|#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTll_pfpuppi_pTmiss_Nj_eq_0", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}=0"),
          "|#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTllj_pfpuppi_pTmiss_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "|#Delta#phi(#vec{p}_{T}^{ll+j1}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTlljj_pfpuppi_pTmiss_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2"),
          "|#Delta#phi(#vec{p}_{T}^{ll+j1j2}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTlljets_pfpuppi_pTmiss_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2"),
          "|#Delta#phi(#vec{p}_{T}^{ll+jets}, #vec{p}_{T}^{miss,PUPPI})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "resolution_pfpuppi_pTmiss", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{miss,PUPPI}/p_{T}^{miss,true} - 1", "",
          50, -1., 4.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "reco_pfpuppi_pTmiss_significance", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{miss,PUPPI} significance", "",
          50, 0., 200.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mTZZ_pfpuppi", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "m_{T}^{ZZ,PUPPI} (GeV)", "",
          300, 0., 3000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mZZ_plus_pfpuppi", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "m^{ZZ,PUPPI} (GeV)", "",
          300, 0., 3000.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_pfchs_pTmiss_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTll_pfchs_pTmiss", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "|#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTll_pfchs_pTmiss_Nj_eq_0", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}=0"),
          "|#Delta#phi(#vec{p}_{T}^{ll}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTllj_pfchs_pTmiss_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "|#Delta#phi(#vec{p}_{T}^{ll+j1}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTlljj_pfchs_pTmiss_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2"),
          "|#Delta#phi(#vec{p}_{T}^{ll+j1j2}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "dPhi_pTlljets_pfchs_pTmiss_Nj_ge_2", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=2"),
          "|#Delta#phi(#vec{p}_{T}^{ll+jets}, #vec{p}_{T}^{miss,CHS})|", "",
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "resolution_pfchs_pTmiss", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{miss,CHS}/p_{T}^{miss,true} - 1", "",
          50, -1., 4.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "reco_pfchs_pTmiss_significance", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{miss,CHS} significance", "",
          50, 0., 200.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mTZZ_pfchs", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "m_{T}^{ZZ,CHS} (GeV)", "",
          300, 0., 3000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mZZ_plus_pfchs", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "m^{ZZ,CHS} (GeV)", "",
          300, 0., 3000.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "ml1j1_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "m_{l1j1}", "",
          100, 0., 1000.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mlj1_min_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "m_{lj1} (min. over l_{1}, l_{2})", "",
          100, 0., 1000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mlj1_closest_lj_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "m_{lj1} (closest lepton)", "",
          100, 0., 1000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mlj_min_best_b_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "m_{lj} (highest medium b-tag, min. over l_{1}, l_{2})", "",
          100, 0., 1000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mlj_closest_l_best_b_Nj_ge_1", strgenmetbin.Data()), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data(), "N_{j}>=1"),
          "m_{lj} (highest medium b-tag, closest lepton)", "",
          100, 0., 1000.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "pTl1", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{l1}", "",
          200, 0., 800.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "pTl2", strgenmetbin.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), strgenmetbinlabel.Data()),
          "p_{T}^{l2}", "",
          200, 0., 800.,
          channeldir
        );
      }

      foutput->cd();
    }

    // Configure histograms
    sample.setup();

    float xsec=-1;
    double sum_wgts=0;
    double sum_ee=0, sum_mumu=0, sum_emu=0;
    double sum_ee_selected=0, sum_mumu_selected=0, sum_emu_selected=0;
    const int nevents = sample_tree.getSelectedNEvents();
    for (int ev=0; ev<nevents; ev++){
      HelperFunctions::progressbar(ev, nevents);
      //if (ev>1000) break;
      sample_tree.getSelectedEvent(ev);
      if (ev==0){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
      }

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();

      float wgt = genInfo->getGenWeight(true);
      sum_wgts += wgt;

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
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=0.4184) ak4jets_tight_btagged.push_back(jet);
        }
      }

      //MELAout << "MET values (PFPUPPI, PFCHS) = ( " << pfpuppimet->met() << ", " << pfchsmet->met() << " )" << endl;

      //eventFilter.constructFilters();

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

        if (is_ee) sum_ee += wgt;
        else if (is_mumu) sum_mumu += wgt;
        else if (is_emu) sum_emu += wgt;

        float mll = theChosenDilepton->m();
        float pTll = theChosenDilepton->pt();
        float gen_pTmiss = genmet.Pt(); if (gen_pTmiss==0.f) gen_pTmiss = 1e-5;

        ParticleObject* leadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(0) : theChosenDilepton->daughter(1));
        ParticleObject* subleadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(1) : theChosenDilepton->daughter(0));
        ParticleObject::LorentzVector_t p4_l1 = leadingLepton->p4();
        ParticleObject::LorentzVector_t p4_l2 = subleadingLepton->p4();
        float pTl1 = leadingLepton->pt();
        float pTl2 = subleadingLepton->pt();

        size_t n_ak4jets_tight = ak4jets_tight.size();
        size_t n_ak4jets_tight_btagged = ak4jets_tight_btagged.size();

        ParticleObject::LorentzVector_t p4_l1j1 = p4_l1;
        ParticleObject::LorentzVector_t p4_j;
        ParticleObject::LorentzVector_t p4_jj;
        ParticleObject::LorentzVector_t p4_alljets;
        AK4JetObject* highest_btag_jet = nullptr;
        AK4JetObject* secondHighest_btag_jet = nullptr;
        float highest_btagvalue = -99;
        AK4JetObject* highest_nonbtagged_jet = nullptr;
        AK4JetObject* secondHighest_nonbtagged_jet = nullptr;
        float highest_nonbtagged_jet_btagvalue = -99;
        for (size_t ijet=0; ijet<n_ak4jets_tight; ijet++){
          AK4JetObject* jet = ak4jets_tight.at(ijet);
          if (ijet<1){
            p4_j = p4_j + jet->p4();
            p4_l1j1 = p4_l1j1 + jet->p4();
          }
          if (ijet<2) p4_jj = p4_jj + jet->p4();
          p4_alljets = p4_alljets + jet->p4();

          float btagval = jet->getBtagValue();
          if (btagval<0.4184 && std::abs(jet->eta())<2.5){
            if (!highest_nonbtagged_jet || highest_nonbtagged_jet_btagvalue<btagval){
              secondHighest_nonbtagged_jet = highest_nonbtagged_jet;
              highest_nonbtagged_jet = jet;
              highest_nonbtagged_jet_btagvalue = btagval;
            }
          }
          if (btagval>0.){
            if (!highest_btag_jet || highest_btagvalue<btagval){
              secondHighest_btag_jet = highest_btag_jet;
              highest_btag_jet = jet;
              highest_btagvalue = btagval;
            }
          }
        }
        float mjj = p4_jj.M();
        float mjets = p4_alljets.M();
        float pZmiss_approx = -(theChosenDilepton->p4()+p4_jj).Z();
        float ml1j1 = p4_l1j1.M();
        float dEta_j1j2=-99; if (n_ak4jets_tight>=2) dEta_j1j2 = (p4_j.eta() - ak4jets_tight.at(1)->eta());

        float dR_highest_btagval_jets = -1;
        if (secondHighest_btag_jet) dR_highest_btagval_jets = reco::deltaR(highest_btag_jet->p4(), secondHighest_btag_jet->p4());
        float dR_highest_btagval_nonbtagged_jets = -1;
        if (secondHighest_nonbtagged_jet) dR_highest_btagval_nonbtagged_jets = reco::deltaR(highest_nonbtagged_jet->p4(), secondHighest_nonbtagged_jet->p4());

        float m_lj1_min = -1;
        float m_lj1_closest = -1;
        ParticleObject* lepton_closest_to_j1 = nullptr;
        if (n_ak4jets_tight>0){
          lepton_closest_to_j1 = (reco::deltaR(p4_l1, p4_j)<reco::deltaR(p4_l2, p4_j) ? leadingLepton : subleadingLepton);
          m_lj1_min = std::min((p4_l1 + p4_j).M(), (p4_l2 + p4_j).M());
          m_lj1_closest = (lepton_closest_to_j1->p4() + p4_j).M();
        }

        float m_lj_min_best_b = -1;
        float m_lj_closest_l_best_b = -1;
        ParticleObject* lepton_closest_to_best_b = nullptr;
        if (highest_nonbtagged_jet){
          lepton_closest_to_best_b = (reco::deltaR(p4_l1, highest_nonbtagged_jet->p4())<reco::deltaR(p4_l2, highest_nonbtagged_jet->p4()) ? leadingLepton : subleadingLepton);
          m_lj_min_best_b = std::min((p4_l1 + highest_nonbtagged_jet->p4()).M(), (p4_l2 + highest_nonbtagged_jet->p4()).M());
          m_lj_closest_l_best_b = (lepton_closest_to_best_b->p4() + highest_nonbtagged_jet->p4()).M();
        }

        float abs_dPhi_min_pTj_pfpuppi_pTmiss = TMath::Pi();
        for (AK4JetObject* jet:ak4jets_tight){
          float dphi_tmp; HelperFunctions::deltaPhi(float(jet->phi()), float(pfpuppimet->phi()), dphi_tmp);
          abs_dPhi_min_pTj_pfpuppi_pTmiss = std::min(abs_dPhi_min_pTj_pfpuppi_pTmiss, std::abs(dphi_tmp));
        }
        float dPhi_pTll_pfpuppi_pTmiss = theChosenDilepton->deltaPhi(pfpuppimet->phi());
        float dPhi_pTllj_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_j).Phi()), float(pfpuppimet->phi()), dPhi_pTllj_pfpuppi_pTmiss);
        float dPhi_pTlljj_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_jj).Phi()), float(pfpuppimet->phi()), dPhi_pTlljj_pfpuppi_pTmiss);
        float dPhi_pTlljets_pfpuppi_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfpuppimet->phi()), dPhi_pTlljets_pfpuppi_pTmiss);
        float reco_pfpuppi_pTmiss = pfpuppimet->pt();
        float reco_pfpuppi_pTmiss_significance = pfpuppimet->extras.metSignificance;
        float resolution_pfpuppi_pTmiss = reco_pfpuppi_pTmiss/gen_pTmiss - 1.;
        // Compute ZZ-style masses
        float mTZZ_pfpuppi = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(reco_pfpuppi_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfpuppimet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t pfpuppi_p4_ZZplusapprox; pfpuppi_p4_ZZplusapprox = ParticleObject::PolarLorentzVector_t(reco_pfpuppi_pTmiss, -std::asinh(pZmiss_approx/reco_pfpuppi_pTmiss), pfpuppimet->phi(), PDGHelpers::Zmass);
        float mZZ_plus_pfpuppi = (pfpuppi_p4_ZZplusapprox + theChosenDilepton->p4()).M();

        float abs_dPhi_min_pTj_pfchs_pTmiss = TMath::Pi();
        for (AK4JetObject* jet:ak4jets_tight){
          float dphi_tmp; HelperFunctions::deltaPhi(float(jet->phi()), float(pfchsmet->phi()), dphi_tmp);
          abs_dPhi_min_pTj_pfchs_pTmiss = std::min(abs_dPhi_min_pTj_pfchs_pTmiss, std::abs(dphi_tmp));
        }
        float dPhi_pTll_pfchs_pTmiss = theChosenDilepton->deltaPhi(pfchsmet->phi());
        float dPhi_pTllj_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_j).Phi()), float(pfchsmet->phi()), dPhi_pTllj_pfchs_pTmiss);
        float dPhi_pTlljj_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_jj).Phi()), float(pfchsmet->phi()), dPhi_pTlljj_pfchs_pTmiss);
        float dPhi_pTlljets_pfchs_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfchsmet->phi()), dPhi_pTlljets_pfchs_pTmiss);
        float reco_pfchs_pTmiss = pfchsmet->pt();
        float reco_pfchs_pTmiss_significance = pfchsmet->extras.metSignificance;
        float resolution_pfchs_pTmiss = reco_pfchs_pTmiss/gen_pTmiss - 1.;
        // Compute ZZ-style masses
        float mTZZ_pfchs = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(reco_pfchs_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfchsmet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t pfchs_p4_ZZplusapprox; pfchs_p4_ZZplusapprox = ParticleObject::PolarLorentzVector_t(reco_pfchs_pTmiss, -std::asinh(pZmiss_approx/reco_pfchs_pTmiss), pfchsmet->phi(), PDGHelpers::Zmass);
        float mZZ_plus_pfchs = (pfchs_p4_ZZplusapprox + theChosenDilepton->p4()).M();

        // Cuts
        bool pass_pTl1 = pTl1>=25.;
        bool pass_pTl2 = pTl1>=20.;
        bool pass_pTll = pTll>=20.;
        bool pass_Nb_veto = n_ak4jets_tight_btagged==0;
        bool pass_mll = (doZZWW==0 && mll>=70. && mll<105.) || (doZZWW==1 && mll>=105.);
        bool pass_dPhi_j_pTmiss = /*n_ak4jets_tight==0 || abs_dPhi_min_pTj_pfpuppi_pTmiss>=0.3*/true;
        float pTmiss_over_pTll_ratio_thr = 0.4;
        if (doZZWW==0){
          pTmiss_over_pTll_ratio_thr = std::min(0.6f, std::max(0.4f, float(pTll<210.f ? (pTll - 40.f)/(210.f - 40.f)*0.1f + 0.4f : (pTll - 210.f)/(400.f - 210.f)*0.1 + 0.5f)));
        }
        bool pass_pTmiss_over_pTll_ratio = (reco_pfpuppi_pTmiss/pTll>pTmiss_over_pTll_ratio_thr && reco_pfpuppi_pTmiss/pTll<1./pTmiss_over_pTll_ratio_thr);
        bool pass_pTmiss_significance = reco_pfpuppi_pTmiss_significance>40. || doZZWW==1;

        if (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio){
          if (is_ee) sum_ee_selected += wgt;
          else if (is_mumu) sum_mumu_selected += wgt;
          else if (is_emu) sum_emu_selected += wgt;
        }

        // Fill histograms
        // Enclosed around braces to localize it_hist
        {
          auto it_hist = sample.hlist_1D.begin();
          for (int ichannel=-1; ichannel<nchannels-1; ichannel++){
            bool isCorrectChannel = (
              ichannel==-1
              || (ichannel==0 && is_ee)
              || (ichannel==1 && is_mumu)
              || (ichannel==2 && is_emu)
              );

            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(gen_pTmiss, wgt); it_hist++;

            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(mll, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(pTll, wgt); it_hist++;

            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(n_ak4jets_tight, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(n_ak4jets_tight_btagged, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(mjj, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(mjets, wgt); it_hist++;

            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfpuppi_pTmiss, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfpuppi_pTmiss_significance, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss)) it_hist->hist.Fill(reco_pfpuppi_pTmiss/pTll, wgt); it_hist++;

            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfchs_pTmiss, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfchs_pTmiss_significance, wgt); it_hist++;
            if (isCorrectChannel && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss)) it_hist->hist.Fill(reco_pfchs_pTmiss/pTll, wgt); it_hist++;
            for (size_t igenmetbin=0; igenmetbin<=genmetthresholds.size(); igenmetbin++){
              float genmetlow = (igenmetbin==0 ? -1 : genmetthresholds.at(igenmetbin-1));
              float genmethigh = (igenmetbin==genmetthresholds.size() ? -1 : genmetthresholds.at(igenmetbin));
              bool doFill = !(
                (genmetlow>=0. && gen_pTmiss<genmetlow)
                ||
                (genmethigh>=0. && gen_pTmiss>=genmethigh)
                );

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(std::abs(dEta_j1j2), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && !pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && secondHighest_btag_jet) it_hist->hist.Fill(std::abs(dR_highest_btagval_jets), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && secondHighest_nonbtagged_jet) it_hist->hist.Fill(std::abs(dR_highest_btagval_nonbtagged_jets), wgt); it_hist++;

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=1) it_hist->hist.Fill(abs_dPhi_min_pTj_pfpuppi_pTmiss, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(std::abs(dPhi_pTll_pfpuppi_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight==0) it_hist->hist.Fill(std::abs(dPhi_pTll_pfpuppi_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=1) it_hist->hist.Fill(std::abs(dPhi_pTllj_pfpuppi_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(std::abs(dPhi_pTlljj_pfpuppi_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(std::abs(dPhi_pTlljets_pfpuppi_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(resolution_pfpuppi_pTmiss, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfpuppi_pTmiss_significance, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(mTZZ_pfpuppi, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(mZZ_plus_pfpuppi, wgt); it_hist++;

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=1) it_hist->hist.Fill(abs_dPhi_min_pTj_pfchs_pTmiss, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(std::abs(dPhi_pTll_pfchs_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight==0) it_hist->hist.Fill(std::abs(dPhi_pTll_pfchs_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=1) it_hist->hist.Fill(std::abs(dPhi_pTllj_pfchs_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(std::abs(dPhi_pTlljj_pfchs_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=2) it_hist->hist.Fill(std::abs(dPhi_pTlljets_pfchs_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(resolution_pfchs_pTmiss, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(reco_pfchs_pTmiss_significance, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(mTZZ_pfchs, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(mZZ_plus_pfchs, wgt); it_hist++;

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && n_ak4jets_tight>=1) it_hist->hist.Fill(ml1j1, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && lepton_closest_to_j1) it_hist->hist.Fill(m_lj1_min, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && lepton_closest_to_j1) it_hist->hist.Fill(m_lj1_closest, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && lepton_closest_to_best_b) it_hist->hist.Fill(m_lj_min_best_b, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio) && lepton_closest_to_best_b) it_hist->hist.Fill(m_lj_closest_l_best_b, wgt); it_hist++;


              if (isCorrectChannel && doFill && (pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(pTl1, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_significance && pass_dPhi_j_pTmiss && pass_pTmiss_over_pTll_ratio)) it_hist->hist.Fill(pTl2, wgt); it_hist++;
            }
          } // End loop over channels
        } // End fill

      }

    } // End loop over events

    constexpr bool normByXsec=true;
    for (auto& hh:sample.hlist_1D){
      TH1F& hist = hh.hist;
      if (normByXsec){
        if (sum_wgts>0.) hist.Scale(xsec*1000.*59.7/sum_wgts);
      }
      else{
        double integral=0;
        for (int ix=1; ix<=hist.GetNbinsX(); ix++){
          double bc = hist.GetBinContent(ix);
          integral += bc;
        }
        if (integral>0.) hist.Scale(1./integral);
      }
    }
    sample.writeHistograms();

    MELAout
      << "xsec: " << xsec << '\n'
      << "sum_wgts: " << sum_wgts << '\n'
      << "xsec_ee: " << sum_ee * xsec / sum_wgts << '\n'
      << "xsec_mumu: " << sum_mumu * xsec / sum_wgts << '\n'
      << "xsec_emu: " << sum_emu * xsec / sum_wgts << '\n'
      << "xsec_ee_selected: " << sum_ee_selected * xsec / sum_wgts << '\n'
      << "xsec_mumu_selected: " << sum_mumu_selected * xsec / sum_wgts << '\n'
      << "xsec_emu_selected: " << sum_emu_selected * xsec / sum_wgts << '\n'
      << "frac_ee: " << sum_ee / sum_wgts << '\n'
      << "frac_mumu: " << sum_mumu / sum_wgts << '\n'
      << "frac_emu: " << sum_emu / sum_wgts << '\n'
      << "frac_ee_selected: " << sum_ee_selected / sum_wgts << '\n'
      << "frac_mumu_selected: " << sum_mumu_selected / sum_wgts << '\n'
      << "frac_emu_selected: " << sum_emu_selected / sum_wgts << '\n'
      << endl;

    MELAout.close();
    foutput->Close();
  } // End loop over samples
}


void plotMET_ZZCuts(int doZZWW, TString strdate="", bool useLogY=false, bool isStacked=false, int nfoci=0, int ifocus=0){
  if (nfoci>0 && (ifocus>=nfoci || ifocus<0)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  //SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  //constexpr int nchannels = 4;

  TString const cinput_main = "output/" + strdate + (doZZWW==0 ? "/ZZCuts" : "/WWCuts");
  TString coutput_main = cinput_main + "/Plots";
  if (isStacked) coutput_main += "/Stacked";
  else coutput_main += "/Unstacked";
  if (useLogY) coutput_main += "/LogY";
  else coutput_main += "/Linear";
  if (nfoci>0) coutput_main += Form("/Divide%i/Focus%i", nfoci, ifocus);
  {
    // Special case to copy index.php if you have one
    std::vector<TString> tmplist;
    HelperFunctions::splitOptionRecursive(coutput_main, tmplist, '/');
    TString indexDir = "${CMSSW_BASE}/src/CMSDataTools/AnalysisTree/data/plotting/index.php";
    HostHelpers::ExpandEnvironmentVariables(indexDir);
    if (HostHelpers::FileReadable(indexDir)){
      MELAout << "Attempting to copy index.php" << endl;
      TString tmpdir = tmplist.at(0) + '/';
      for (size_t idir=1; idir<tmplist.size(); idir++){
        tmpdir = tmpdir + tmplist.at(idir) + '/';
        TString tmpCmd = "cp ~/public_html/index.pages.php ";
        tmpCmd += tmpdir + "index.php";
        MELAout << "Copying index.php into " << tmpdir << endl;
        HostHelpers::ExecuteCommand(tmpCmd);
      }
    }
  }


  gSystem->mkdir(coutput_main, true);

  std::vector<std::string> sampleList;
  sampleList.emplace_back("DY_M10-50");
  sampleList.emplace_back("DY_M50");
  sampleList.emplace_back("TT2L2Nu");
  sampleList.emplace_back("WW2L2Nu");
  sampleList.emplace_back("ZZ2L2Nu");
  sampleList.emplace_back("ggHWW_M125");
  sampleList.emplace_back("ggHWW_M500");
  sampleList.emplace_back("ggHWW_M3000");
  sampleList.emplace_back("ggHZZ_M200");
  sampleList.emplace_back("ggHZZ_M500");
  sampleList.emplace_back("ggHZZ_M3000");

  std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;

  bool firstFile = true;
  std::vector<TString> hnames;
  std::vector<TFile*> finputlist;
  for (auto const& sample:sampleList){
    TFile* finput = TFile::Open(cinput_main + '/' + sample.data() + ".root", "read");
    finputlist.push_back(finput);

    std::vector<TH1F*> hlist;
    HelperFunctions::extractHistogramsFromDirectory(finput, hlist);
    sample_hist_map[sample] = hlist;

    if (firstFile){
      for (TH1F* hh:hlist) hnames.push_back(hh->GetName());
      firstFile=false;
    }
  }

  for (size_t iplot=0; iplot<hnames.size(); iplot++){
    int nbins=sample_hist_map[sampleList.front()].at(0)->GetNbinsX();
    int binstart=1;
    int binend=nbins;
    if (nfoci>0){
      int bininc = nbins / nfoci;
      binstart = (ifocus * bininc) + 1;
      binend = (ifocus==nfoci-1 ? nbins : binstart + bininc);
    }
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);
      hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinLowEdge(binstart), hist->GetXaxis()->GetBinUpEdge(binend));
    }
    if (isStacked){
      for (size_t is=0; is<sampleList.size(); is++){
        auto& sample = sampleList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        bool isSignal = sample.find("ggH")!=std::string::npos;
        if (!isSignal){
          for (size_t js=is+1; js<sampleList.size(); js++){
            auto& sample_j = sampleList.at(js);
            TH1F* hist_j = sample_hist_map[sample_j].at(iplot);
            bool isSignal_j = sample_j.find("ggH")!=std::string::npos;
            if (isSignal_j) continue;
            hist->Add(hist_j);
            //MELAout << "Adding " << sample_j << " to " << sample << endl;
          }
        }
        /*
        else{
          TH1F* accbkg = sample_hist_map[sampleList.front()].at(iplot);
          double intbkg = accbkg->Integral(binstart, binend);
          double intsig = hist->Integral(binstart, binend);
          if (intsig>0. && intbkg>0.) hist->Scale(intbkg / intsig);
        }
        */
      }
    }
    else{
      for (size_t is=0; is<sampleList.size(); is++){
        auto& sample = sampleList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        double inthist = hist->Integral(binstart, binend);
        hist->Scale(1. / inthist);
      }
    }

    size_t nplottables=0;
    double ymin = 0;
    double ymax = -1;
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (sample.find("ZZ")!=std::string::npos && TString(hist->GetName()).Contains("emu")) continue;
      if (doZZWW==0 && TString(hist->GetName()).Contains("emu")) continue;
      if (sample.find("HZZ")!=std::string::npos && doZZWW==1) continue;
      if (sample.find("HWW")!=std::string::npos && doZZWW==0) continue;
      nplottables++;

      for (int ix=1; ix<=hist->GetNbinsX(); ix++){
        double bc = hist->GetBinContent(ix);
        double be = hist->GetBinError(ix);
        if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
        ymax = std::max(ymax, bc+be);
        if (useLogY && bc>0.){
          if (ymin<=0.) ymin = bc;
          else ymin = std::min(ymin, bc);
        }
      }
    }
    ymax *= (useLogY ? 15. : 1.5);
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);
      hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    TString canvasname = hnames.at(iplot);
    TString canvasname_core = (useLogY ? (isStacked ? "cLogY_Stacked_" : "cLogY_") : (isStacked ? "c_Stacked_" : "c_"));
    if (nfoci>0) canvasname_core += Form("Focus_%i_of_%i_", ifocus, nfoci);
    HelperFunctions::replaceString(canvasname, "h1D_", "");
    canvasname = canvasname_core + canvasname;
    TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
    canvas->cd();
    gStyle->SetOptStat(0);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(2);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLeftMargin(0.17);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.07);
    canvas->SetBottomMargin(0.13);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);
    if (useLogY) canvas->SetLogy();

    TLegend* legend = new TLegend(
      0.55,
      0.90-0.10/4.*2.*float(nplottables),
      0.90,
      0.90
    );
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColor(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    TText* text;

    TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);
    text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
    text->SetTextSize(0.044);
    text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", 59.7, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (sample.find("ZZ")!=std::string::npos && TString(hist->GetName()).Contains("emu")) continue;

      std::vector<TString> tmplist;
      TString htitle = hist->GetTitle();
      HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
      TString hlabel = tmplist.front();

      hist->SetTitle("");
      legend->AddEntry(hist, hlabel, "l");

      if (firstHist){
        for (size_t is=1; is<tmplist.size(); is++) selectionList.push_back(tmplist.at(is));
        hist->Draw("hist");
        firstHist = false;
      }
      else{
        hist->Draw("histsame");
      }
    }

    legend->Draw("same");
    pt->Draw();

    std::vector<TPaveText*> ptSelectionList;
    {
      float pt_ymax = 0.90;
      float pt_dy = 0.05;
      for (auto const& strSel:selectionList){
        TPaveText* ptSel = nullptr;
        ptSel = new TPaveText(0.20, pt_ymax - pt_dy, 0.50, pt_ymax, "brNDC");
        ptSel->SetBorderSize(0);
        ptSel->SetFillStyle(0);
        ptSel->SetTextAlign(12);
        ptSel->SetTextFont(42);
        ptSel->SetTextSize(0.045);
        text = ptSel->AddText(0.025, 0.45, strSel);
        text->SetTextSize(0.0315);
        ptSel->Draw();

        ptSelectionList.push_back(ptSel);
        pt_ymax -= pt_dy;
      }
    }

    canvas->RedrawAxis();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(coutput_main + "/" + canvasname + ".pdf");
    canvas->SaveAs(coutput_main + "/" + canvasname + ".png");

    for (auto*& ptSel:ptSelectionList) delete ptSel;
    delete pt;
    delete legend;
    canvas->Close();
  }

  for (TFile*& finput:finputlist) finput->Close();
}
