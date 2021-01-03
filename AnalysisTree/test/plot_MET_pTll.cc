#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TF2.h"


struct HistogramObject_2D{
  TString name;
  TString title;
  TString xlabel;
  TString ylabel;
  TString zlabel;
  int nbinsx;
  float xlow;
  float xhigh;
  int nbinsy;
  float ylow;
  float yhigh;
  TDirectory* srcDir;

  TH2F hist;

  HistogramObject_2D();
  HistogramObject_2D(
    TString name_,
    TString title_,
    TString xlabel_,
    TString ylabel_,
    TString zlabel_,
    int nbinsx_,
    float xlow_,
    float xhigh_,
    int nbinsy_,
    float ylow_,
    float yhigh_,
    TDirectory* srcDir_=nullptr
  );
  HistogramObject_2D(HistogramObject_2D const& other);

  void write();
};

HistogramObject_2D::HistogramObject_2D() :
  nbinsx(-1),
  xlow(0),
  xhigh(0),
  nbinsy(-1),
  ylow(0),
  yhigh(0),
  srcDir(nullptr)
{}
HistogramObject_2D::HistogramObject_2D(
  TString name_,
  TString title_,
  TString xlabel_,
  TString ylabel_,
  TString zlabel_,
  int nbinsx_,
  float xlow_,
  float xhigh_,
  int nbinsy_,
  float ylow_,
  float yhigh_,
  TDirectory* srcDir_
) :
  name(name_),
  title(title_),
  xlabel(xlabel_),
  ylabel(ylabel_),
  zlabel(zlabel_),
  nbinsx(nbinsx_),
  xlow(xlow_),
  xhigh(xhigh_),
  nbinsy(nbinsy_),
  ylow(ylow_),
  yhigh(yhigh_),
  srcDir(srcDir_),
  hist(name, title, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh)
{
  hist.GetXaxis()->SetTitle(xlabel);
  hist.GetYaxis()->SetTitle(ylabel);
  hist.GetZaxis()->SetTitle(zlabel);

  MELAout << "Created histogram " << hist.GetName() << " [" << hist.GetTitle() << "]" << endl;
}
HistogramObject_2D::HistogramObject_2D(HistogramObject_2D const& other) :
  name(other.name),
  title(other.title),
  xlabel(other.xlabel),
  ylabel(other.ylabel),
  zlabel(other.zlabel),
  nbinsx(other.nbinsx),
  xlow(other.xlow),
  xhigh(other.xhigh),
  nbinsy(other.nbinsy),
  ylow(other.ylow),
  yhigh(other.yhigh),
  srcDir(other.srcDir),
  hist(other.hist)
{}
void HistogramObject_2D::write(){
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

  std::vector<HistogramObject_2D> hlist_2D;

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
  hlist_2D(other.hlist_2D)
{}
void SampleSpecs::setup(){
  for (auto& hh:hlist_2D){
    TH2F& hist = hh.hist;

    hist.Sumw2();
    hist.SetOption("colz");

    hist.SetLineColor(props.color);
    hist.SetMarkerColor(props.color);
    /*
    if (mass<0){
      hist.SetFillColor(props.color);
      hist.SetFillStyle(3018);
    }
    */

    hist.SetLineWidth(props.width);
    hist.SetLineStyle(props.dashtype);

    hist.GetXaxis()->SetNdivisions(505);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetLabelOffset(0.007);
    hist.GetXaxis()->SetLabelSize(0.0315);
    hist.GetXaxis()->SetTitleSize(0.04);
    hist.GetXaxis()->SetTitleOffset(0.9);
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetNdivisions(505);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelOffset(0.007);
    hist.GetYaxis()->SetLabelSize(0.0315);
    hist.GetYaxis()->SetTitleSize(0.04);
    hist.GetYaxis()->SetTitleOffset(1.1);
    hist.GetYaxis()->SetTitleFont(42);
    hist.GetZaxis()->SetNdivisions(505);
    hist.GetZaxis()->SetLabelFont(42);
    hist.GetZaxis()->SetLabelOffset(0.007);
    hist.GetZaxis()->SetLabelSize(0.0315);
    hist.GetZaxis()->SetTitleSize(0.04);
    hist.GetZaxis()->SetTitleOffset(1.1);
    hist.GetZaxis()->SetTitleFont(42);
  }
}
void SampleSpecs::writeHistograms(){ for (auto& hh:hlist_2D) hh.write(); }

struct CutSpecs{
  TString cutvar;
  TString cutvarlabel;
  bool doCutLow;
  bool doCutHigh;
  float cutlow;
  float cuthigh;

  CutSpecs(
    TString cutvar_="",
    TString cutvarlabel_="",
    bool doCutLow_=false,
    bool doCutHigh_=false,
    float cutlow_=-1,
    float cuthigh_=-1
  );
  CutSpecs(const CutSpecs& other);

  bool testCut(float var) const;

  TString getTitle() const;
  TString getLabel() const;

  static TString getCutValLabelString(float cutval);
  static TString getCutValTitleString(float cutval);
};
CutSpecs::CutSpecs(
  TString cutvar_,
  TString cutvarlabel_,
  bool doCutLow_,
  bool doCutHigh_,
  float cutlow_,
  float cuthigh_
) : 
  cutvar(cutvar_),
  cutvarlabel(cutvarlabel_),
  doCutLow(doCutLow_),
  doCutHigh(doCutHigh_),
  cutlow(cutlow_),
  cuthigh(cuthigh_)
{}
CutSpecs::CutSpecs(const CutSpecs& other) :
  cutvar(other.cutvar),
  cutvarlabel(other.cutvarlabel),
  doCutLow(other.doCutLow),
  doCutHigh(other.doCutHigh),
  cutlow(other.cutlow),
  cuthigh(other.cuthigh)
{}
TString CutSpecs::getLabel() const{
  if (!doCutLow && !doCutHigh) return Form("Inclusive %s", cutvarlabel.Data());
  else if (doCutLow && doCutHigh && cutlow==cuthigh) return Form("%s=%s", cutvarlabel.Data(), CutSpecs::getCutValLabelString(cutlow).Data());
  else if (doCutLow && doCutHigh) return Form("%s #in [%s, %s)", cutvarlabel.Data(), CutSpecs::getCutValLabelString(cutlow).Data(), CutSpecs::getCutValLabelString(cuthigh).Data());
  else if (doCutLow) return Form("%s #geq %s", cutvarlabel.Data(), CutSpecs::getCutValLabelString(cutlow).Data());
  else return Form("%s < %s", cutvarlabel.Data(), CutSpecs::getCutValLabelString(cuthigh).Data());
}
TString CutSpecs::getTitle() const{
  if (!doCutLow && !doCutHigh) return Form("%s_inclusive", cutvar.Data());
  else if (doCutLow && doCutHigh && cutlow==cuthigh) return Form("%s_eq_%s", cutvar.Data(), CutSpecs::getCutValTitleString(cutlow).Data());
  else if (doCutLow && doCutHigh) return Form("%s_in_%s_%s", cutvar.Data(), CutSpecs::getCutValTitleString(cutlow).Data(), CutSpecs::getCutValTitleString(cuthigh).Data());
  else if (doCutLow) return Form("%s_ge_%s", cutvar.Data(), CutSpecs::getCutValTitleString(cutlow).Data());
  else return Form("%s_lt_%s", cutvar.Data(), CutSpecs::getCutValTitleString(cuthigh).Data());
}
TString CutSpecs::getCutValLabelString(float cutval){
  float decimals = std::abs(cutval - float((int) cutval));
  if (decimals == 0.f) return Form("%.0f", cutval);
  int base10exponent = std::ceil(std::abs(std::log10(decimals)));
  TString strprintf = Form("%s%i%s", "$.", base10exponent, "f");
  return Form(strprintf.Data(), cutval);
}
TString CutSpecs::getCutValTitleString(float cutval){
  TString label = getCutValLabelString(cutval);
  HelperFunctions::replaceString(label, ".", "p");
  return label;
}
bool CutSpecs::testCut(float var) const{
  if (!doCutLow && !doCutHigh) return true;
  else if (doCutLow && doCutHigh && cutlow==cuthigh) return (var==cutlow);
  else if (doCutLow && doCutHigh) return (var>=cutlow && var<cuthigh);
  else if (doCutLow) return (var>=cutlow);
  else return (var<cuthigh);
}



void getChannelTitleLabel(int ichannel, TString& title, TString& label){
  if (ichannel==0){
    title = "OS_ee";
    label = "OS ee";
  }
  else if (ichannel==1){
    title = "OS_mumu";
    label = "OS #mu#mu";
  }
  else if (ichannel==2){
    title = "OS_emu";
    label = "OS e#mu";
  }
  else if (ichannel==3){
    title = "OS_eeORmumu";
    label = "OS ee or #mu#mu";
  }
  else{
    title="AllChannels";
    label = "OS ee, #mu#mu, or e#mu";
  }
}

void getHistograms(int doZZWW, int procsel, TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  constexpr int nchannels = 4; // ichannel=-1, 0, 1, 2, 3 for any, ee, mumu, emu, ee+mumu

  TString const coutput_main = "output/" + strdate + "/MET_vs_pTll" + (doZZWW==0 ? "/ZZCuts" : "/WWCuts");
  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", ((procsel>=11) ? "191212" : "store:200101"));

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);
  float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle
    }
  );

  std::vector<SampleSpecs> sampleList;
  //sampleList.emplace_back("DY_M10-50", "DY ll (m_{ll}=10-50 GeV)", "DY_2l_M_10to50", -1, HistogramProperties((int) kGreen+2, 1, 2));
  //sampleList.emplace_back("DY_M50", "DY ll (m_{ll}>50 GeV)", "DY_2l_M_50", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_70-100", "DY ll (m_{ll}>50 GeV, H_{T}: 70-100 GeV)", "DY_2l_M_50_HT_70-100", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_100-200", "DY ll (m_{ll}>50 GeV, H_{T}: 100-200 GeV)", "DY_2l_M_50_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_200-400", "DY ll (m_{ll}>50 GeV, H_{T}: 200-400 GeV)", "DY_2l_M_50_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_400-600", "DY ll (m_{ll}>50 GeV, H_{T}: 400-600 GeV)", "DY_2l_M_50_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_600-800", "DY ll (m_{ll}>50 GeV, H_{T}: 600-800 GeV)", "DY_2l_M_50_HT_600-800", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_800-1200", "DY ll (m_{ll}>50 GeV, H_{T}: 800-1200 GeV)", "DY_2l_M_50_HT_800-1200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_1200-2500", "DY ll (m_{ll}>50 GeV, H_{T}: 1200-2500 GeV)", "DY_2l_M_50_HT_1200-2500", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_2500-inf", "DY ll (m_{ll}>50 GeV, H_{T}>2500 GeV)", "DY_2l_M_50_HT_2500-inf", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("TT2L2Nu", "t#bar{t} ll", "TT_2l2nu", -1, HistogramProperties((int) kOrange-3, 1, 2));
  sampleList.emplace_back("ZZ2L2Nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("WW2L2Nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  //sampleList.emplace_back("ST_t-channel_top_5f", "Single t", "ST_t-channel_top_5f", -1, HistogramProperties((int) kTeal-1, 1, 2));
  //sampleList.emplace_back("ST_t-channel_antitop_5f", "Single #bar{t}", "ST_t-channel_antitop_5f", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("ggZZ_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kAzure-2, 1, 2));
  sampleList.emplace_back("ggZZ_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", 125, HistogramProperties((int) kViolet, 7, 2));
  sampleList.emplace_back("VBF_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kCyan+2, 1, 2));
  sampleList.emplace_back("VBF_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1000_POWHEG", 125, HistogramProperties((int) kRed, 7, 2));

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  //PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  std::vector< std::vector<CutSpecs> > cutsets;
  // Nj==0, Nb==0
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 0, 0
  );
  // Nj==1, Nb==0
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 1, 1
  );
  // Nj==2, Nb==0
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 2, 2
  );

  for (size_t isample=0; isample<sampleList.size(); isample++){
    if (procsel>=0 && isample!=static_cast<size_t>(procsel)) continue;

    auto& sample = sampleList.at(isample);

    std::vector<TString> sampledirs;
    SampleHelpers::constructSamplesList(sample.path, theGlobalSyst, sampledirs);
    if (sampledirs.size()>1){
      MELAout << "Size > 1 not implemented yet!" << endl;
      continue;
    }
    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sampledirs.front()), EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());

    TString stroutput = Form("%s/%s%s", coutput_main.Data(), sample.name.data(), ".root");
    TString stroutput_txt = Form("%s/%s%s", coutput_main.Data(), sample.name.data(), ".txt");
    TFile* foutput = TFile::Open(stroutput, "recreate");
    MELAout.open(stroutput_txt.Data());

    // Get cross section
    sample_tree.bookBranch<float>("xsec", 0.f);

    // Configure handlers
    simEventHandler.bookBranches(&sample_tree);
    simEventHandler.wrapTree(&sample_tree);

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    //photonHandler.bookBranches(&sample_tree);
    //photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    foutput->cd();

    for (int ichannel=3; ichannel<nchannels; ichannel++){
      TString strChannel, strChannelLabel;
      getChannelTitleLabel(ichannel, strChannel, strChannelLabel);

      TDirectory* channeldir = foutput->mkdir(strChannel);
      channeldir->cd();

      for (auto const& cutset:cutsets){
        TString cutlabel, cuttitle;
        for (auto it_cut = cutset.cbegin(); it_cut != cutset.cend(); it_cut++){
          if (it_cut == cutset.cbegin()){
            cutlabel = it_cut->getLabel();
            cuttitle = it_cut->getTitle();
          }
          else{
            cutlabel = cutlabel + '|' + it_cut->getLabel();
            cuttitle = cuttitle + '_' + it_cut->getTitle();
          }
        }

        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s", strChannel.Data(), "puppimet_pTmiss_over_pTll_VS_pTll", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{ll} (GeV)", "p_{T}^{miss,PUPPI} / p_{T}^{ll}", "a.u.",
          100, 0., 1000.,
          40, 0., 4.,
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s", strChannel.Data(), "pfmet_pTmiss_over_pTll_VS_pTll", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{ll} (GeV)", "p_{T}^{miss,PF} / p_{T}^{ll}", "a.u.",
          100, 0., 1000.,
          40, 0., 4.,
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s", strChannel.Data(), "puppimet_pTmiss_over_pTlljets_VS_pTlljets", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{ll+jets} (GeV)", "p_{T}^{miss,PUPPI} / p_{T}^{ll+jets}", "a.u.",
          100, 0., 1000.,
          40, 0., 4.,
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s", strChannel.Data(), "pfmet_pTmiss_over_pTlljets_VS_pTlljets", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{ll+jets} (GeV)", "p_{T}^{miss,PF} / p_{T}^{ll+jets}", "a.u.",
          100, 0., 1000.,
          40, 0., 4.,
          channeldir
        );

        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_puppimet_pTmiss_pTll", cuttitle.Data(), "puppimet_pTmiss_ge_METthr"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PUPPI}>85 GeV"),
          "p_{T}^{ll} (GeV)", "|#Delta#phi(#vec{p}_{T}^{miss,PUPPI}, #vec{p}_{T}^{ll})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_pfmet_pTmiss_pTll", cuttitle.Data(), "pfmet_pTmiss_ge_METthr"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PF}>85 GeV"),
          "p_{T}^{ll} (GeV)", "|#Delta#phi(#vec{p}_{T}^{miss,PF}, #vec{p}_{T}^{ll})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_puppimet_pTmiss_pTlljets", cuttitle.Data(), "puppimet_pTmiss_ge_METthr_pTll_ge_35"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PUPPI}>85 GeV, p_{T}^{ll}>35 GeV"),
          "p_{T}^{ll+jets} (GeV)", "|#Delta#phi(#vec{p}_{T}^{miss,PUPPI}, #vec{p}_{T}^{ll+jets})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_pfmet_pTmiss_pTlljets", cuttitle.Data(), "pfmet_pTmiss_ge_METthr_pTll_ge_35"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PF}>85 GeV, p_{T}^{ll}>35 GeV"),
          "p_{T}^{ll+jets} (GeV)", "|#Delta#phi(#vec{p}_{T}^{miss,PF}, #vec{p}_{T}^{ll+jets})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );

        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_puppimet_pTmiss", cuttitle.Data(), "puppimet_pTmiss_ge_METthr"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PUPPI}>85 GeV"),
          "p_{T}^{ll} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PUPPI})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_pfmet_pTmiss", cuttitle.Data(), "pfmet_pTmiss_ge_METthr"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PF}>85 GeV"),
          "p_{T}^{ll} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PF})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_puppimet_pTmiss", cuttitle.Data(), "puppimet_pTmiss_ge_METthr_dPhi_pTlljets_puppimet_pTmiss_ge_2p6"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PUPPI}>85 GeV|#Delta#phi(#vec{p}_{T}^{miss,PUPPI}, #vec{p}_{T}^{ll+jets})>2.6"),
          "p_{T}^{ll} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PUPPI})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "dPhi_min_pTj_pfmet_pTmiss", cuttitle.Data(), "pfmet_pTmiss_ge_METthr_dPhi_pTlljets_pfmet_pTmiss_ge_2p6"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PF}>85 GeV|#Delta#phi(#vec{p}_{T}^{miss,PF}, #vec{p}_{T}^{ll+jets})>2.6"),
          "p_{T}^{ll} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PF})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );

        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "minimizedpTj_VS_dPhi_min_pTj_puppimet_pTmiss", cuttitle.Data(), "puppimet_pTmiss_ge_METthr_dPhi_pTlljets_puppimet_pTmiss_ge_2p6"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PUPPI}>85 GeV|#Delta#phi(#vec{p}_{T}^{miss,PUPPI}, #vec{p}_{T}^{ll+jets})>2.6"),
          "Minimized p_{T}^{j} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PUPPI})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
          channeldir
        );
        sample.hlist_2D.emplace_back(
          Form("h2D_%s_%s_%s_%s", strChannel.Data(), "minimizedpTj_VS_dPhi_min_pTj_pfmet_pTmiss", cuttitle.Data(), "pfmet_pTmiss_ge_METthr_dPhi_pTlljets_pfmet_pTmiss_ge_2p6"), Form("%s|%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data(), "p_{T}^{miss,PF}>85 GeV|#Delta#phi(#vec{p}_{T}^{miss,PF}, #vec{p}_{T}^{ll+jets})>2.6"),
          "Minimized p_{T}^{j} (GeV)", "Min. |#Delta#phi(#vec{p}_{T}^{j}, #vec{p}_{T}^{miss,PF})|", "a.u.",
          100, 0., 1000.,
          30, 0., +TMath::Pi(),
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
      if ((sample.name == "TT2L2Nu" || sample.name == "ZZ2L2Nu") && ev%10!=0) continue; // Take every 10 events in ttbar
      sample_tree.getSelectedEvent(ev);
      if (ev==0){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
      }

      simEventHandler.constructSimEvent(theGlobalSyst);

      genInfoHandler.constructGenInfo(theGlobalSyst);
      auto const& genInfo = genInfoHandler.getGenInfo();

      float genwgt = genInfo->getGenWeight(true);
      float puwgt = simEventHandler.getPileUpWeight();
      float wgt = genwgt * puwgt;
      sum_wgts += wgt;

      eventFilter.constructFilters();
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float me_wgt=1, cps_wgt=1;
      if (sample.name == "ggZZ_BSI"){ me_wgt = genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]; cps_wgt = genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      else if (sample.name == "ggZZ_Sig"){ me_wgt = genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]; cps_wgt = genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      else if (sample.name == "VBF_BSI"){ me_wgt = genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]; cps_wgt = genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      else if (sample.name == "VBF_Sig"){ me_wgt = genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]; cps_wgt = genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      wgt *= me_wgt*cps_wgt;

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      //photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, nullptr/*&photonHandler*/);

      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      //auto const& photons = photonHandler.getProducts();

      jetHandler.constructJetMET(&simEventHandler, theGlobalSyst, &muons, &electrons, nullptr/*&photons*/);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();
      ParticleObject::LorentzVector_t genmet;
      genmet = ParticleObject::PolarLorentzVector_t(genInfo->extras.genmet_met, 0, genInfo->extras.genmet_metPhi, 0);

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);
        }
      }

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      //MELAout << "Ndileptons: " << dileptons.size() << " | pTs = ";
      //for (auto const& dilepton:dileptons) MELAout << dilepton->pt() << " ";
      //MELAout << endl;

      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }

      float triggerWeight = eventFilter.getTriggerWeight(triggerCheckList);
      wgt *= triggerWeight;
      if (wgt==0.f) continue;

      if (theChosenDilepton && nTightDilep == 1){
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
          if (btagval>=0. && btagval<btagvalue_thr){
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
        float pZmiss_approx = theChosenDilepton->p4().Z();
        float etamiss_approx = theChosenDilepton->eta();
        float ml1j1 = p4_l1j1.M();
        float dEta_j1j2=-99; if (n_ak4jets_tight>=2) dEta_j1j2 = (p4_j.eta() - ak4jets_tight.at(1)->eta());
        float pTllj1 = (theChosenDilepton->p4() + p4_j).pt();
        float pTllj1j2 = (theChosenDilepton->p4() + p4_jj).pt();
        float pTlljets = (theChosenDilepton->p4() + p4_alljets).pt();

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

        float min_abs_dPhi_j_puppimet = TMath::Pi();
        float min_abs_dPhi_j_pfmet = TMath::Pi();
        float pTjet_min_abs_dPhi_j_puppimet = -1;
        float pTjet_min_abs_dPhi_j_pfmet = -1;
        for (AK4JetObject* jet:ak4jets_tight){
          float dphi_tmp;
          HelperFunctions::deltaPhi(float(jet->phi()), float(puppimet->phi()), dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          if (dphi_tmp<min_abs_dPhi_j_puppimet) pTjet_min_abs_dPhi_j_puppimet = jet->pt();
          min_abs_dPhi_j_puppimet = std::min(min_abs_dPhi_j_puppimet, dphi_tmp);

          HelperFunctions::deltaPhi(float(jet->phi()), float(pfmet->phi()), dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          if (dphi_tmp<min_abs_dPhi_j_pfmet) pTjet_min_abs_dPhi_j_pfmet = jet->pt();
          min_abs_dPhi_j_pfmet = std::min(min_abs_dPhi_j_pfmet, dphi_tmp);
        }
        float dPhi_pTll_puppimet_pTmiss = theChosenDilepton->deltaPhi(puppimet->phi());
        float dPhi_pTll_pfmet_pTmiss = theChosenDilepton->deltaPhi(pfmet->phi());
        float dPhi_pTllj_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_j).Phi()), float(puppimet->phi()), dPhi_pTllj_puppimet_pTmiss);
        float dPhi_pTlljj_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_jj).Phi()), float(puppimet->phi()), dPhi_pTlljj_puppimet_pTmiss);
        float dPhi_pTlljets_puppimet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(puppimet->phi()), dPhi_pTlljets_puppimet_pTmiss);
        float dPhi_pTlljets_pfmet_pTmiss; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfmet->phi()), dPhi_pTlljets_pfmet_pTmiss);
        float reco_puppimet_pTmiss = puppimet->pt();
        float reco_puppimet_pTmiss_significance = puppimet->extras.metSignificance;
        float resolution_puppimet_pTmiss = reco_puppimet_pTmiss/gen_pTmiss - 1.;
        // Compute ZZ-style masses
        float mTZZ_puppimet = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(reco_puppimet_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + puppimet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t puppimet_p4_ZZplusapprox; puppimet_p4_ZZplusapprox = ParticleObject::PolarLorentzVector_t(reco_puppimet_pTmiss, etamiss_approx, puppimet->phi(), PDGHelpers::Zmass);
        float mZZ_plus_puppimet = (puppimet_p4_ZZplusapprox + theChosenDilepton->p4()).M();

        float reco_pfmet_pTmiss = pfmet->pt();

        // Cuts
        bool pass_Nb_veto = n_ak4jets_tight_btagged==0;
        bool pass_pTl1 = pTl1>=25.;
        bool pass_pTl2 = pTl2>=20.;
        bool pass_mll = (doZZWW==0 && mll>=81. && mll<101.) || (doZZWW==1 && mll>=101.);
        bool pass_pTll = (pTll>=35.);
        bool pass_puppimet_thr = (reco_puppimet_pTmiss>=85. && reco_puppimet_pTmiss/pTlljets>=std::pow(85. / pTlljets, 1.5));
        bool pass_pfmet_thr = (reco_pfmet_pTmiss>=85. && reco_pfmet_pTmiss/pTlljets>=std::pow(85. / pTlljets, 1.5));
        bool pass_puppimet_dPhilljets_thr = (std::abs(dPhi_pTlljets_puppimet_pTmiss)>=2.6);
        bool pass_pfmet_dPhilljets_thr = (std::abs(dPhi_pTlljets_pfmet_pTmiss)>=2.6);
        bool pass_min_abs_dPhi_j_puppimet = (min_abs_dPhi_j_puppimet>=0.6);
        bool pass_min_abs_dPhi_j_pfmet = (min_abs_dPhi_j_pfmet>=0.6);

        if (!pass_Nb_veto) continue;

        // Fill histograms
        // Enclosed around braces to localize it_hist
        {
          auto it_hist = sample.hlist_2D.begin();
          for (int ichannel=3; ichannel<nchannels; ichannel++){
            bool isCorrectChannel = (
              ichannel==-1
              || (ichannel==0 && is_ee)
              || (ichannel==1 && is_mumu)
              || (ichannel==2 && is_emu)
              || (ichannel==3 && (is_ee || is_mumu))
              );

            for (auto const& cutset:cutsets){
              bool doFill=true;
              for (auto it_cut = cutset.cbegin(); it_cut != cutset.cend(); it_cut++){
                TString const& cutvar = it_cut->cutvar;
                float cutval=0;
                if (cutvar == "genmet") cutval = gen_pTmiss;
                else if (cutvar == "Nj") cutval = n_ak4jets_tight;
                else if (cutvar == "Nb") cutval = n_ak4jets_tight_btagged;
                else if (cutvar == "puppimet") cutval = reco_puppimet_pTmiss;
                doFill &= it_cut->testCut(cutval);
              }

              // pT ratios
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_min_abs_dPhi_j_puppimet) it_hist->hist.Fill(pTll, reco_puppimet_pTmiss/pTll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_min_abs_dPhi_j_pfmet) it_hist->hist.Fill(pTll, reco_pfmet_pTmiss/pTll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && (pass_pTl1 || n_ak4jets_tight==0) && pass_mll && pass_min_abs_dPhi_j_puppimet) it_hist->hist.Fill(pTlljets, reco_puppimet_pTmiss/pTlljets, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && (pass_pTl1 || n_ak4jets_tight==0) && pass_mll && pass_min_abs_dPhi_j_pfmet) it_hist->hist.Fill(pTlljets, reco_pfmet_pTmiss/pTlljets, wgt); it_hist++;

              // dPhi MET - pTll+x
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_puppimet_thr && pass_min_abs_dPhi_j_puppimet) it_hist->hist.Fill(pTll, std::abs(dPhi_pTll_puppimet_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_pfmet_thr && pass_min_abs_dPhi_j_pfmet) it_hist->hist.Fill(pTll, std::abs(dPhi_pTll_pfmet_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && (pass_pTl1 || n_ak4jets_tight==0) && pass_mll && pass_puppimet_thr && pass_min_abs_dPhi_j_puppimet) it_hist->hist.Fill(pTlljets, std::abs(dPhi_pTlljets_puppimet_pTmiss), wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && (pass_pTl1 || n_ak4jets_tight==0) && pass_mll && pass_pfmet_thr && pass_min_abs_dPhi_j_pfmet) it_hist->hist.Fill(pTlljets, std::abs(dPhi_pTlljets_pfmet_pTmiss), wgt); it_hist++;

              // dPhi pTj - MET
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_puppimet_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTll, min_abs_dPhi_j_puppimet, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_pfmet_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTll, min_abs_dPhi_j_pfmet, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_puppimet_thr && pass_puppimet_dPhilljets_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTll, min_abs_dPhi_j_puppimet, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_pfmet_thr && pass_pfmet_dPhilljets_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTll, min_abs_dPhi_j_pfmet, wgt); it_hist++;

              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_puppimet_thr && pass_puppimet_dPhilljets_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTjet_min_abs_dPhi_j_puppimet, min_abs_dPhi_j_puppimet, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_pTl1 && pass_pTl2 && pass_mll && pass_pfmet_thr && pass_pfmet_dPhilljets_thr && n_ak4jets_tight>0) it_hist->hist.Fill(pTjet_min_abs_dPhi_j_pfmet, min_abs_dPhi_j_pfmet, wgt); it_hist++;
            }
          } // End loop over channels
        } // End fill
      }

    } // End loop over events

    for (auto& hh:sample.hlist_2D){
      TH2F& hist = hh.hist;
      if (sum_wgts>0.) hist.Scale(xsec*1000.*59.7/sum_wgts);
    }
    sample.writeHistograms();

    MELAout.close();
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples
}


void plot_MET_pTll(int doZZWW, int iproc, bool doCondX, TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  //SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  //constexpr int nchannels = 4;

  TString const cinput_main = "output/" + strdate + "/MET_vs_pTll" + (doZZWW==0 ? "/ZZCuts" : "/WWCuts");
  TString coutput_main = cinput_main + "/Plots";
  gSystem->mkdir(coutput_main, true);

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

  std::vector<std::string> sampleList;
  //sampleList.emplace_back("DY_M10-50");
  //sampleList.emplace_back("DY_M50");
  sampleList.emplace_back("DY_2l_M_50_HT");
  sampleList.emplace_back("TT2L2Nu");
  sampleList.emplace_back("ZZ2L2Nu");
  sampleList.emplace_back("WW2L2Nu");
  sampleList.emplace_back("ggZZ_Sig");
  sampleList.emplace_back("VBF_Sig");

  std::vector< std::vector<TH2F*> > sample_hist_list;
  std::vector<std::string> samplePlottedList;

  bool firstFile = true;
  std::vector<TString> hnames;
  std::vector<TFile*> finputlist;
  int jproc=0;
  for (auto const& sample:sampleList){
    if (iproc<0 || jproc==iproc){
      TFile* finput = TFile::Open(cinput_main + '/' + sample.data() + ".root", "read");
      finputlist.push_back(finput);

      std::vector<TH2F*> hlist;
      HelperFunctions::extractHistogramsFromDirectory(finput, hlist);
      if (firstFile){
        for (TH2F* hh:hlist) hnames.push_back(hh->GetName());
      }
      if (hnames.size() == hlist.size()){
        sample_hist_list.push_back(hlist);
        samplePlottedList.push_back(sample);
      }

      if (firstFile) firstFile=false;
    }
    jproc++;
  }

  for (size_t iplot=0; iplot<hnames.size(); iplot++){
    int nbinsx=sample_hist_list[0][iplot]->GetNbinsX();
    int nbinsy=sample_hist_list[0][iplot]->GetNbinsY();
    std::vector<TH2F*> hlist;
    for (auto& it:sample_hist_list) hlist.push_back(it.at(iplot));
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& hist = hlist.at(is);
      auto& sample = samplePlottedList.at(is);

      
      TString strxtitle = hist->GetXaxis()->GetTitle();
      TString strytitle = hist->GetYaxis()->GetTitle();
      /*
      // FIXME: Temporary fix, fix pTll+jets label
      if (strytitle.Contains("p_{T}^{ll+jets}") && strxtitle == "p_{T}^{ll} (GeV)") hist->GetXaxis()->SetTitle("p_{T}^{ll+jets} (GeV)");
      strxtitle = hist->GetXaxis()->GetTitle();
      */

      if (doCondX) HelperFunctions::conditionalizeHistogram<TH2F>(hist, 0, nullptr, false);

      float zmin, zmax;
      HelperFunctions::findBinContentRange(hist, zmin, zmax, false, false, true);
      hist->GetZaxis()->SetRangeUser(0., zmax);
      //if (strxtitle.Contains("p_{T}")) hist->GetXaxis()->SetRangeUser(0., 200);
      HelperFunctions::wipeOverUnderFlows(hist, false, true);

      hist->GetXaxis()->CenterTitle();
      hist->GetYaxis()->CenterTitle();
      hist->GetZaxis()->CenterTitle();
      hist->GetZaxis()->SetNdivisions(505);
      hist->GetZaxis()->SetLabelFont(42);
      hist->GetZaxis()->SetLabelOffset(0.007);
      hist->GetZaxis()->SetLabelSize(0.0315);
      hist->GetZaxis()->SetTitleSize(0.04);
      hist->GetZaxis()->SetTitleOffset(1.1);
      hist->GetZaxis()->SetTitleFont(42);

      TString canvasname = hnames.at(iplot) + "_" + sample.data();
      HelperFunctions::replaceString(canvasname, "h2D_", "");
      if (doCondX) canvasname += "_CondX";
      if (!canvasname.Contains("eeORmumu")) continue;
      TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 1000, 800);
      canvas->cd();
      gStyle->SetOptStat(0);
      canvas->SetFillColor(0);
      canvas->SetBorderMode(0);
      canvas->SetBorderSize(2);
      canvas->SetTickx(1);
      canvas->SetTicky(1);
      canvas->SetLeftMargin(0.136);
      canvas->SetRightMargin(0.24);
      canvas->SetTopMargin(0.07);
      canvas->SetBottomMargin(0.13);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend* legend = new TLegend(
        0.55,
        0.90-0.10/4.*2.,
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

      TPaveText* pt = new TPaveText(0.12, 0.93, 0.68, 1, "brNDC");
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
      text = pt->AddText(0.83, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      TF1* f_metcut85 = nullptr;
      TF1* f_metcut85oblique = nullptr;
      if (strytitle.Contains(" / p_{T}")){
        f_metcut85 = new TF1("metcut85", "85/x", 1e-5, hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1));
        f_metcut85oblique = new TF1("metcut85oblique", "pow(85/x, 1.5)", 1e-5, hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1));
      }
      TLine* l_ptllcut35 = nullptr;
      if (strxtitle.Contains("p_{T}^{ll}")) l_ptllcut35 = new TLine(35, hist->GetYaxis()->GetBinLowEdge(1), 35, hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY()+1));
      
      {
        std::vector<TString> tmplist;
        TString htitle = hist->GetTitle();
        HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
        TString hlabel = tmplist.front();

        hist->SetTitle("");
        legend->AddEntry(hist, hlabel, "l");
      }

      hist->Draw("colz");
      if (f_metcut85){
        f_metcut85->SetLineColor(kRed);
        f_metcut85->SetLineWidth(2);
        f_metcut85->SetLineStyle(2);
        f_metcut85->Draw("csame");
      }
      if (f_metcut85oblique){
        f_metcut85oblique->SetLineColor(kViolet);
        f_metcut85oblique->SetLineWidth(2);
        f_metcut85oblique->SetLineStyle(2);
        f_metcut85oblique->Draw("csame");
      }
      if (l_ptllcut35){
        l_ptllcut35->SetLineColor(kBlue);
        l_ptllcut35->SetLineWidth(2);
        l_ptllcut35->SetLineStyle(2);
        l_ptllcut35->Draw("same");
      }

      //legend->Draw("same");
      pt->Draw();

      canvas->RedrawAxis();
      canvas->Modified();
      canvas->Update();
      canvas->SaveAs(coutput_main + "/" + canvasname + ".pdf");
      canvas->SaveAs(coutput_main + "/" + canvasname + ".png");

      delete pt;
      delete legend;
      delete l_ptllcut35;
      delete f_metcut85oblique;
      delete f_metcut85;
      canvas->Close();
    }
  }

  for (TFile*& finput:finputlist) finput->Close();
}
