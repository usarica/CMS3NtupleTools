#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


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
    if (mass<0){
      hist.SetFillColor(props.color);
      hist.SetFillStyle(3018);
    }

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
  if (cutvar.Contains("flag")) return cutvarlabel;
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
  TString strprintf = Form("%s%i%s", "%.", base10exponent, "f");
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

void getHistograms(int doZZWW, int procsel, bool doOldSelection=true, bool usePuppiMETForSelection=false, TString strdate="");
void getHistograms(int doZZWW, int procsel, bool doOldSelection=true, TString strdate=""){
  getHistograms(doZZWW, procsel, doOldSelection, !doOldSelection, strdate);
}

void getHistograms(int doZZWW, int procsel, bool doOldSelection, bool usePuppiMETForSelection, TString strdate){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  constexpr int nchannels = 4; // ichannel=0, 1, 2, 3 for ee, mumu, emu, ee+mumu

  TString const coutput_main = "output/" + strdate + (doZZWW==0 ? "/ZZCuts" : "/WWCuts") + (doOldSelection ? "/OldCuts" : "/NewCuts") + (usePuppiMETForSelection ? "/PUPPIMETCuts" : "/PFMETCuts");

  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", ((procsel>=11) ? "191212" : "hadoop:200101"));

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);
  float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  std::vector<std::string> triggerCheckList = doZZWW==0 ?
    OffshellTriggerHelpers::getHLTMenus(
      {
        OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle,
        OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle
      }
    )
    :
    OffshellTriggerHelpers::getHLTMenus(
      {
        OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kMuEle,
        OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kSingleEle
      }
    )
    ;

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("DY_2l_M_50_HT_70-100", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_70-100", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_100-200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_200-400", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_400-600", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_600-800", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_600-800", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_800-1200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_800-1200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_1200-2500", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_1200-2500", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_2500-inf", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_2500-inf", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("TT2L2Nu", "t#bar{t} ll", "TT_2l2nu", -1, HistogramProperties((int) kOrange-3, 1, 2));
  sampleList.emplace_back("ZZ2L2Nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("WW2L2Nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
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

  // ME block
  CMS3MELAHelpers::GMECBlock MEblock;
  std::vector<std::string> MElist{
    "Name:JJVBF_SIG_ghv1_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJVBF MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1",
    "Name:JJQCD_SIG_ghg2_1_JHUGen_JECNominal Alias:<Name> Process:HSMHiggs Production:JJQCD MatrixElement:JHUGen Cluster:J2JECNominal DefaultME:-1 Options:AddPConst=1"
  };
  CMS3MELAHelpers::setupMela(SampleHelpers::theDataYear, 125., TVar::ERROR); // Sets up MELA only once
  MEblock.buildMELABranches(MElist, false);
  Discriminant* DjjVBF = DiscriminantClasses::constructKDFromType(
    DiscriminantClasses::kDjjVBF,
    ANALYSISTREEPKGDATAPATH+"RecoMEConstants/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth"
  );
  CMS3MELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_Stable);
  MELAout << "getHistograms: DjjVBF is built!" << endl;

  std::vector< std::vector<CutSpecs> > cutsets;
  /*
  // Cuts on gen. MET
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(2);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    false, true, -1, 2
  );
  cutsets.back().emplace_back(
    "genmet", "p_{T}^{miss,true}",
    false, true, -1, 50
  );

  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(2);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    false, true, -1, 2
  );
  cutsets.back().emplace_back(
    "genmet", "p_{T}^{miss,true}",
    true, true, 50, 250
  );

  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(2);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    false, true, -1, 2
  );
  cutsets.back().emplace_back(
    "genmet", "p_{T}^{miss,true}",
    true, false, 250, -1
  );
  */

  // Nj<3
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 0, 3
  );
  // Nj==0
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 0, 0
  );
  // Nj==1
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 1, 1
  );
  // Nj==2
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(1);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 2, 2
  );
  // Nj==2, VBF
  cutsets.push_back(std::vector<CutSpecs>());
  cutsets.back().reserve(2);
  cutsets.back().emplace_back(
    "Nj", "N_{j}",
    true, true, 2, 2
  );
  cutsets.back().emplace_back(
    "VBF_flag", "VBF category",
    true, true, 1, 1
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

    for (int ichannel=(doZZWW==0 ? 3 : 0); ichannel<nchannels; ichannel++){
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

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "genmet_pTmiss", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{miss,true} (GeV)", "",
          80, 0., 800.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mll", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "m_{ll} (GeV)", "",
          40, 0., 200.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "pTll", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{ll} (GeV)", "",
          120, 0., 600.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "pTl1", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{l1}", "",
          50, 0., 500.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "pTl2", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "p_{T}^{l2}", "",
          40, 0., 400.,
          channeldir
        );

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "n_ak4jets_tight", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "N_{j}", "",
          5, 0, 5,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "n_ak4jets_tight_btagged", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "N_{b}", "",
          5, 0, 5,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mjj", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "m_{jj} (GeV)", "",
          50, 0., 1000.,
          channeldir
        );
        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "abs_dEta_jj", cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
          "|#eta_{jj}|", "",
          static_cast<int>(AK4JetSelectionHelpers::etaThr_skim_tight*2./0.2), 0., AK4JetSelectionHelpers::etaThr_skim_tight*2.,
          channeldir
        );

        for (unsigned char imet=0; imet<2; imet++){ // Loop over MET types
          TString strMETname = (imet==0 ? "pfmet" : "puppimet");
          TString strMETtitle = (imet==0 ? "PF" : "PUPPI");

          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("min_abs_dPhi_j_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("Min. |#phi_{j} - #phi_{miss}^{%s})|", strMETtitle.Data()), "",
            15, 0., +TMath::Pi(),
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("abs_dPhi_ll_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("|#phi_{ll} - #phi_{miss}^{%s})|", strMETtitle.Data()), "",
            15, 0., +TMath::Pi(),
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("abs_dPhi_lljets_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("|#phi_{ll+jets} - #phi_{miss}^{%s})|", strMETtitle.Data()), "",
            15, 0., +TMath::Pi(),
            channeldir
          );

          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("%s_pTmiss", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("p_{T}^{miss,%s} (GeV)", strMETtitle.Data()), "",
            32, 0., 800.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("%s_pTmiss_over_pTll", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("p_{T}^{miss,%s}/p_{T}^{ll}", strMETtitle.Data()), "",
            20, 0., 5.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("%s_pTmiss_over_pTlljets", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("p_{T}^{miss,%s}/p_{T}^{ll+jets}", strMETtitle.Data()), "",
            20, 0., 5.,
            channeldir
          );

          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("%s_pTmiss_over_genmet_pTmiss", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("p_{T}^{miss,%s}/p_{T}^{miss,true}", strMETtitle.Data()), "",
            50, 0., 5.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("abs_dPhi_%s_genmet", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("|#phi_{miss}^{%s} - #phi_{miss}^{true}|", strMETtitle.Data()), "",
            30, 0., +TMath::Pi(),
            channeldir
          );

          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("mTZZ_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("m_{T}^{ZZ,%s} (GeV)", strMETtitle.Data()), "",
            95, 150., 3000.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("mZZ_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("m^{ZZ,%s} (GeV)", strMETtitle.Data()), "",
            95, 150., 3000.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), Form("DjjVBF_%s", strMETname.Data()), cuttitle.Data()), Form("%s|%s|%s", sample.label.data(), strChannelLabel.Data(), cutlabel.Data()),
            Form("D_{jj}^{VBF,%s} (GeV)", strMETtitle.Data()), "",
            20, 0., 1.,
            channeldir
          );
        } // End loop over MET type

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
      if (sample.name == "TT2L2Nu" && ev%10!=0) continue; // Take every 10 events in ttbar
      if (sample.name == "ZZ2L2Nu" && ev%18!=0) continue; // Take every 10 events in ttbar
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
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, /*&photonHandler*/nullptr);

      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      //auto const& photons = photonHandler.getProducts();

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, /*&photons*/nullptr);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& pfmet = jetHandler.getPFMET();
      float genmet_pTmiss = genInfo->extras.genmet_met;
      float genmet_phi = genInfo->extras.genmet_metPhi;

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

      if (theChosenDilepton && nTightDilep == 1 && check_NoExtraLeptons(muons, electrons, theChosenDilepton, doOldSelection)){
        bool is_ee=false, is_mumu=false, is_emu=false;
        if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -121) is_ee=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -143) is_emu=true;
        else if (theChosenDilepton->daughter(0)->pdgId() * theChosenDilepton->daughter(1)->pdgId() == -169) is_mumu=true;

        if (is_ee) sum_ee += wgt;
        else if (is_mumu) sum_mumu += wgt;
        else if (is_emu) sum_emu += wgt;

        float mll = theChosenDilepton->m();
        float pTll = theChosenDilepton->pt();

        ParticleObject* leadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(0) : theChosenDilepton->daughter(1));
        ParticleObject* subleadingLepton = (theChosenDilepton->daughter(0)->pt() > theChosenDilepton->daughter(1)->pt() ? theChosenDilepton->daughter(1) : theChosenDilepton->daughter(0));
        ParticleObject::LorentzVector_t p4_l1 = leadingLepton->p4();
        ParticleObject::LorentzVector_t p4_l2 = subleadingLepton->p4();
        float pTl1 = leadingLepton->pt();
        float pTl2 = subleadingLepton->pt();

        size_t n_ak4jets_tight = ak4jets_tight.size();
        size_t n_ak4jets_tight_btagged = ak4jets_tight_btagged.size();

        ParticleObject::LorentzVector_t p4_j;
        ParticleObject::LorentzVector_t p4_jj;
        ParticleObject::LorentzVector_t p4_alljets;
        for (size_t ijet=0; ijet<n_ak4jets_tight; ijet++){
          AK4JetObject* jet = ak4jets_tight.at(ijet);
          if (ijet<1) p4_j = p4_j + jet->p4();
          if (ijet<2) p4_jj = p4_jj + jet->p4();
          p4_alljets = p4_alljets + jet->p4();
        }
        float mjj = p4_jj.M();
        float etamiss_approx = theChosenDilepton->eta();
        float abs_dEta_j1j2=-1; if (n_ak4jets_tight>=2) abs_dEta_j1j2 = std::abs(p4_j.eta() - ak4jets_tight.at(1)->eta());
        float pTlljets = (theChosenDilepton->p4() + p4_alljets).pt();

        float min_abs_dPhi_j_puppimet = TMath::Pi();
        float min_abs_dPhi_j_pfmet = TMath::Pi();
        for (AK4JetObject* jet:ak4jets_tight){
          float dphi_tmp;
          
          HelperFunctions::deltaPhi(float(jet->phi()), float(puppimet->phi()), dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          min_abs_dPhi_j_puppimet = std::min(min_abs_dPhi_j_puppimet, dphi_tmp);
          
          HelperFunctions::deltaPhi(float(jet->phi()), float(pfmet->phi()), dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          min_abs_dPhi_j_pfmet = std::min(min_abs_dPhi_j_pfmet, dphi_tmp);
        }
        float abs_dPhi_ll_puppimet = theChosenDilepton->deltaPhi(puppimet->phi()); abs_dPhi_ll_puppimet = std::abs(abs_dPhi_ll_puppimet);
        float abs_dPhi_lljets_puppimet; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(puppimet->phi()), abs_dPhi_lljets_puppimet); abs_dPhi_lljets_puppimet = std::abs(abs_dPhi_lljets_puppimet);
        float puppimet_pTmiss = puppimet->pt();
        float puppimet_pTmiss_over_genmet_pTmiss = puppimet_pTmiss/genmet_pTmiss;
        float abs_dPhi_puppimet_genmet; HelperFunctions::deltaPhi(float(puppimet->phi()), genmet_phi, abs_dPhi_puppimet_genmet); abs_dPhi_puppimet_genmet = std::abs(abs_dPhi_puppimet_genmet);
        float puppimet_pTmiss_over_pTll = puppimet_pTmiss/pTll;
        float puppimet_pTmiss_over_pTlljets = puppimet_pTmiss/pTlljets;
        float mTZZ_puppimet = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(puppimet_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + puppimet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t puppimet_p4_ZZapprox; puppimet_p4_ZZapprox = ParticleObject::PolarLorentzVector_t(puppimet_pTmiss, etamiss_approx, puppimet->phi(), PDGHelpers::Zmass);
        float mZZ_puppimet = (puppimet_p4_ZZapprox + theChosenDilepton->p4()).M();

        float abs_dPhi_ll_pfmet = theChosenDilepton->deltaPhi(pfmet->phi()); abs_dPhi_ll_pfmet = std::abs(abs_dPhi_ll_pfmet);
        float abs_dPhi_lljets_pfmet; HelperFunctions::deltaPhi(float((theChosenDilepton->p4() + p4_alljets).Phi()), float(pfmet->phi()), abs_dPhi_lljets_pfmet); abs_dPhi_lljets_pfmet = std::abs(abs_dPhi_lljets_pfmet);
        float pfmet_pTmiss = pfmet->pt();
        float pfmet_pTmiss_over_genmet_pTmiss = pfmet_pTmiss/genmet_pTmiss;
        float abs_dPhi_pfmet_genmet; HelperFunctions::deltaPhi(float(pfmet->phi()), genmet_phi, abs_dPhi_pfmet_genmet); abs_dPhi_pfmet_genmet = std::abs(abs_dPhi_pfmet_genmet);
        float pfmet_pTmiss_over_pTll = pfmet_pTmiss/pTll;
        float pfmet_pTmiss_over_pTlljets = pfmet_pTmiss/pTlljets;
        float mTZZ_pfmet = sqrt(pow(sqrt(pow(pTll, 2) + pow(mll, 2)) + sqrt(pow(pfmet_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfmet->p4()).Pt(), 2));
        ParticleObject::LorentzVector_t pfmet_p4_ZZapprox; pfmet_p4_ZZapprox = ParticleObject::PolarLorentzVector_t(pfmet_pTmiss, etamiss_approx, pfmet->phi(), PDGHelpers::Zmass);
        float mZZ_pfmet = (pfmet_p4_ZZapprox + theChosenDilepton->p4()).M();

        // Cuts
        bool pass_Nb_veto = check_Nb_veto(n_ak4jets_tight_btagged);
        bool pass_pTl1 = check_pTl1(pTl1, doOldSelection);
        bool pass_pTl2 = check_pTl2(pTl2, doOldSelection);
        bool pass_mll = check_ml1(mll, doOldSelection, doZZWW==0);
        bool pass_pTll = check_pTll(pTll, doOldSelection);

        bool pass_pTmiss_thr = check_pTmiss((!usePuppiMETForSelection ? pfmet_pTmiss : puppimet_pTmiss), doOldSelection, doZZWW==0);
        bool pass_pTmiss_over_pTlljets = check_pTmiss_over_pTlljets((!usePuppiMETForSelection ? pfmet_pTmiss : puppimet_pTmiss), pTlljets, doOldSelection, doZZWW==0);
        bool pass_dPhi_pTlljets_pTmiss = check_dPhi_pTlljets_pTmiss((!usePuppiMETForSelection ? abs_dPhi_lljets_pfmet : abs_dPhi_lljets_puppimet), doOldSelection);
        bool pass_dPhi_pTll_pTmiss = check_dPhi_pTll_pTmiss((!usePuppiMETForSelection ? abs_dPhi_ll_pfmet : abs_dPhi_ll_puppimet), doOldSelection);
        bool pass_min_abs_dPhi_j_met = check_min_abs_dPhi_pTj_pTmiss((!usePuppiMETForSelection ? min_abs_dPhi_j_pfmet : min_abs_dPhi_j_puppimet), doOldSelection);

        if (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met){
          if (is_ee) sum_ee_selected += wgt;
          else if (is_mumu) sum_mumu_selected += wgt;
          else if (is_emu) sum_emu_selected += wgt;
        }


        float DjjVBF_puppimet=-1;
        float DjjVBF_pfmet=-1;
        if (n_ak4jets_tight==2){
          std::unordered_map<std::string, float> ME_values;

          SimpleParticleCollection_t daughters_puppimet;
          daughters_puppimet.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(puppimet_p4_ZZapprox)));
          SimpleParticleCollection_t daughters_pfmet;
          daughters_pfmet.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(pfmet_p4_ZZapprox)));

          SimpleParticleCollection_t associated;
          associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(ak4jets_tight.at(0)->p4())));
          associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(ak4jets_tight.at(1)->p4())));

          CMS3MELAHelpers::melaHandle->setInputEvent(&daughters_puppimet, &associated, nullptr, false);
          MEblock.computeMELABranches();
          MEblock.pushMELABranches();
          MEblock.getBranchValues(ME_values); // Record the MEs into the EDProducer product
          CMS3MELAHelpers::melaHandle->resetInputEvent();
          DjjVBF->update({ ME_values["p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"], ME_values["p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"] }, mZZ_puppimet);
          DjjVBF_puppimet = *DjjVBF;

          CMS3MELAHelpers::melaHandle->setInputEvent(&daughters_pfmet, &associated, nullptr, false);
          MEblock.computeMELABranches();
          MEblock.pushMELABranches();
          MEblock.getBranchValues(ME_values); // Record the MEs into the EDProducer product
          CMS3MELAHelpers::melaHandle->resetInputEvent();
          DjjVBF->update({ ME_values["p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"], ME_values["p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"] }, mZZ_pfmet);
          DjjVBF_pfmet = *DjjVBF;
        }
        bool pass_VBF_selection = check_VBF_category((!usePuppiMETForSelection ? DjjVBF_pfmet : DjjVBF_puppimet), ak4jets_tight, theChosenDilepton, doOldSelection);

        // Fill histograms
        // Enclosed around braces to localize it_hist
        {
          auto it_hist = sample.hlist_1D.begin();
          for (int ichannel=(doZZWW==0 ? 3 : 0); ichannel<nchannels; ichannel++){
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
                if (cutvar == "Nj") cutval = n_ak4jets_tight;
                else if (cutvar == "puppimet") cutval = puppimet_pTmiss;
                else if (cutvar == "VBF_flag") cutval = pass_VBF_selection;
                doFill &= it_cut->testCut(cutval);
              }

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(genmet_pTmiss, wgt); it_hist++;

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(pTll, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(pTl1, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(pTl2, wgt); it_hist++;

              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(n_ak4jets_tight, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(n_ak4jets_tight_btagged, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met) && n_ak4jets_tight>=2) it_hist->hist.Fill(mjj, wgt); it_hist++;
              if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met) && n_ak4jets_tight>=2) it_hist->hist.Fill(abs_dEta_j1j2, wgt); it_hist++;

              for (unsigned char imet=0; imet<2; imet++){
                float const& min_abs_dPhi_j_met = (imet==0 ? min_abs_dPhi_j_pfmet : min_abs_dPhi_j_puppimet);
                float const& abs_dPhi_ll_met = (imet==0 ? abs_dPhi_ll_pfmet : abs_dPhi_ll_puppimet);
                float const& abs_dPhi_lljets_met = (imet==0 ? abs_dPhi_lljets_pfmet : abs_dPhi_lljets_puppimet);

                float const& met_pTmiss = (imet==0 ? pfmet_pTmiss : puppimet_pTmiss);
                float const& met_pTmiss_over_pTll = (imet==0 ? pfmet_pTmiss_over_pTll : puppimet_pTmiss_over_pTll);
                float const& met_pTmiss_over_pTlljets = (imet==0 ? pfmet_pTmiss_over_pTlljets : puppimet_pTmiss_over_pTlljets);

                float const& met_pTmiss_over_genmet_pTmiss = (imet==0 ? pfmet_pTmiss_over_genmet_pTmiss : puppimet_pTmiss_over_genmet_pTmiss);
                float const& abs_dPhi_met_genmet = (imet==0 ? abs_dPhi_pfmet_genmet : abs_dPhi_puppimet_genmet);

                float const& mTZZ_met = (imet==0 ? mTZZ_pfmet : mTZZ_puppimet);
                float const& mZZ_met = (imet==0 ? mZZ_pfmet : mZZ_puppimet);
                float const& DjjVBF_met = (imet==0 ? DjjVBF_pfmet : DjjVBF_puppimet);

                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss) && n_ak4jets_tight>=1) it_hist->hist.Fill(min_abs_dPhi_j_met, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(abs_dPhi_ll_met, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_min_abs_dPhi_j_met) && n_ak4jets_tight>=1) it_hist->hist.Fill(abs_dPhi_lljets_met, wgt); it_hist++;

                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(met_pTmiss, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(met_pTmiss_over_pTll, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met) && n_ak4jets_tight>=1) it_hist->hist.Fill(met_pTmiss_over_pTlljets, wgt); it_hist++;

                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(met_pTmiss_over_genmet_pTmiss, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(abs_dPhi_met_genmet, wgt); it_hist++;

                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(mTZZ_met, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met)) it_hist->hist.Fill(mZZ_met, wgt); it_hist++;
                if (isCorrectChannel && doFill && (pass_pTl1 && pass_pTl2 && pass_pTll && pass_Nb_veto && pass_mll && pass_pTmiss_thr && pass_pTmiss_over_pTlljets && pass_dPhi_pTll_pTmiss && pass_dPhi_pTlljets_pTmiss && pass_min_abs_dPhi_j_met) && n_ak4jets_tight>=2) it_hist->hist.Fill(DjjVBF_met, wgt); it_hist++;
              }
            }
          } // End loop over channels
        } // End fill
      }

    } // End loop over events

    constexpr bool normByXsec=true;
    for (auto& hh:sample.hlist_1D){
      TH1F& hist = hh.hist;
      if (normByXsec){
        if (sum_wgts>0.) hist.Scale(xsec*1000.*lumi/sum_wgts);
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

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples

  delete DjjVBF;
}

void makePlots(
  int doZZWW, bool doOldSelection, bool usePuppiMETForSelection, TString strdate="",
  bool useLogY=true, bool isStacked=true, int nfoci=0, int ifocus=0
){
  if (nfoci>0 && (ifocus>=nfoci || ifocus<0)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  //SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  //constexpr int nchannels = 4;

  TString const cinput_main = "output/" + strdate + (doZZWW==0 ? "/ZZCuts" : "/WWCuts") + (doOldSelection ? "/OldCuts" : "/NewCuts") + (usePuppiMETForSelection ? "/PUPPIMETCuts" : "/PFMETCuts");;
  TString coutput_main = cinput_main + "/Plots";
  if (isStacked) coutput_main += "/Stacked";
  else coutput_main += "/Unstacked";
  if (useLogY) coutput_main += "/LogY";
  else coutput_main += "/Linear";
  if (nfoci>0) coutput_main += Form("/Divide%i/Focus%i", nfoci, ifocus);
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
  //sampleList.emplace_back("DY_2l_M_10to50");
  //sampleList.emplace_back("DY_2l_M_50");
  sampleList.emplace_back("DY_2l_M_50_HT");
  sampleList.emplace_back("TT2L2Nu");
  sampleList.emplace_back("ZZ2L2Nu");
  sampleList.emplace_back("WW2L2Nu");
  sampleList.emplace_back("ggZZ_BSI");
  sampleList.emplace_back("VBF_BSI");
  sampleList.emplace_back("ggZZ_Sig");
  sampleList.emplace_back("VBF_Sig");

  std::vector<std::string> samplePlottedList;

  std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;

  bool firstFile = true;
  std::vector<TString> hnames;
  std::vector<TFile*> finputlist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample.data() + ".root";
    if (!HostHelpers::FileReadable(cinput)) continue;

    TFile* finput = TFile::Open(cinput, "read");
    finputlist.push_back(finput);

    std::vector<TH1F*> hlist;
    HelperFunctions::extractHistogramsFromDirectory(finput, hlist);
    if (firstFile){
      for (TH1F* hh:hlist) hnames.push_back(hh->GetName());
    }
    if (hnames.size() == hlist.size()){
      sample_hist_map[sample] = hlist;
      samplePlottedList.push_back(sample);
    }

    if (firstFile) firstFile=false;
  }

  for (size_t iplot=0; iplot<hnames.size(); iplot++){
    int nbins=sample_hist_map[samplePlottedList.front()].at(0)->GetNbinsX();
    int binstart=1;
    int binend=nbins;
    if (nfoci>0){
      int bininc = nbins / nfoci;
      binstart = (ifocus * bininc) + 1;
      binend = (ifocus==nfoci-1 ? nbins : binstart + bininc);
    }
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);
      hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinLowEdge(binstart), hist->GetXaxis()->GetBinUpEdge(binend));
    }
    if (isStacked){
      for (size_t is=0; is<samplePlottedList.size(); is++){
        auto& sample = samplePlottedList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        bool isSignal = sample.find("Sig")!=std::string::npos;
        for (size_t js=is+1; js<samplePlottedList.size(); js++){
          auto& sample_j = samplePlottedList.at(js);
          TH1F* hist_j = sample_hist_map[sample_j].at(iplot);
          bool isSignal_j = sample_j.find("Sig")!=std::string::npos;
          if (isSignal_j!=isSignal) continue;
          hist->Add(hist_j);
        }
      }
    }
    else{
      for (size_t is=0; is<samplePlottedList.size(); is++){
        auto& sample = samplePlottedList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        double inthist = hist->Integral(binstart, binend);
        hist->Scale(1. / inthist);
      }
    }

    size_t nplottables=0;
    double ymin = 0;
    double ymax = -1;
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (doZZWW==0 && TString(hist->GetName()).Contains("emu")) continue;
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
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
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
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (doZZWW==0 && TString(hist->GetName()).Contains("emu")) continue;

      std::vector<TString> tmplist;
      TString htitle = hist->GetTitle();
      HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
      TString hlabel = tmplist.front();

      hist->SetTitle("");
      if (sample.find("Sig")!=string::npos){
        hist->SetFillStyle(3354);
        hist->SetFillColor(hist->GetLineColor());
      }
      else{
        hist->SetLineColor(kBlack);
        hist->SetFillStyle(1001);
      }
      if (sample=="DY_2l_M_50_HT") hlabel = "DY (HT-binned)";
      legend->AddEntry(hist, hlabel, "f");

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
