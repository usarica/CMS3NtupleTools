#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace reco;


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
  srcDir(srcDir_)
{
  TDirectory* curdir = gDirectory;
  if (srcDir) srcDir->cd();
  hist = TH1F(name, title, nbins, xlow, xhigh);
  curdir->cd();
  hist.GetXaxis()->SetTitle(xlabel);
  hist.GetYaxis()->SetTitle(ylabel);

  cout << "Created histogram " << hist.GetName() << " [" << hist.GetTitle() << "]" << endl;
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

  void setupHistogram(TH1F& hist) const;
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
void HistogramProperties::setupHistogram(TH1F& hist) const{
  hist.SetMarkerColor(color);
  hist.SetLineColor(color);
  hist.SetLineStyle(dashtype);
  hist.SetLineWidth(width);
}


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

  void setupHistogram(TH1F& hist) const;
  void setupHistogram(TH1F* const& hist) const{ if (hist) setupHistogram(*hist); }
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

    setupHistogram(hist);

    if (mass<0){
      hist.SetFillColor(props.color);
      hist.SetFillStyle(3018);
    }
  }
}
void SampleSpecs::setupHistogram(TH1F& hist) const{
  hist.Sumw2();

  props.setupHistogram(hist);

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

void getCutSets(std::vector< std::vector<CutSpecs> >& cutsets){
  for (int ijet=0; ijet<=3; ijet++){
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_loose_extraLeptons_flag", "No extra loose l",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_stdveto_lepId_flag", "Veto stop, l",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_stdveto_anyId_flag", "Veto stop, any",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_stdveto_anyId_or_tightLeptons_flag", "Veto stop, any, or tight leptons",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_veto_lepId_flag", "New veto, l",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_veto_lepId_hadId_flag", "New veto, l+ch",
      true, true, 0, 0
    );
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, (ijet==3 ? false : true), (ijet==3 ? 3 : ijet), (ijet==3 ? -1 : ijet)
    );
    cutsets.back().emplace_back(
      "event_has_isotracks_veto_lepId_hadId_or_tightLeptons_flag", "New veto, l+ch+tight leptons",
      true, true, 0, 0
    );
  }
}


void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts="hist",
  bool adjustYLow=false,
  float factorYHigh=-1
);


using namespace SystematicsHelpers;
void getTrees(
  int procsel, int ichunk, int nchunks, TString strdate,
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (procsel<0) return;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main =
    "output/TTBarClosure/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
  gSystem->mkdir(coutput_main, true);

  SampleHelpers::configure("2018", "hadoop:200326");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btag_loose_thr = BtagHelpers::getBtagWP(false);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

   auto triggerPropsCheckList_OSDF = TriggerHelpers::getHLTMenuProperties(
    {
      TriggerHelpers::kMuEle,
      TriggerHelpers::kSingleEle, TriggerHelpers::kSingleMu
    }
  );
  auto triggerPropsCheckList_OSSF = TriggerHelpers::getHLTMenuProperties(
    {
      TriggerHelpers::kDoubleEle, TriggerHelpers::kSingleEle,
      TriggerHelpers::kDoubleMu, TriggerHelpers::kSingleMu
    }
  );

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 0
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_70-100", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_70-100", -1, HistogramProperties((int) kCyan, 1, 2)); // 5
  sampleList.emplace_back("DY_2l_M_50_HT_100-200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_200-400", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_400-600", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_600-800", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_600-800", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_800-1200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_800-1200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_1200-2500", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_1200-2500", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_2500-inf", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_2500-inf", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_10to50_ext", "DY ll (m_{ll}: [10, 50] GeV)", "DY_2l_M_10to50_ext", -1, HistogramProperties((int) kCyan, 1, 2)); // 13
  sampleList.emplace_back("DY_2l_M_50_ext", "DY ll (m_{ll}>50 GeV)", "DY_2l_M_50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("2018A", "Observed (2018A)", "Run2018A", -1, HistogramProperties((int) kCyan, 1, 2)); // 15
  sampleList.emplace_back("2018B", "Observed (2018B)", "Run2018B", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("2018C", "Observed (2018C)", "Run2018C", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("2018D", "Observed (2018D)", "Run2018D", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 19
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 23
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M125_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 27
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M125_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 29
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 31
  sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M125_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 35
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_WW2L2Nu_M125_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));

  if ((unsigned int) procsel>=sampleList.size()) return;
  auto const& sample = sampleList.at(procsel);
  TString const& strSampleSet = sample.path;
  if (
    sample.name.find("Sig")!=std::string::npos
    ||
    sample.name.find("Bkg")!=std::string::npos
    ||
    sample.name.find("BSI")!=std::string::npos
    ) SampleHelpers::configure("2018", "hadoop:200313");
  else{
    SampleHelpers::setInputDirectory("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/DileptonSkims/2018");
  }

  // Get sample specifications
  bool const isData = SampleHelpers::checkSampleIsData(strSampleSet);
  if (isData && theGlobalSyst!=sNominal) return;
  if (isData && nchunks>0) return;

  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;

  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  TString stroutput = Form("%s/%s", coutput_main.Data(), sample.name.data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
  stroutput += Form("_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  TTree* tout = new TTree("SkimTree", "");
#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
  // Event variables
  BRANCH_COMMAND(float, event_wgt);
  BRANCH_COMMAND(float, event_wgt_undo_eff);
  BRANCH_COMMAND(float, pTmiss);
  BRANCH_COMMAND(float, phimiss);
  BRANCH_COMMAND(float, mTZZ);
  BRANCH_COMMAND(float, mZZ);
  BRANCH_COMMAND(float, mjj);
  BRANCH_COMMAND(unsigned int, event_nvtxs_good);
  BRANCH_COMMAND(unsigned int, event_Nphotons);
  BRANCH_COMMAND(unsigned int, event_Njets);
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_loose);
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_medium);
  BRANCH_COMMAND(unsigned int, event_Njets20);
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_loose);
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_medium);
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers);
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers);
  BRANCH_COMMAND(bool, is_ee);
  BRANCH_COMMAND(bool, is_mumu);
  BRANCH_COMMAND(bool, is_emu);
  BRANCH_COMMAND(bool, has_electrons_inHEM1516);
  BRANCH_COMMAND(bool, has_photons_inHEM1516);
  BRANCH_COMMAND(bool, has_ak4jets_inHEM1516);
  BRANCH_COMMAND(bool, has_ak8jets_inHEM1516);
  // LL
  BRANCH_COMMAND(float, pt_ll);
  BRANCH_COMMAND(float, eta_ll);
  BRANCH_COMMAND(float, phi_ll);
  BRANCH_COMMAND(float, mass_ll);
  // L1, L2
  BRANCH_COMMAND(int, id_l1);
  BRANCH_COMMAND(float, pt_l1);
  BRANCH_COMMAND(float, eta_l1);
  BRANCH_COMMAND(float, eff_l1);
  BRANCH_COMMAND(float, efferr_l1);
  BRANCH_COMMAND(int, id_l2);
  BRANCH_COMMAND(float, pt_l2);
  BRANCH_COMMAND(float, eta_l2);
  BRANCH_COMMAND(float, eff_l2);
  BRANCH_COMMAND(float, efferr_l2);
  // J1
  BRANCH_COMMAND(float, pt_j1);
  BRANCH_COMMAND(float, eta_j1);
  BRANCH_COMMAND(float, phi_j1);
  BRANCH_COMMAND(float, mass_j1);
  // J2
  BRANCH_COMMAND(float, pt_j2);
  BRANCH_COMMAND(float, eta_j2);
  BRANCH_COMMAND(float, phi_j2);
  BRANCH_COMMAND(float, mass_j2);
#undef BRANCH_COMMAND

  bool isFirstInputFile=true;
  for (auto const& sname:sampledirs){
    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);

    // Set data tracking options
    eventFilter.setTrackDataEvents(isData);
    eventFilter.setCheckUniqueDataEvent(isData && !isFirstInputFile);

    const int nEntries = sample_tree.getSelectedNEvents();

    float sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      TString firstFile = SampleHelpers::getDatasetFirstFileName(sname);
      TFile* fFirstFile = TFile::Open(firstFile, "read");
      TH2F* hCounters = (TH2F*) fFirstFile->Get("cms3ntuple/Counters");
      if (hCounters) sum_wgts = hCounters->GetBinContent(1, 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp));
      else{
        MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
        for (int ev=0; ev<nEntries; ev++){
          HelperFunctions::progressbar(ev, nEntries);
          sample_tree.getSelectedEvent(ev);

          genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
          auto const& genInfo = genInfoHandler.getGenInfo();
          float genwgt = genInfo->getGenWeight(true);

          simEventHandler.constructSimEvent(theGlobalSyst);
          float puwgt = simEventHandler.getPileUpWeight();
          sum_wgts += genwgt * puwgt;
        }
      }
      fFirstFile->Close();
    }

    // Configure handlers
    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    isotrackHandler.bookBranches(&sample_tree);
    isotrackHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    foutput->cd();

    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events from " << sample_tree.sampleIdentifier << ", starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    size_t n_evts_acc=0;
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      //if (ev>100) break;

      if (!isData && firstEvent){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        xsec *= 1000.;
      }
      if (firstEvent) firstEvent=false;

      event_wgt = xsec * (isData ? 1.f : lumi) / sum_wgts;

      if (!isData){
        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();
        event_wgt *= genInfo->getGenWeight(true);

        if (event_wgt==0.f) continue;

        if (
          sample.name == "ggZZ_2l2nu_BSI"
          ||
          sample.name == "ggWW_2l2nu_BSI"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (
          sample.name == "ggZZ_2l2nu_Sig"
          ||
          sample.name == "ggWW_2l2nu_Sig"
          ||
          sample.name == "ggZZ_2l2nu_Sig_Onshell"
          ||
          sample.name == "ggWW_2l2nu_Sig_Onshell"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (
          sample.name == "ggZZ_2l2nu_BSI_g4"
          ||
          sample.name == "ggWW_2l2nu_BSI_g4"
          ){ event_wgt *= (genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz4_1_MCFM"]*2.5502 + genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*(-2.55052+std::pow(2.55052, 2)) + genInfo->extras.LHE_ME_weights["p_Gen_GG_BKG_MCFM"]*(-2.5502+1.))*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (
          sample.name == "ggZZ_2l2nu_Sig_g4"
          ||
          sample.name == "ggWW_2l2nu_Sig_g4"
          ||
          sample.name == "ggZZ_2l2nu_Sig_Onshell_g4"
          ||
          sample.name == "ggWW_2l2nu_Sig_Onshell_g4"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*std::pow(2.55052, 2)*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (
          sample.name == "VBFZZ_2l2nu_BSI"
          ||
          sample.name == "VBFWW_2l2nu_BSI"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
        else if (
          sample.name == "VBFZZ_2l2nu_Sig"
          ||
          sample.name == "VBFWW_2l2nu_Sig"
          ||
          sample.name == "VBFZZ_2l2nu_Sig_Onshell"
          ||
          sample.name == "VBFWW_2l2nu_Sig_Onshell"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
        else if (
          sample.name == "VBFZZ_2l2nu_BSI_g4"
          ||
          sample.name == "VBFWW_2l2nu_BSI_g4"
          ){
          event_wgt *= (
            genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*std::pow(2.55052, 2)
            + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv4_1_MCFM"]*(-std::pow(2.55052, 2)+std::pow(2.55052, 4))
            + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BKG_MCFM"]*(-std::pow(2.55052, 2)+1.)
            )*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"];
        }
        else if (
          sample.name == "VBFZZ_2l2nu_Sig_g4"
          ||
          sample.name == "VBFWW_2l2nu_Sig_g4"
          ||
          sample.name == "VBFZZ_2l2nu_Sig_Onshell_g4"
          ||
          sample.name == "VBFWW_2l2nu_Sig_Onshell_g4"
          ){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv4_1_MCFM"]*std::pow(2.55052, 4)*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }

        if (event_wgt==0.f) continue;

        simEventHandler.constructSimEvent(theGlobalSyst);
        event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();

        if (event_wgt==0.f) continue;
      }
      //MELAout << "Pass line " << __LINE__ << endl;

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      auto const& photons = photonHandler.getProducts();
      event_Nphotons = 0;
      for (auto const& part:photons){
        if (ParticleSelectionHelpers::isTightParticle(part)) event_Nphotons++;
      }

      isotrackHandler.constructIsotracks(&muons, &electrons);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      auto pfmet_p4 = pfmet->p4(true, true, true);
      pTmiss = pfmet_p4.Pt();
      phimiss = pfmet_p4.Phi();

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      event_wgt_OSSF_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_OSSF, &muons, &electrons, nullptr, nullptr, nullptr, nullptr);
      event_wgt_OSDF_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_OSDF, &muons, &electrons, nullptr, nullptr, nullptr, nullptr);
      if ((event_wgt_OSSF_triggers + event_wgt_OSDF_triggers) == 0.f) continue;

      has_electrons_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, &electrons, nullptr, nullptr, nullptr);
      has_photons_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, nullptr, &photons, nullptr, nullptr);
      has_ak4jets_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, nullptr);
      has_ak8jets_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, nullptr, &ak8jets);

      mjj = -1;
      pt_j1 = -1;
      eta_j1 = 0;
      phi_j1 = 0;
      mass_j1 = 0;
      pt_j2 = -1;
      eta_j2 = 0;
      phi_j2 = 0;
      mass_j2 = 0;
      event_Njets_btagged_loose = event_Njets_btagged_medium = event_Njets20_btagged_loose = event_Njets20_btagged_medium = event_Njets = event_Njets20 = 0;
      ParticleObject::LorentzVector_t ak4jets_sump4;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets_btagged_loose++;
          if (jet->getBtagValue()>=btag_medium_thr) event_Njets_btagged_medium++;
          ak4jets_sump4 += jet->p4();
        }
        if (
          jet->testSelectionBit(AK4JetSelectionHelpers::kTightId)
          &&
          (!applyPUIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kPUJetId))
          &&
          (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
          &&
          jet->pt()>=20.f && fabs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
          ){
          event_Njets20++;
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets20_btagged_loose++;
          if (jet->getBtagValue()>=btag_medium_thr) event_Njets20_btagged_medium++;
        }
      }
      event_Njets = ak4jets_tight.size();
      if (event_Njets>0){
        auto const& jet1 = ak4jets_tight.front();
        pt_j1 = jet1->pt();
        eta_j1 = jet1->eta();
        phi_j1 = jet1->phi();
        mass_j1 = jet1->m();
        if (event_Njets>1){
          auto const& jet2 = ak4jets_tight.at(1);
          pt_j2 = jet2->pt();
          eta_j2 = jet2->eta();
          phi_j2 = jet2->phi();
          mass_j2 = jet2->m();

          mjj = (jet1->p4() + jet2->p4()).M();
        }
      }

      bool pass_dPhi_jet_met = true;
      for (auto const& jet:ak4jets_tight){
        float dphi_tmp;
        HelperFunctions::deltaPhi(float(jet->phi()), phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
        if (dphi_tmp<=0.25){ pass_dPhi_jet_met=false; break; }
      }
      if (!pass_dPhi_jet_met) continue;

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
      if (!theChosenDilepton || nTightDilep>1) continue;

      bool pass_dPhi_ll_met = true;
      {
        float abs_dPhi_ll_pfmet = theChosenDilepton->deltaPhi(phimiss); abs_dPhi_ll_pfmet = std::abs(abs_dPhi_ll_pfmet);
        pass_dPhi_ll_met = abs_dPhi_ll_pfmet>1.0;
      }
      if (!pass_dPhi_ll_met) continue;

      bool pass_dPhi_lljets_met = true;
      {
        float abs_dPhi_lljets_pfmet;
        HelperFunctions::deltaPhi(float((theChosenDilepton->p4()+ak4jets_sump4).Phi()), phimiss, abs_dPhi_lljets_pfmet);
        abs_dPhi_lljets_pfmet = std::abs(abs_dPhi_lljets_pfmet);
        pass_dPhi_lljets_met = abs_dPhi_lljets_pfmet>2.5;
      }
      if (!pass_dPhi_lljets_met) continue;

      std::vector<ParticleObject*> dileptonDaughters = theChosenDilepton->getDaughters();
      ParticleObjectHelpers::sortByGreaterPt(dileptonDaughters);
      ParticleObject* leadingLepton = dileptonDaughters.front();
      ParticleObject* subleadingLepton = dileptonDaughters.back();

      is_ee = is_mumu = is_emu=false;
      if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -121) is_ee=true;
      else if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -143) is_emu=true;
      else if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -169) is_mumu=true;
      if (
        ((is_ee || is_mumu) && event_wgt_OSSF_triggers==0.f)
        ||
        (is_emu && event_wgt_OSDF_triggers==0.f)
        ) continue;

      sample_tree.getVal("vtxs_nvtxs_good", event_nvtxs_good);

      mass_ll = theChosenDilepton->m();
      pt_ll = theChosenDilepton->pt();
      eta_ll = theChosenDilepton->eta();
      phi_ll = theChosenDilepton->phi();

      mTZZ = sqrt(pow(sqrt(pow(pt_ll, 2) + pow(mass_ll, 2)) + sqrt(pow(pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfmet_p4).Pt(), 2));

      float etamiss_approx = theChosenDilepton->eta();
      ParticleObject::LorentzVector_t pfmet_p4_approx; pfmet_p4_approx = ParticleObject::PolarLorentzVector_t(pTmiss, etamiss_approx, phimiss, PDGHelpers::Zmass);
      ParticleObject::LorentzVector_t pfmet_ZZ_p4_approx = pfmet_p4_approx + theChosenDilepton->p4();
      mZZ = pfmet_ZZ_p4_approx.M();

      id_l1 = leadingLepton->pdgId();
      pt_l1 = leadingLepton->pt();
      eta_l1 = leadingLepton->eta();
      id_l2 = subleadingLepton->pdgId();
      pt_l2 = subleadingLepton->pt();
      eta_l2 = subleadingLepton->eta();

      //if (std::abs(id_l1)==11) electronSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (ElectronObject const*) leadingLepton, isData, false);
      //else muonSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (MuonObject const*) leadingLepton, isData, false);
      //if (std::abs(id_l2)==11) electronSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (ElectronObject const*) subleadingLepton, isData, false);
      //else muonSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (MuonObject const*) subleadingLepton, isData, false);

      tout->Fill();
      n_evts_acc++;
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc << " / " << (ev_end - ev_start) << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;
  } // End loop over samples list

  // Write output
  foutput->WriteTObject(tout);
  foutput->Close();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void count(
  int flagSB=0,
  bool useMediumBtagForOutSel=false,
  TString strdate="",
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200313");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString strSideband = "SBupdn";
  if (flagSB==1) strSideband="onlySBup";
  else if (flagSB==-1) strSideband="onlySBdn";

  TString const cinput_main =
    "output/TTBarClosure/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
  TString coutput_main =
    "output/TTBarClosure/SkimTrees/" + strdate + "/Plots"
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + (useMediumBtagForOutSel ? "MediumBtagOut" : "LooseBtagOut") + "_" + strSideband + "/Counts";
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = Form("%s/Integrals", coutput_main.Data());
  if (std::abs(flagSB)==1) stroutput_txt += Form("_%s", strSideband.Data());
  stroutput_txt += ".txt";
  MELAout.open(stroutput_txt.Data());

  TDirectory* curdir = gDirectory;

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

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(float, phimiss) \
  BRANCH_COMMAND(float, mTZZ) \
  BRANCH_COMMAND(float, mZZ) \
  BRANCH_COMMAND(float, mjj) \
  BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
  BRANCH_COMMAND(unsigned int, event_Nphotons) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_Njets20) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_medium) \
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers) \
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers) \
  BRANCH_COMMAND(bool, is_ee) \
  BRANCH_COMMAND(bool, is_mumu) \
  BRANCH_COMMAND(bool, is_emu) \
  BRANCH_COMMAND(bool, has_electrons_inHEM1516) \
  BRANCH_COMMAND(bool, has_photons_inHEM1516) \
  BRANCH_COMMAND(bool, has_ak4jets_inHEM1516) \
  BRANCH_COMMAND(bool, has_ak8jets_inHEM1516) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(int, id_l1) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(int, id_l2) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, pt_j1) \
  BRANCH_COMMAND(float, eta_j1) \
  BRANCH_COMMAND(float, phi_j1) \
  BRANCH_COMMAND(float, mass_j1) \
  BRANCH_COMMAND(float, pt_j2) \
  BRANCH_COMMAND(float, eta_j2) \
  BRANCH_COMMAND(float, phi_j2) \
  BRANCH_COMMAND(float, mass_j2)
#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("Data_2018", "Observed", "2018", -1, HistogramProperties((int) kBlack, 1, 2));
  //
  sampleList.emplace_back("Total_noHiggs_noZpeak", "Total (no Higgs, no Z)", "Total_noHiggs_noZpeak", -1, HistogramProperties((int) kBlue, 1, 2));
  sampleList.emplace_back("Total_noHiggs", "Total (no Higgs)", "Total_noHiggs", -1, HistogramProperties((int) kBlue, 1, 2));
  sampleList.emplace_back("Total_noZpeak_g1", "Total (SM Higgs, no Z)", "Total_noZpeak_g1", -1, HistogramProperties((int) kGreen+2, 2, 2));
  sampleList.emplace_back("Total_g1", "Total (SM Higgs)", "Total", -1, HistogramProperties((int) kGreen+2, 2, 2));
  sampleList.emplace_back("Total_noZpeak_g4", "Total (PS Higgs, no Z)", "Total_noZpeak_g4", -1, HistogramProperties((int) kViolet, 7, 2));
  sampleList.emplace_back("Total_g4", "Total (PS Higgs)", "Total_g4", -1, HistogramProperties((int) kViolet, 7, 2));
  //
  /*
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kOrange-3, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  /*
  sampleList.emplace_back("DY_2l_M_50_HT", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI", -1, HistogramProperties((int) kViolet, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI", -1, HistogramProperties((int) kBlue, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kViolet, 7, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kBlue, 7, 2));
  //
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell", "gg#rightarrowWW sig. on-shell (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_Onshell", -1, HistogramProperties((int) kGreen+2, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell", "EW WW+jj sig. on-shell (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_Onshell", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell_g4", "gg#rightarrowWW sig. on-shell (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_Onshell_g4", -1, HistogramProperties((int) kGreen+2, 7, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell_g4", "EW WW+jj sig. on-shell (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_Onshell_g4", -1, HistogramProperties((int) kYellow-3, 7, 2));

  curdir->cd();

  std::vector<TString> const channel_labels={ "ee", "mumu", "emu" };
  std::vector<TString> const channel_hlabels={ "ee", "#mu#mu", "e#mu" };
  std::vector<TString> const in_out_labels={ "in", "out" };
  constexpr int nbins_Nj=4;

  constexpr float inf_pTmiss = 40;
  constexpr float sup_pTmiss = 130;
  constexpr int nbins_pTmiss = (sup_pTmiss-inf_pTmiss)/5.;


  std::vector<SampleSpecs> sampleProcessedList; sampleProcessedList.reserve(sampleList.size());
  std::vector<SampleSpecs> sampleExtraList; sampleExtraList.reserve(sampleList.size());
  std::vector<TChain*> treelist;
  for (auto const& sample:sampleList){
    if (sample.name.find("Total")!=std::string::npos){
      sampleExtraList.push_back(sample);
      continue;
    }
    TChain* tin = new TChain("SkimTree");
    TString cinput = cinput_main + '/' + sample.path.data() + Form("*%s.root", SystematicsHelpers::getSystName(theGlobalSyst).data());
    if (
      sample.name.find("Sig")!=std::string::npos
      ||
      sample.name.find("BSI")!=std::string::npos
      ||
      sample.name.find("Bkg")!=std::string::npos
      ) cinput = cinput_main + '/' + sample.path.data() + Form("_%s.root", SystematicsHelpers::getSystName(theGlobalSyst).data());
    tin->Add(cinput);
    if (tin->GetEntries()==0){
      delete tin;
      continue;
    }
#define BRANCH_COMMAND(type, name) tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    treelist.push_back(tin);
    MELAout << "Extracted input file " << cinput << endl;

    sampleProcessedList.push_back(sample);
  }
  const size_t nsamples = sampleProcessedList.size();
  curdir->cd();

  TFile* foutput = TFile::Open(coutput_main + "/histograms_all.root", "recreate");

  std::vector<TH1F> h_pTmiss[3][nbins_Nj][2]; // Channel (ee, mumu, emu), Njets, in/out-like selection
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      for (size_t io=0; io<in_out_labels.size(); io++) h_pTmiss[ic][ijet][io].reserve(sampleProcessedList.size() + sampleExtraList.size());
    }
  }
  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << isample << " (" << sampleProcessedList.at(isample).name << ")" << endl;

    TH1F* hfill_pTmiss[3][nbins_Nj][2];
    TString strhnamecore = Form("%s_pTmiss", sampleProcessedList.at(isample).name.data());
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);
        for (size_t io=0; io<in_out_labels.size(); io++){
          TString hnametmp = Form("%s_%s_%s_%s", strhnamecore.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data());
          MELAout << "\t- Creating pTmiss histogram " << hnametmp << endl;
          h_pTmiss[ic][ijet][io].emplace_back(
            hnametmp,
            sampleProcessedList.at(isample).label.data(),
            nbins_pTmiss, inf_pTmiss, sup_pTmiss
          );
          hfill_pTmiss[ic][ijet][io] = &(h_pTmiss[ic][ijet][io].back());
          hfill_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} threshold (GeV)");
          hfill_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("Events (p_{T}^{miss} > thr.)");
          sampleProcessedList.at(isample).setupHistogram(*(hfill_pTmiss[ic][ijet][io]));
        }
      }
    }

    float xsecScale = 1;
    if (sampleProcessedList.at(isample).name.find("VBFWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.064714*0.00847/0.00814; // BR

    if (sampleProcessedList.at(isample).name.find("ggZZ")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.541; // BR scale
    if (sampleProcessedList.at(isample).name.find("ggWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.616; // BR scale
    if (sampleProcessedList.at(isample).name.find("VBFZZ")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.312; // BR scale
    if (sampleProcessedList.at(isample).name.find("VBFWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.616*0.407/0.412; // BR scale

    auto const& tin = treelist.at(isample);
    int nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (has_ak4jets_inHEM1516 || has_ak8jets_inHEM1516) continue;

      if (pt_ll<55.f || pt_l1<25.f || pt_l2<25.f) continue;
      int index_channel = 0*is_ee + 1*is_mumu + 2*is_emu;

      if (sampleProcessedList.at(isample).name.find("ggWW_2l2nu_BSI")!=std::string::npos && std::abs(event_wgt)>150.) continue;
      if (sampleProcessedList.at(isample).name.find("ggWW_2l2nu_BSI_g4")!=std::string::npos && std::abs(event_wgt)>200.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_BSI")!=std::string::npos && std::abs(event_wgt)>500.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_BSI_g4")!=std::string::npos && std::abs(event_wgt)>30000.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_Sig_Onshell_g4")!=std::string::npos && std::abs(event_wgt)>200.) continue;
      event_wgt *= xsecScale;
      if (index_channel<2) event_wgt *= event_wgt_OSSF_triggers;
      else event_wgt *= event_wgt_OSDF_triggers;

      if (flagSB!=1 && mass_ll<40.f) continue;
      if (flagSB!=-1 && mass_ll>=200.f) continue;

      int index_Njets = std::min(nbins_Nj-1, static_cast<int>(event_Njets));
      int index_INorOUT=-1;
      if (std::abs(mass_ll-91.2f)<15.f){
        index_INorOUT = 0;
        if (pTmiss<125.f) continue;
        if (event_Njets_btagged_loose!=0) continue;
      }
      else{
        index_INorOUT = 1;
        if (pTmiss<40.f) continue;
        if ((event_Njets>0 && (useMediumBtagForOutSel ? event_Njets_btagged_medium : event_Njets_btagged_loose)==0) || (event_Njets==0 && (useMediumBtagForOutSel ? event_Njets20_btagged_medium : event_Njets20_btagged_loose)!=0)) continue;
      }

      //MELAout << "Attempting to fill (" << index_channel << ", " << index_Njets << ", " << index_INorOUT << endl;
      hfill_pTmiss[index_channel][index_Njets][index_INorOUT]->Fill(pTmiss, event_wgt);
    }
  }

  for (auto const& sample:sampleExtraList){
    bool const isNoZ = sample.name.find("noZpeak")!=std::string::npos;
    bool const isNoHiggs = sample.name.find("noHiggs")!=std::string::npos;
    bool const isSMHiggs = sample.name.find("_g1")!=std::string::npos;
    bool const isPSHiggs = sample.name.find("_g4")!=std::string::npos;

    TString strhnamecore = Form("%s_pTmiss", sample.name.data());
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);
        for (size_t io=0; io<in_out_labels.size(); io++){
          h_pTmiss[ic][ijet][io].emplace_back(
            Form("%s_%s_%s_%s", strhnamecore.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()),
            sample.label.data(),
            nbins_pTmiss, inf_pTmiss, sup_pTmiss
          );
          TH1F& hfill = h_pTmiss[ic][ijet][io].back();
          sample.setupHistogram(hfill);

          for (size_t isample=0; isample<nsamples; isample++){
            auto const& tmpsample = sampleProcessedList.at(isample);
            if (tmpsample.name.find("Data")!=std::string::npos) continue;

            bool const hasZ = !(
              tmpsample.name.find("DY")==std::string::npos
              &&
              tmpsample.name.find("ZZ_2l2nu")==std::string::npos
              &&
              tmpsample.name.find("WZ_3l")==std::string::npos
              );
            bool const hasHiggs = !(
              tmpsample.name.find("BSI")==std::string::npos
              &&
              tmpsample.name.find("Sig")==std::string::npos
              );
            bool const hasOnshellHiggs = tmpsample.name.find("Onshell")!=std::string::npos;
            bool const hasHiggsBSI = hasHiggs && tmpsample.name.find("BSI")!=std::string::npos;
            bool const hasHiggsSig = hasHiggs && tmpsample.name.find("Sig")!=std::string::npos;
            bool const hasSMHiggs = (hasHiggsBSI || hasHiggsSig) && tmpsample.name.find("_g4")==std::string::npos;
            bool const hasPSHiggs = (hasHiggsBSI || hasHiggsSig) && tmpsample.name.find("_g4")!=std::string::npos;

            if (isNoZ && hasZ) continue;
            if (isNoHiggs && hasHiggs) continue;
            if (isSMHiggs && hasHiggsBSI && !hasSMHiggs) continue;
            if (isSMHiggs && hasHiggsSig && !hasSMHiggs) continue;
            if (isSMHiggs && hasHiggsSig && !hasOnshellHiggs) continue;
            if (isPSHiggs && hasHiggsBSI && !hasPSHiggs) continue;
            if (isPSHiggs && hasHiggsSig && !hasPSHiggs) continue;
            if (isPSHiggs && hasHiggsSig && !hasOnshellHiggs) continue;

            TH1F& htmp = h_pTmiss[ic][ijet][io].at(isample);
            hfill.Add(&htmp, 1);

            //MELAout << "\t- " << hfill.GetName() << " += " <<  htmp.GetName() << endl;
          }
        }
      }
    }
  }

  // Convert and save raw pTmiss histograms
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      for (size_t io=0; io<in_out_labels.size(); io++){
        for (auto& htmp:h_pTmiss[ic][ijet][io]){
          int const nbinsx = htmp.GetNbinsX()+1;
          for (int ix=0; ix<=nbinsx+1; ix++){
            double yield[2]={ 0 }; // Value, error
            yield[0] = HelperFunctions::getHistogramIntegralAndError(&htmp, ix, nbinsx, false, &(yield[1]));
            htmp.SetBinContent(ix, yield[0]);
            htmp.SetBinError(ix, yield[1]);
          }
          foutput->WriteTObject(&htmp);
        }
      }
    }
  }


  std::vector<TH1F*> h_alpha_predicted[2][nbins_Nj]; // Channel (ee, mumu), Njets
  std::vector<TH1F*> h_yield_predicted[2][nbins_Nj]; // Channel (ee, mumu), Njets
  std::vector<TH1F*> h_yield_expected[2][nbins_Nj]; // Channel (ee, mumu), Njets
  std::vector<TH1F*> h_yield_diff_norm[2][nbins_Nj]; // Channel (ee, mumu), Njets
  std::vector<TH1F*> h_yield_predicted_fracerr[2][nbins_Nj]; // Channel (ee, mumu), Njets
  for (size_t ic=0; ic<2; ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      size_t iproc=0;
      for (auto const& htmp:h_pTmiss[ic][ijet][0]){
        TString hname_alpha_predicted = htmp.GetName(); HelperFunctions::replaceString(hname_alpha_predicted, "pTmiss", "alpha_predicted");
        h_alpha_predicted[ic][ijet].push_back((TH1F*) htmp.Clone(hname_alpha_predicted)); h_alpha_predicted[ic][ijet].back()->Reset("ICESM"); h_alpha_predicted[ic][ijet].back()->GetYaxis()->SetTitle("#alpha (p_{T}^{miss} > thr.)");
        TString hname_yield_predicted = htmp.GetName(); HelperFunctions::replaceString(hname_yield_predicted, "pTmiss", "yield_predicted");
        h_yield_predicted[ic][ijet].push_back((TH1F*) htmp.Clone(hname_yield_predicted)); h_yield_predicted[ic][ijet].back()->Reset("ICESM"); h_yield_predicted[ic][ijet].back()->GetYaxis()->SetTitle("Events (p_{T}^{miss} > thr.)");
        TString hname_yield_expected = htmp.GetName(); HelperFunctions::replaceString(hname_yield_expected, "pTmiss", "yield_expected");
        h_yield_expected[ic][ijet].push_back((TH1F*) htmp.Clone(hname_yield_expected)); h_yield_expected[ic][ijet].back()->Reset("ICESM"); h_yield_expected[ic][ijet].back()->GetYaxis()->SetTitle("Events (p_{T}^{miss} > thr.)");

        TString hname_yield_diff_norm = htmp.GetName(); HelperFunctions::replaceString(hname_yield_diff_norm, "pTmiss", "yield_diff_norm");
        h_yield_diff_norm[ic][ijet].push_back((TH1F*) htmp.Clone(hname_yield_diff_norm)); h_yield_diff_norm[ic][ijet].back()->Reset("ICESM"); h_yield_diff_norm[ic][ijet].back()->GetYaxis()->SetTitle("(Exp. - pred.)/#sqrt{pred} (p_{T}^{miss} > thr.)");
        TString hname_yield_predicted_fracerr = htmp.GetName(); HelperFunctions::replaceString(hname_yield_predicted_fracerr, "pTmiss", "yield_predicted_fracerr");
        h_yield_predicted_fracerr[ic][ijet].push_back((TH1F*) htmp.Clone(hname_yield_predicted_fracerr)); h_yield_predicted_fracerr[ic][ijet].back()->Reset("ICESM"); h_yield_predicted_fracerr[ic][ijet].back()->GetYaxis()->SetTitle("Calc. frac. stat. err. (p_{T}^{miss} > thr.)");

        int const nbinsx = htmp.GetNbinsX()+1;
        int const bin125 = htmp.GetXaxis()->FindBin(125.0001);

        // The actual yield we expect
        double yield_channel_in[2]={ htmp.GetBinContent(bin125), htmp.GetBinError(bin125) }; // Value, error

        // 'in'-like yield from emu
        double yield_emu_in[2]={ h_pTmiss[2][ijet][0].at(iproc).GetBinContent(bin125), h_pTmiss[2][ijet][0].at(iproc).GetBinError(bin125) }; // Value, error

        for (int ix=0; ix<=nbinsx+1; ix++){
          h_yield_expected[ic][ijet].back()->SetBinContent(ix, yield_channel_in[0]);
          h_yield_expected[ic][ijet].back()->SetBinError(ix, yield_channel_in[1]);

          // The yield we estimate from emu and ee/mumu 'out'-like selections
          double yield_channel_out[2]={ h_pTmiss[ic][ijet][1].at(iproc).GetBinContent(ix), h_pTmiss[ic][ijet][1].at(iproc).GetBinError(ix) }; // Value, error
          double yield_emu_out[2]={ h_pTmiss[2][ijet][1].at(iproc).GetBinContent(ix), h_pTmiss[2][ijet][1].at(iproc).GetBinError(ix) }; // Value, error

          double alpha_predicted[3]={ 0 };
          alpha_predicted[0] = yield_channel_out[0] / yield_emu_out[0];
          alpha_predicted[1] = alpha_predicted[0] * std::sqrt(std::pow(yield_channel_out[1]/yield_channel_out[0], 2) + std::pow(yield_emu_out[1]/yield_emu_out[0], 2));
          alpha_predicted[2] = alpha_predicted[0] * std::sqrt(1./std::abs(yield_channel_out[0]) + 1./std::abs(yield_emu_out[0]));
          h_alpha_predicted[ic][ijet].back()->SetBinContent(ix, alpha_predicted[0]);
          h_alpha_predicted[ic][ijet].back()->SetBinError(ix, alpha_predicted[1]);

          double yield_predicted[3]={ 0 };
          yield_predicted[0] = yield_emu_in[0] * alpha_predicted[0];
          yield_predicted[1] = yield_predicted[0] * std::sqrt(std::pow(yield_emu_in[1]/yield_emu_in[0], 2) + std::pow(alpha_predicted[1]/alpha_predicted[0], 2));
          yield_predicted[2] = yield_predicted[0] * std::sqrt(1./std::abs(yield_emu_in[0]) + std::pow(alpha_predicted[2]/alpha_predicted[0], 2));
          h_yield_predicted[ic][ijet].back()->SetBinContent(ix, yield_predicted[0]);
          h_yield_predicted[ic][ijet].back()->SetBinError(ix, yield_predicted[1]);

          //h_yield_diff_norm[ic][ijet].back()->SetBinContent(ix, (yield_predicted[0] - yield_channel_in[0])/std::sqrt(std::pow(yield_predicted[1], 2) + std::pow(yield_channel_in[1], 2)));
          //h_yield_diff_norm[ic][ijet].back()->SetBinContent(ix, (yield_predicted[0] - yield_channel_in[0])/std::sqrt(std::pow(yield_predicted[1], 2)));
          h_yield_diff_norm[ic][ijet].back()->SetBinContent(ix, (yield_predicted[0] - yield_channel_in[0])/std::sqrt(std::abs(yield_predicted[0])));

          h_yield_predicted_fracerr[ic][ijet].back()->SetBinContent(ix, yield_predicted[2]/std::abs(yield_predicted[0]));
        }
        iproc++;
      }
    }
  }

  // Make plots of alpha and predictions
  for (size_t ic=0; ic<2; ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);

      TString hlabel_Nj;
      if (ijet==nbins_Nj-1) hlabel_Nj = Form("N_{j}#geq%i", ijet);
      else hlabel_Nj = Form("N_{j}=%i", ijet);

      TString selectionLabels = channel_hlabels.at(ic) + ", " + hlabel_Nj;
      if (useMediumBtagForOutSel) selectionLabels = selectionLabels + "|'Out'-like medium b";
      else selectionLabels = selectionLabels + "|'Out'-like loose b";

      std::vector<TH1F*> hlist;
      std::vector<TString> hlabels;

      // alpha
      {
        for (auto const& htmp:h_alpha_predicted[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("Total_noHiggs_noZpeak") || hname.Contains("Total_noZpeak") || hname.Contains("Data")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        for (auto const& htmp:h_alpha_predicted[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("TT_2l2nu")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        makePlot(
          coutput_main, lumi,
          Form("c_alpha_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
          hlist, hlabels,
          selectionLabels,
          "e1p", true, 1.2
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }

      // Plot predicted yields
      {
        for (auto const& htmp:h_yield_predicted[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("Total_noHiggs_noZpeak") || hname.Contains("Total_noZpeak") || hname.Contains("Data")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        for (auto const& htmp:h_yield_predicted[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("TT_2l2nu")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        for (auto const& htmp:h_yield_expected[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("Total_noHiggs_noZpeak") || hname.Contains("Total_noZpeak") || hname.Contains("TT_2l2nu")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back("");
            for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++) hlist.back()->SetBinError(ix, 0);
          }
        }
        makePlot(
          coutput_main, lumi,
          Form("c_yield_predicted_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
          hlist, hlabels,
          selectionLabels,
          "e1p", true
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }

      // Plot prediction diff_norm
      {
        for (auto const& htmp:h_yield_diff_norm[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("Total_noHiggs_noZpeak") || hname.Contains("Total_noZpeak")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        for (auto const& htmp:h_yield_diff_norm[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("TT_2l2nu")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        makePlot(
          coutput_main, lumi,
          Form("c_yield_diff_norm_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
          hlist, hlabels,
          selectionLabels,
          "e1p", true
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }

      // Plot fractional errors
      {
        for (auto const& htmp:h_yield_predicted_fracerr[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("Total_noHiggs_noZpeak") || hname.Contains("Total_noZpeak")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        for (auto const& htmp:h_yield_predicted_fracerr[ic][ijet]){
          TString hname = htmp->GetName();
          if (hname.Contains("TT_2l2nu")){
            hlist.push_back((TH1F*) htmp->Clone(Form("copy_%s", htmp->GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp->GetTitle());
          }
        }
        makePlot(
          coutput_main, lumi,
          Form("c_yield_predicted_fracerr_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
          hlist, hlabels,
          selectionLabels,
          "e1p", true, 1.2
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }

    }
  }

  // Delete prediction histograms and close the output file
  for (size_t ic=0; ic<2; ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      for (auto& htmp:h_alpha_predicted[ic][ijet]){ foutput->WriteTObject(htmp); delete htmp; }
      for (auto& htmp:h_yield_predicted[ic][ijet]){ foutput->WriteTObject(htmp); delete htmp; }
      for (auto& htmp:h_yield_expected[ic][ijet]){ foutput->WriteTObject(htmp); delete htmp; }
      for (auto& htmp:h_yield_diff_norm[ic][ijet]){ foutput->WriteTObject(htmp); delete htmp; }
      for (auto& htmp:h_yield_predicted_fracerr[ic][ijet]){ foutput->WriteTObject(htmp); delete htmp; }
    }
  }
  foutput->Close();

  for (auto& tin:treelist) delete tin;
  MELAout.close();
#undef BRANCH_COMMANDS
}


void plot_Step0(
  int flagSB=0,
  bool useMediumBtagForOutSel=false,
  TString strdate="",
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200313");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString strSideband = "SBupdn";
  if (flagSB==1) strSideband="onlySBup";
  else if (flagSB==-1) strSideband="onlySBdn";

  TString const cinput_main =
    "output/TTBarClosure/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
  TString coutput_main =
    "output/TTBarClosure/SkimTrees/" + strdate + "/Plots"
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + (useMediumBtagForOutSel ? "MediumBtagOut" : "LooseBtagOut") + "_" + strSideband + "/Step0";
  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;

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

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(float, phimiss) \
  BRANCH_COMMAND(float, mTZZ) \
  BRANCH_COMMAND(float, mZZ) \
  BRANCH_COMMAND(float, mjj) \
  BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
  BRANCH_COMMAND(unsigned int, event_Nphotons) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_Njets20) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged_medium) \
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers) \
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers) \
  BRANCH_COMMAND(bool, is_ee) \
  BRANCH_COMMAND(bool, is_mumu) \
  BRANCH_COMMAND(bool, is_emu) \
  BRANCH_COMMAND(bool, has_electrons_inHEM1516) \
  BRANCH_COMMAND(bool, has_photons_inHEM1516) \
  BRANCH_COMMAND(bool, has_ak4jets_inHEM1516) \
  BRANCH_COMMAND(bool, has_ak8jets_inHEM1516) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(int, id_l1) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(int, id_l2) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, pt_j1) \
  BRANCH_COMMAND(float, eta_j1) \
  BRANCH_COMMAND(float, phi_j1) \
  BRANCH_COMMAND(float, mass_j1) \
  BRANCH_COMMAND(float, pt_j2) \
  BRANCH_COMMAND(float, eta_j2) \
  BRANCH_COMMAND(float, phi_j2) \
  BRANCH_COMMAND(float, mass_j2)
#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("Data_2018", "Observed", "2018", -1, HistogramProperties((int) kBlack, 1, 2));
  //
  sampleList.emplace_back("ggWW_2l2nu_BSI", "Total, f_{ai}=0, #Gamma_{H}=#Gamma_{H}^{SM}", "ggWW_2l2nu_BSI", -1, HistogramProperties((int) kViolet, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI", "Total, f_{ai}=0, #Gamma_{H}=#Gamma_{H}^{SM}", "VBFWW_2l2nu_BSI", -1, HistogramProperties((int) kBlue, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "Total, f_{a3}=1, #Gamma_{H}=#Gamma_{H}^{SM}", "ggWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kViolet, 7, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "Total, f_{a3}=1, #Gamma_{H}=#Gamma_{H}^{SM}", "VBFWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kBlue, 7, 2));
  //
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell", "On-shell H, f_{ai}=0, #Gamma_{H}=#Gamma_{H}^{SM}", "ggWW_2l2nu_Sig_Onshell", -1, HistogramProperties((int) kGreen+2, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell", "On-shell H, f_{ai}=0, #Gamma_{H}=#Gamma_{H}^{SM}", "VBFWW_2l2nu_Sig_Onshell", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_Onshell_g4", "On-shell H, f_{ai}=1, #Gamma_{H}=#Gamma_{H}^{SM}", "ggWW_2l2nu_Sig_Onshell_g4", -1, HistogramProperties((int) kGreen+2, 7, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_Onshell_g4", "On-shell H, f_{ai}=1, #Gamma_{H}=#Gamma_{H}^{SM}", "VBFWW_2l2nu_Sig_Onshell_g4", -1, HistogramProperties((int) kYellow-3, 7, 2));
  //
  /*
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kOrange-3, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  /*
  sampleList.emplace_back("DY_2l", "DY ll", "DY_2l_M_10to50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */

  curdir->cd();

  std::vector<TString> const channel_labels={ "ee", "mumu", "emu" };
  std::vector<TString> const channel_hlabels={ "ee", "#mu#mu", "e#mu" };
  std::vector<TString> const in_out_labels={ "in", "out" };
  constexpr int nbins_Nj=4;

  constexpr float inf_pTmiss = 40;
  constexpr float sup_pTmiss = 355;
  constexpr int nbins_pTmiss = (sup_pTmiss-inf_pTmiss)/15.;

  constexpr float inf_mll = 40;
  constexpr float sup_mll = 200;
  constexpr int nbins_mll = (sup_mll-inf_mll)/10.;

  constexpr float inf_pTll = 55;
  constexpr float sup_pTll = 255;
  constexpr int nbins_pTll = (sup_pTll-inf_pTll)/20.;

  constexpr float inf_mTZZ = 150;
  constexpr float sup_mTZZ = 500;
  constexpr int nbins_mTZZ = (sup_mTZZ-inf_mTZZ)/25.;

  constexpr float inf_nvtxs_good = 0;
  constexpr float sup_nvtxs_good = 100;
  constexpr int nbins_nvtxs_good = (sup_nvtxs_good-inf_nvtxs_good)/4.;

  std::vector<SampleSpecs> sampleProcessedList; sampleProcessedList.reserve(sampleList.size());
  std::vector<SampleSpecs> sampleExtraList; sampleExtraList.reserve(sampleList.size());
  std::vector<TChain*> treelist;
  for (auto const& sample:sampleList){
    if (sample.name.find("Total")!=std::string::npos){
      sampleExtraList.push_back(sample);
      continue;
    }
    TChain* tin = new TChain("SkimTree");
    TString cinput = cinput_main + '/' + sample.path.data() + Form("*%s.root", SystematicsHelpers::getSystName(theGlobalSyst).data());
    if (
      sample.name.find("Sig")!=std::string::npos
      ||
      sample.name.find("BSI")!=std::string::npos
      ||
      sample.name.find("Bkg")!=std::string::npos
      ) cinput = cinput_main + '/' + sample.path.data() + Form("_%s.root", SystematicsHelpers::getSystName(theGlobalSyst).data());
    tin->Add(cinput);
    if (tin->GetEntries()==0){
      delete tin;
      continue;
    }
#define BRANCH_COMMAND(type, name) tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    treelist.push_back(tin);
    MELAout << "Extracted input file " << cinput << endl;

    sampleProcessedList.push_back(sample);
  }
  const size_t nsamples = sampleProcessedList.size();
  curdir->cd();

  std::vector<TH1F> h_pTmiss[3][nbins_Nj]; // Channel (ee, mumu, emu), Njets
  std::vector<TH1F> h_mll[3][nbins_Nj]; // Channel (ee, mumu, emu), Njets
  std::vector<TH1F> h_pTll[3][nbins_Nj]; // Channel (ee, mumu, emu), Njets
  std::vector<TH1F> h_mTZZ[3][nbins_Nj]; // Channel (ee, mumu, emu), Njets
  std::vector<TH1F> h_nvtxs_good[3][nbins_Nj]; // Channel (ee, mumu, emu), Njets
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      h_pTmiss[ic][ijet].reserve(sampleProcessedList.size() + sampleExtraList.size());
      h_mll[ic][ijet].reserve(sampleProcessedList.size() + sampleExtraList.size());
      h_pTll[ic][ijet].reserve(sampleProcessedList.size() + sampleExtraList.size());
      h_mTZZ[ic][ijet].reserve(sampleProcessedList.size() + sampleExtraList.size());
      h_nvtxs_good[ic][ijet].reserve(sampleProcessedList.size() + sampleExtraList.size());
    }
  }
  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << isample << " (" << sampleProcessedList.at(isample).name << ")" << endl;

    TH1F* hfill_pTmiss[3][nbins_Nj];
    TH1F* hfill_mll[3][nbins_Nj];
    TH1F* hfill_pTll[3][nbins_Nj];
    TH1F* hfill_mTZZ[3][nbins_Nj];
    TH1F* hfill_nvtxs_good[3][nbins_Nj];
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);

        h_pTmiss[ic][ijet].emplace_back(
          Form("%s_pTmiss_%s_%s", sampleProcessedList.at(isample).name.data(), channel_labels.at(ic).Data(), label_Nj.Data()),
          sampleProcessedList.at(isample).label.data(),
          nbins_pTmiss, inf_pTmiss, sup_pTmiss
        ); hfill_pTmiss[ic][ijet] = &(h_pTmiss[ic][ijet].back()); hfill_pTmiss[ic][ijet]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); hfill_pTmiss[ic][ijet]->GetYaxis()->SetTitle("Events"); sampleProcessedList.at(isample).setupHistogram(*(hfill_pTmiss[ic][ijet]));
        h_mll[ic][ijet].emplace_back(
          Form("%s_mll_%s_%s", sampleProcessedList.at(isample).name.data(), channel_labels.at(ic).Data(), label_Nj.Data()),
          sampleProcessedList.at(isample).label.data(),
          nbins_mll, inf_mll, sup_mll
        ); hfill_mll[ic][ijet] = &(h_mll[ic][ijet].back()); hfill_mll[ic][ijet]->GetXaxis()->SetTitle("m_{ll} (GeV)"); hfill_mll[ic][ijet]->GetYaxis()->SetTitle("Events"); sampleProcessedList.at(isample).setupHistogram(*(hfill_mll[ic][ijet]));
        h_pTll[ic][ijet].emplace_back(
          Form("%s_pTll_%s_%s", sampleProcessedList.at(isample).name.data(), channel_labels.at(ic).Data(), label_Nj.Data()),
          sampleProcessedList.at(isample).label.data(),
          nbins_pTll, inf_pTll, sup_pTll
        ); hfill_pTll[ic][ijet] = &(h_pTll[ic][ijet].back()); hfill_pTll[ic][ijet]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); hfill_pTll[ic][ijet]->GetYaxis()->SetTitle("Events"); sampleProcessedList.at(isample).setupHistogram(*(hfill_pTll[ic][ijet]));
        h_mTZZ[ic][ijet].emplace_back(
          Form("%s_mTZZ_%s_%s", sampleProcessedList.at(isample).name.data(), channel_labels.at(ic).Data(), label_Nj.Data()),
          sampleProcessedList.at(isample).label.data(),
          nbins_mTZZ, inf_mTZZ, sup_mTZZ
        ); hfill_mTZZ[ic][ijet] = &(h_mTZZ[ic][ijet].back()); hfill_mTZZ[ic][ijet]->GetXaxis()->SetTitle("m_{T}^{ZZ} (GeV)"); hfill_mTZZ[ic][ijet]->GetYaxis()->SetTitle("Events"); sampleProcessedList.at(isample).setupHistogram(*(hfill_mTZZ[ic][ijet]));
        h_nvtxs_good[ic][ijet].emplace_back(
          Form("%s_nvtxs_good_%s_%s", sampleProcessedList.at(isample).name.data(), channel_labels.at(ic).Data(), label_Nj.Data()),
          sampleProcessedList.at(isample).label.data(),
          nbins_nvtxs_good, inf_nvtxs_good, sup_nvtxs_good
        ); hfill_nvtxs_good[ic][ijet] = &(h_nvtxs_good[ic][ijet].back()); hfill_nvtxs_good[ic][ijet]->GetXaxis()->SetTitle("N_{vtx}^{good}"); hfill_nvtxs_good[ic][ijet]->GetYaxis()->SetTitle("Events"); sampleProcessedList.at(isample).setupHistogram(*(hfill_nvtxs_good[ic][ijet]));
      }
    }

    float xsecScale = 1;
    if (sampleProcessedList.at(isample).name.find("VBFWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.064714*0.00847/0.00814; // BR

    if (sampleProcessedList.at(isample).name.find("ggZZ")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.541; // BR scale
    if (sampleProcessedList.at(isample).name.find("ggWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.616; // BR scale
    if (sampleProcessedList.at(isample).name.find("VBFZZ")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.312; // BR scale
    if (sampleProcessedList.at(isample).name.find("VBFWW")!=std::string::npos && sampleProcessedList.at(isample).name.find("Onshell")==std::string::npos) xsecScale *= 0.616*0.407/0.412; // BR scale

    auto const& tin = treelist.at(isample);
    int nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (has_ak4jets_inHEM1516 || has_ak8jets_inHEM1516) continue;

      if (pt_ll<55.f || pt_l1<25.f || pt_l2<25.f) continue;
      int index_channel = 0*is_ee + 1*is_mumu + 2*is_emu;

      if (sampleProcessedList.at(isample).name.find("ggWW_2l2nu_BSI")!=std::string::npos && std::abs(event_wgt)>150.) continue;
      if (sampleProcessedList.at(isample).name.find("ggWW_2l2nu_BSI_g4")!=std::string::npos && std::abs(event_wgt)>200.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_BSI")!=std::string::npos && std::abs(event_wgt)>500.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_BSI_g4")!=std::string::npos && std::abs(event_wgt)>30000.) continue;
      if (sampleProcessedList.at(isample).name.find("VBFWW_2l2nu_Sig_Onshell_g4")!=std::string::npos && std::abs(event_wgt)>200.) continue;
      event_wgt *= xsecScale;
      if (index_channel<2) event_wgt *= event_wgt_OSSF_triggers;
      else event_wgt *= event_wgt_OSDF_triggers;

      int index_Njets = std::min(nbins_Nj-1, static_cast<int>(event_Njets));
      int index_INorOUT=-1;
      if (std::abs(mass_ll-91.2f)<15.f){
        index_INorOUT = 0;
        if (event_Njets_btagged_loose!=0) continue;
        if (sampleProcessedList.at(isample).name.find("Data")!=std::string::npos && index_channel<2) continue;
      }
      else{
        index_INorOUT = 1;
        if ((event_Njets>0 && (useMediumBtagForOutSel ? event_Njets_btagged_medium : event_Njets_btagged_loose)==0) || (event_Njets==0 && (useMediumBtagForOutSel ? event_Njets20_btagged_medium : event_Njets20_btagged_loose)!=0)) continue;
      }

      if (
        !(
        (index_INorOUT==0 && pTmiss<125.f)
        ||
        (index_INorOUT==1 && pTmiss<85.f)
        )
        ) hfill_mll[index_channel][index_Njets]->Fill(mass_ll, event_wgt);

      if (flagSB==0 && std::abs(mass_ll-91.2f)>=15.f) continue;
      else if (flagSB==1 && !(mass_ll>=91.2f + 15.f && mass_ll<200.f)) continue;
      else if (flagSB==-1 && !(mass_ll<91.2f - 15.f && mass_ll>=40.f)) continue;

      hfill_pTmiss[index_channel][index_Njets]->Fill(pTmiss, event_wgt);

      if (
        (index_INorOUT==0 && pTmiss<125.f)
        ||
        (index_INorOUT==1 && pTmiss<85.f)
        ) continue;

      hfill_pTll[index_channel][index_Njets]->Fill(pt_ll, event_wgt);
      hfill_mTZZ[index_channel][index_Njets]->Fill(mTZZ, event_wgt);
      hfill_nvtxs_good[index_channel][index_Njets]->Fill(event_nvtxs_good, event_wgt);
    }
  }

  // Make stacked histograms to plot
  {
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        for (unsigned int ip=0; ip<h_pTmiss[ic][ijet].size(); ip++){
          TH1F& htmp_pTmiss_i = h_pTmiss[ic][ijet].at(ip);
          TH1F& htmp_pTll_i = h_pTll[ic][ijet].at(ip);
          TH1F& htmp_mTZZ_i = h_mTZZ[ic][ijet].at(ip);
          TH1F& htmp_mll_i = h_mll[ic][ijet].at(ip);
          TH1F& htmp_nvtxs_good_i = h_nvtxs_good[ic][ijet].at(ip);
          TString hname_pTmiss_i = htmp_pTmiss_i.GetName();

          if (hname_pTmiss_i.Contains("Data")) continue; // Do not stack anything on data

          bool const isHiggs_i = hname_pTmiss_i.Contains("BSI") || hname_pTmiss_i.Contains("Sig");
          bool const isOnshellHiggs_i = hname_pTmiss_i.Contains("Onshell");
          bool const isHiggsBSI_i = isHiggs_i && hname_pTmiss_i.Contains("BSI");
          bool const isHiggsSig_i = isHiggs_i && hname_pTmiss_i.Contains("Sig");
          bool const isSMHiggs_i = (isHiggsBSI_i || isHiggsSig_i) && !hname_pTmiss_i.Contains("_g4");
          bool const isPSHiggs_i = (isHiggsBSI_i || isHiggsSig_i) && hname_pTmiss_i.Contains("_g4");

          for (unsigned int jp=ip+1; jp<h_pTmiss[ic][ijet].size(); jp++){
            TH1F& htmp_pTmiss_j = h_pTmiss[ic][ijet].at(jp);
            TH1F& htmp_pTll_j = h_pTll[ic][ijet].at(jp);
            TH1F& htmp_mTZZ_j = h_mTZZ[ic][ijet].at(jp);
            TH1F& htmp_mll_j = h_mll[ic][ijet].at(jp);
            TH1F& htmp_nvtxs_good_j = h_nvtxs_good[ic][ijet].at(jp);
            TString hname_pTmiss_j = htmp_pTmiss_j.GetName();

            bool const isHiggs_j = hname_pTmiss_j.Contains("BSI") || hname_pTmiss_j.Contains("Sig");
            bool const isOnshellHiggs_j = hname_pTmiss_j.Contains("Onshell");
            bool const isHiggsBSI_j = isHiggs_j && hname_pTmiss_j.Contains("BSI");
            bool const isHiggsSig_j = isHiggs_j && hname_pTmiss_j.Contains("Sig");
            bool const isSMHiggs_j = (isHiggsBSI_j || isHiggsSig_j) && !hname_pTmiss_j.Contains("_g4");
            bool const isPSHiggs_j = (isHiggsBSI_j || isHiggsSig_j) && hname_pTmiss_j.Contains("_g4");

            bool doAdd_j_to_i = (
              !isHiggs_j
              ||
              (isSMHiggs_i && isSMHiggs_j)
              ||
              (isPSHiggs_i && isPSHiggs_j)
              );

            if (doAdd_j_to_i){
              htmp_pTmiss_i.Add(&htmp_pTmiss_j);
              htmp_pTll_i.Add(&htmp_pTll_j);
              htmp_mTZZ_i.Add(&htmp_mTZZ_j);
              htmp_mll_i.Add(&htmp_mll_j);
              htmp_nvtxs_good_i.Add(&htmp_nvtxs_good_j);
            }
          }

          if (!isHiggs_i){
            htmp_pTmiss_i.SetFillColor(htmp_pTmiss_i.GetLineColor()); htmp_pTmiss_i.SetFillStyle(3018);
            htmp_pTll_i.SetFillColor(htmp_pTll_i.GetLineColor()); htmp_pTll_i.SetFillStyle(3018);
            htmp_mTZZ_i.SetFillColor(htmp_mTZZ_i.GetLineColor()); htmp_mTZZ_i.SetFillStyle(3018);
            htmp_mll_i.SetFillColor(htmp_pTmiss_i.GetLineColor()); htmp_mll_i.SetFillStyle(3018);
            htmp_nvtxs_good_i.SetFillColor(htmp_nvtxs_good_i.GetLineColor()); htmp_nvtxs_good_i.SetFillStyle(3018);
          }
        }

        // Now make plots
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);

        TString hlabel_Nj;
        if (ijet==nbins_Nj-1) hlabel_Nj = Form("N_{j}#geq%i", ijet);
        else hlabel_Nj = Form("N_{j}=%i", ijet);

        TString selectionLabels;

        std::vector<TH1F*> hlist;
        std::vector<TString> hlabels;

        // mll
        selectionLabels = channel_hlabels.at(ic) + ", " + hlabel_Nj;
        if (useMediumBtagForOutSel) selectionLabels = selectionLabels + "|'In' b tag loose, 'out' medium";
        else selectionLabels = selectionLabels + "|'In' and 'out' b tag loose";
        selectionLabels = selectionLabels + "|p_{T}^{miss}>125, 85 GeV ('in', 'out')";
        if (flagSB==0){
          for (auto const& htmp:h_mll[ic][ijet]){
            TString hname = htmp.GetName();
            if (hname.Contains("VBFWW_2l2nu_BSI")) continue;
            if (hname.Contains("VBFWW_2l2nu_Sig_Onshell")) continue;
            //if (hname.Contains("Onshell")) continue;
            hlist.push_back((TH1F*) htmp.Clone(Form("copy_%s", htmp.GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp.GetTitle());
            if (!hname.Contains("Data")){ for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++)hlist.back()->SetBinError(ix, 0); }
          }
          makePlot(
            coutput_main, lumi,
            Form("c_mll_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
            hlist, hlabels,
            selectionLabels,
            "e1p", false, 2.0
          );
          for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
        }

        // Refresh labels
        selectionLabels = channel_hlabels.at(ic) + ", " + hlabel_Nj;
        if (flagSB==0) selectionLabels += "|'In' selection";
        else if (flagSB==1) selectionLabels += "|'Out'-like selection, upper SB";
        else if (flagSB==-1) selectionLabels += "|'Out'-like selection, lower SB";
        if (flagSB!=0){
          if (useMediumBtagForOutSel) selectionLabels = selectionLabels + "|Medium b tag";
          else selectionLabels = selectionLabels + "|Loose b tag";
        }
        else selectionLabels = selectionLabels + "|Loose b tag";

        // pTmiss
        {
          for (auto const& htmp:h_pTmiss[ic][ijet]){
            TString hname = htmp.GetName();
            if (hname.Contains("VBFWW_2l2nu_BSI")) continue;
            if (hname.Contains("VBFWW_2l2nu_Sig_Onshell")) continue;
            //if (hname.Contains("Onshell")) continue;
            hlist.push_back((TH1F*) htmp.Clone(Form("copy_%s", htmp.GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp.GetTitle());
            if (!hname.Contains("Data")){ for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++)hlist.back()->SetBinError(ix, 0); }
          }
          makePlot(
            coutput_main, lumi,
            Form("c_pTmiss_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
            hlist, hlabels,
            selectionLabels,
            "e1p", false, 1.5
          );
          for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
        }

        // pTll
        {
          for (auto const& htmp:h_pTll[ic][ijet]){
            TString hname = htmp.GetName();
            if (hname.Contains("VBFWW_2l2nu_BSI")) continue;
            if (hname.Contains("VBFWW_2l2nu_Sig_Onshell")) continue;
            //if (hname.Contains("Onshell")) continue;
            hlist.push_back((TH1F*) htmp.Clone(Form("copy_%s", htmp.GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp.GetTitle());
            if (!hname.Contains("Data")){ for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++)hlist.back()->SetBinError(ix, 0); }
          }
          makePlot(
            coutput_main, lumi,
            Form("c_pTll_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
            hlist, hlabels,
            selectionLabels,
            "e1p", false, 1.5
          );
          for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
        }

        // nvtxs_good
        {
          for (auto const& htmp:h_nvtxs_good[ic][ijet]){
            TString hname = htmp.GetName();
            if (hname.Contains("VBFWW_2l2nu_BSI")) continue;
            if (hname.Contains("VBFWW_2l2nu_Sig_Onshell")) continue;
            //if (hname.Contains("Onshell")) continue;
            hlist.push_back((TH1F*) htmp.Clone(Form("copy_%s", htmp.GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp.GetTitle());
            if (!hname.Contains("Data")){ for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++)hlist.back()->SetBinError(ix, 0); }
          }
          makePlot(
            coutput_main, lumi,
            Form("c_nvtxs_good_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
            hlist, hlabels,
            selectionLabels,
            "e1p", false, 1.5
          );
          for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
        }

        // mTZZ
        if (flagSB==0){
          for (auto const& htmp:h_mTZZ[ic][ijet]){
            TString hname = htmp.GetName();
            if (hname.Contains("VBFWW_2l2nu_BSI")) continue;
            if (hname.Contains("VBFWW_2l2nu_Sig_Onshell")) continue;
            //if (hname.Contains("Onshell")) continue;
            hlist.push_back((TH1F*) htmp.Clone(Form("copy_%s", htmp.GetName()))); hlist.back()->SetTitle("");
            hlabels.push_back(htmp.GetTitle());
            if (!hname.Contains("Data")){ for (int ix=0; ix<=hlist.back()->GetNbinsX()+1; ix++)hlist.back()->SetBinError(ix, 0); }
          }
          makePlot(
            coutput_main, lumi,
            Form("c_mTZZ_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data()),
            hlist, hlabels,
            selectionLabels,
            "e1p", false, 1.5
          );
          for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
        }

      }
    }
  }

  for (auto& tin:treelist) delete tin;
#undef BRANCH_COMMANDS
}


void get_Step1(TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200313");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString const cinput_main = "output/TTBarClosure/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Step1";
  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;
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

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(float, phimiss) \
  BRANCH_COMMAND(float, mTZZ) \
  BRANCH_COMMAND(float, mZZ) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
  BRANCH_COMMAND(unsigned int, event_Njets20) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged) \
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers) \
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers) \
  BRANCH_COMMAND(bool, is_ee) \
  BRANCH_COMMAND(bool, is_mumu) \
  BRANCH_COMMAND(bool, is_emu) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(int, id_l1) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(float, eff_l1) \
  BRANCH_COMMAND(int, id_l2) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, eff_l2)
#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kGray, 1, 2));
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kRed, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("DY_2l_M_10to50_ext", "DY", "DY_2l_M_10to50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_ext", "DY", "DY_2l_M_50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  /*
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("Total_noZpeak", "Total (no Z)", "total_noZpeak", -1, HistogramProperties((int) kGray, 1, 2));
  sampleList.emplace_back("Total", "Total", "total", -1, HistogramProperties((int) kGray, 1, 2));

  curdir->cd();

  std::vector<SampleSpecs> sampleProcessedList; sampleProcessedList.reserve(sampleList.size());
  std::vector<TChain*> treelist;
  for (auto const& sample:sampleList){
    if (sample.name.find("Total")!=std::string::npos) continue;
    TChain* tin = new TChain("SkimTree");
    TString cinput = cinput_main + '/' + sample.path.data() + "*.root";
    tin->Add(cinput);
    if (tin->GetEntries()==0){
      delete tin;
      continue;
    }
#define BRANCH_COMMAND(type, name) tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    treelist.push_back(tin);
    MELAout << "Extracted input file " << cinput << endl;

    sampleProcessedList.push_back(sample);
  }
  const size_t nsamples = sampleProcessedList.size();
  curdir->cd();

  constexpr int nbins_Nj=4;
  std::vector<TString> const channel_labels={ "ee", "mumu", "emu" };
  std::vector<TString> const in_out_labels={ "in", "out" };
  std::vector<double> const bins_mll={ 40, 50, 60, 70, 76.2, 86.2, 96.2, 106.2, 115, 125, 135, 145, 155, 165, 175, 185, 210 };
  std::vector<double> const bins_pTmiss={ 70, 80, 90, 100, 110, 125, 140, 170, 200, 350 };
  std::vector<double> const bins_pTll={ 55, 65, 75, 85, 95, 115, 135, 190, 250 };
  std::vector<double> const bins_pTl={ 25, 30, 40, 50, 65, 95, 130, 250 };
  std::vector<double> const bins_pTll_over_mll={ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0 };

  TFile* foutput = TFile::Open(coutput_main + "/histograms_all.root", "recreate");

  TDirectory* outdir_total = foutput->mkdir(sampleList.back().name.data());
  outdir_total->cd();
  TH1F* h_total_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_mTZZ[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_pTmiss[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_pTll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_pTe[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_pTmu[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH2F* h_total_pTll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH2F* h_total_pTll_over_mll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);
      for (size_t io=0; io<in_out_labels.size(); io++){
        MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

        h_total_mll[ic][ijet][io] = new TH1F(
          Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data()
        ); h_total_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_mll[ic][ijet][io]);
        h_total_mTZZ[ic][ijet][io] = new TH1F(
          Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          8, 180, 500
        ); h_total_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_total_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_mTZZ[ic][ijet][io]);
        h_total_pTmiss[ic][ijet][io] = new TH1F(
          Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTmiss.size()-1, bins_pTmiss.data()
        ); h_total_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_total_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTmiss[ic][ijet][io]);
        h_total_pTll[ic][ijet][io] = new TH1F(
          Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTll.size()-1, bins_pTll.data()
        ); h_total_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTll[ic][ijet][io]);
        h_total_pTe[ic][ijet][io] = new TH1F(
          Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_total_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTe[ic][ijet][io]);
        h_total_pTmu[ic][ijet][io] = new TH1F(
          Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_total_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTmu[ic][ijet][io]);
        h_total_pTll_vs_mll[ic][ijet][io] = new TH2F(
          Form("pTll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data(), bins_pTll.size()-1, bins_pTll.data()
        ); h_total_pTll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_pTll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_pTll_vs_mll[ic][ijet][io]->Sumw2();
        h_total_pTll_over_mll_vs_mll[ic][ijet][io] = new TH2F(
          Form("pTll_over_mll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data(), bins_pTll_over_mll.size()-1, bins_pTll_over_mll.data()
        ); h_total_pTll_over_mll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_pTll_over_mll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} / m_{ll}"); h_total_pTll_over_mll_vs_mll[ic][ijet][io]->Sumw2();
      }
    }
  }
  foutput->cd();

  TDirectory* outdir_total_noZpeak = foutput->mkdir(sampleList.at(sampleList.size()-2).name.data());
  outdir_total_noZpeak->cd();
  TH1F* h_total_noZpeak_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_noZpeak_mTZZ[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTmiss[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTe[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTmu[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH2F* h_total_noZpeak_pTll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH2F* h_total_noZpeak_pTll_over_mll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);
      for (size_t io=0; io<in_out_labels.size(); io++){
        MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

        h_total_noZpeak_mll[ic][ijet][io] = new TH1F(
          Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data()
        ); h_total_noZpeak_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_noZpeak_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_mll[ic][ijet][io]);
        h_total_noZpeak_mTZZ[ic][ijet][io] = new TH1F(
          Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          8, 180, 500
        ); h_total_noZpeak_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_total_noZpeak_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_mTZZ[ic][ijet][io]);
        h_total_noZpeak_pTmiss[ic][ijet][io] = new TH1F(
          Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTmiss.size()-1, bins_pTmiss.data()
        ); h_total_noZpeak_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_total_noZpeak_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTmiss[ic][ijet][io]);
        h_total_noZpeak_pTll[ic][ijet][io] = new TH1F(
          Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTll.size()-1, bins_pTll.data()
        ); h_total_noZpeak_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_noZpeak_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTll[ic][ijet][io]);
        h_total_noZpeak_pTe[ic][ijet][io] = new TH1F(
          Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_noZpeak_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_total_noZpeak_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTe[ic][ijet][io]);
        h_total_noZpeak_pTmu[ic][ijet][io] = new TH1F(
          Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_noZpeak_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_total_noZpeak_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTmu[ic][ijet][io]);
        h_total_noZpeak_pTll_vs_mll[ic][ijet][io] = new TH2F(
          Form("pTll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data(), bins_pTll.size()-1, bins_pTll.data()
        ); h_total_noZpeak_pTll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_noZpeak_pTll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_noZpeak_pTll_vs_mll[ic][ijet][io]->Sumw2();
        h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io] = new TH2F(
          Form("pTll_over_mll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data(), bins_pTll_over_mll.size()-1, bins_pTll_over_mll.data()
        ); h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} / m_{ll}"); h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io]->Sumw2();
      }
    }
  }
  foutput->cd();

  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << isample << " (" << sampleProcessedList.at(isample).name << ")" << endl;
    auto const& tin = treelist.at(isample);

    foutput->cd();
    TDirectory* outdir = foutput->mkdir(sampleProcessedList.at(isample).name.data());

    TH1F* h_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH1F* h_mTZZ[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH1F* h_pTmiss[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH1F* h_pTll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH1F* h_pTe[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH1F* h_pTmu[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH2F* h_pTll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    TH2F* h_pTll_over_mll_vs_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    outdir->cd();
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);
        for (size_t io=0; io<in_out_labels.size(); io++){
          MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

          h_mll[ic][ijet][io] = new TH1F(
            Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_mll.size()-1, bins_mll.data()
          ); h_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_mll[ic][ijet][io]);
          h_mTZZ[ic][ijet][io] = new TH1F(
            Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            8, 180, 500
          ); h_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_mTZZ[ic][ijet][io]);
          h_pTmiss[ic][ijet][io] = new TH1F(
            Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTmiss.size()-1, bins_pTmiss.data()
          ); h_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTmiss[ic][ijet][io]);
          h_pTll[ic][ijet][io] = new TH1F(
            Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTll.size()-1, bins_pTll.data()
          ); h_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTll[ic][ijet][io]);
          h_pTe[ic][ijet][io] = new TH1F(
            Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTl.size()-1, bins_pTl.data()
          ); h_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTe[ic][ijet][io]);
          h_pTmu[ic][ijet][io] = new TH1F(
            Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTl.size()-1, bins_pTl.data()
          ); h_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTmu[ic][ijet][io]);
          h_pTll_vs_mll[ic][ijet][io] = new TH2F(
            Form("pTll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_mll.size()-1, bins_mll.data(), bins_pTll.size()-1, bins_pTll.data()
          ); h_pTll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_pTll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_pTll_vs_mll[ic][ijet][io]->Sumw2();
          h_pTll_over_mll_vs_mll[ic][ijet][io] = new TH2F(
            Form("pTll_over_mll_vs_mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_mll.size()-1, bins_mll.data(), bins_pTll_over_mll.size()-1, bins_pTll_over_mll.data()
          ); h_pTll_over_mll_vs_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_pTll_over_mll_vs_mll[ic][ijet][io]->GetYaxis()->SetTitle("p_{T}^{ll} / m_{ll}"); h_pTll_over_mll_vs_mll[ic][ijet][io]->Sumw2();
        }
      }
    }
    foutput->cd();

    int nEntries = tin->GetEntries();
    MELAout << "Looping over " << nEntries << " entries now..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);

      if (pt_l1<25.f || pt_l2<25.f) continue;
      if (pt_ll<55.f) continue;
      int index_channel = 0*is_ee + 1*is_mumu + 2*is_emu;

      if (index_channel<2) event_wgt *= event_wgt_OSSF_triggers;
      else event_wgt *= event_wgt_OSDF_triggers;

      if (mass_ll<40.f || mass_ll>=200.f) continue;

      int index_Njets = std::min(nbins_Nj-1, static_cast<int>(event_Njets));
      int index_INorOUT=-1;
      if (std::abs(mass_ll-91.2f)<15.f){
        index_INorOUT = 0;
        if (pTmiss<125.f) continue;
        if (event_Njets_btagged!=0) continue;
      }
      else{
        index_INorOUT = 1;
        if (pTmiss<70.f) continue;
        if ((event_Njets>0 && event_Njets_btagged==0) || (event_Njets==0 && event_Njets20_btagged!=0)) continue;
      }

      h_mll[index_channel][index_Njets][index_INorOUT]->Fill(mass_ll, event_wgt);
      h_mTZZ[index_channel][index_Njets][index_INorOUT]->Fill(mTZZ, event_wgt);
      h_pTmiss[index_channel][index_Njets][index_INorOUT]->Fill(pTmiss, event_wgt);
      h_pTll[index_channel][index_Njets][index_INorOUT]->Fill(pt_ll, event_wgt);
      if (std::abs(id_l1)==11) h_pTe[index_channel][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt);
      else h_pTmu[index_channel][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt);
      if (std::abs(id_l2)==11) h_pTe[index_channel][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt);
      else h_pTmu[index_channel][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt);
      h_pTll_vs_mll[index_channel][index_Njets][index_INorOUT]->Fill(mass_ll, pt_ll, event_wgt);
      h_pTll_over_mll_vs_mll[index_channel][index_Njets][index_INorOUT]->Fill(mass_ll, pt_ll/mass_ll, event_wgt);
    }

    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        for (size_t io=0; io<in_out_labels.size(); io++){
          if (
            sampleProcessedList.at(isample).name.find("DY")==std::string::npos
            &&
            sampleProcessedList.at(isample).name.find("ZZ_2l2nu")==std::string::npos
            &&
            sampleProcessedList.at(isample).name.find("WZ_3l")==std::string::npos
            ){
            h_total_noZpeak_mll[ic][ijet][io]->Add(h_mll[ic][ijet][io], 1.);
            h_total_noZpeak_mTZZ[ic][ijet][io]->Add(h_mTZZ[ic][ijet][io], 1.);
            h_total_noZpeak_pTmiss[ic][ijet][io]->Add(h_pTmiss[ic][ijet][io], 1.);
            h_total_noZpeak_pTll[ic][ijet][io]->Add(h_pTll[ic][ijet][io], 1.);
            h_total_noZpeak_pTe[ic][ijet][io]->Add(h_pTe[ic][ijet][io], 1.);
            h_total_noZpeak_pTmu[ic][ijet][io]->Add(h_pTmu[ic][ijet][io], 1.);
            h_total_noZpeak_pTll_vs_mll[ic][ijet][io]->Add(h_pTll_vs_mll[ic][ijet][io], 1.);
            h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io]->Add(h_pTll_over_mll_vs_mll[ic][ijet][io], 1.);
          }

          outdir->WriteTObject(h_mll[ic][ijet][io]); h_total_mll[ic][ijet][io]->Add(h_mll[ic][ijet][io], 1.); delete h_mll[ic][ijet][io];
          outdir->WriteTObject(h_mTZZ[ic][ijet][io]); h_total_mTZZ[ic][ijet][io]->Add(h_mTZZ[ic][ijet][io], 1.); delete h_mTZZ[ic][ijet][io];
          outdir->WriteTObject(h_pTmiss[ic][ijet][io]); h_total_pTmiss[ic][ijet][io]->Add(h_pTmiss[ic][ijet][io], 1.); delete h_pTmiss[ic][ijet][io];
          outdir->WriteTObject(h_pTll[ic][ijet][io]); h_total_pTll[ic][ijet][io]->Add(h_pTll[ic][ijet][io], 1.); delete h_pTll[ic][ijet][io];
          outdir->WriteTObject(h_pTe[ic][ijet][io]); h_total_pTe[ic][ijet][io]->Add(h_pTe[ic][ijet][io], 1.); delete h_pTe[ic][ijet][io];
          outdir->WriteTObject(h_pTmu[ic][ijet][io]); h_total_pTmu[ic][ijet][io]->Add(h_pTmu[ic][ijet][io], 1.); delete h_pTmu[ic][ijet][io];
          outdir->WriteTObject(h_pTll_vs_mll[ic][ijet][io]); h_total_pTll_vs_mll[ic][ijet][io]->Add(h_pTll_vs_mll[ic][ijet][io], 1.); delete h_pTll_vs_mll[ic][ijet][io];
          outdir->WriteTObject(h_pTll_over_mll_vs_mll[ic][ijet][io]); h_total_pTll_over_mll_vs_mll[ic][ijet][io]->Add(h_pTll_over_mll_vs_mll[ic][ijet][io], 1.); delete h_pTll_over_mll_vs_mll[ic][ijet][io];
        }
      }
    }
    outdir->Close();
  }

  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      for (size_t io=0; io<in_out_labels.size(); io++){
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_mll[ic][ijet][io]); delete h_total_noZpeak_mll[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_mTZZ[ic][ijet][io]); delete h_total_noZpeak_mTZZ[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTmiss[ic][ijet][io]); delete h_total_noZpeak_pTmiss[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTll[ic][ijet][io]); delete h_total_noZpeak_pTll[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTe[ic][ijet][io]); delete h_total_noZpeak_pTe[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTmu[ic][ijet][io]); delete h_total_noZpeak_pTmu[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTll_vs_mll[ic][ijet][io]); delete h_total_noZpeak_pTll_vs_mll[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io]); delete h_total_noZpeak_pTll_over_mll_vs_mll[ic][ijet][io];

        outdir_total->WriteTObject(h_total_mll[ic][ijet][io]); delete h_total_mll[ic][ijet][io];
        outdir_total->WriteTObject(h_total_mTZZ[ic][ijet][io]); delete h_total_mTZZ[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTmiss[ic][ijet][io]); delete h_total_pTmiss[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTll[ic][ijet][io]); delete h_total_pTll[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTe[ic][ijet][io]); delete h_total_pTe[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTmu[ic][ijet][io]); delete h_total_pTmu[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTll_vs_mll[ic][ijet][io]); delete h_total_pTll_vs_mll[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTll_over_mll_vs_mll[ic][ijet][io]); delete h_total_pTll_over_mll_vs_mll[ic][ijet][io];
      }
    }
  }
  outdir_total_noZpeak->Close();
  outdir_total->Close();

  foutput->Close();
  for (auto& tin:treelist) delete tin;
#undef BRANCH_COMMANDS
}

void get_Step2(TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200313");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString const cinput_core = "output/TTBarClosure/SkimTrees/" + strdate;
  TString cinput_main = cinput_core + "/Step1";
  TString coutput_main = cinput_core + "/Step2";
  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;
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

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(float, phimiss) \
  BRANCH_COMMAND(float, mTZZ) \
  BRANCH_COMMAND(float, mZZ) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
  BRANCH_COMMAND(unsigned int, event_Njets20) \
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged) \
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers) \
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers) \
  BRANCH_COMMAND(bool, is_ee) \
  BRANCH_COMMAND(bool, is_mumu) \
  BRANCH_COMMAND(bool, is_emu) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(int, id_l1) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(float, eff_l1) \
  BRANCH_COMMAND(int, id_l2) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, eff_l2)
#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kGray, 1, 2));
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kRed, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("DY_2l_M_10to50_ext", "DY", "DY_2l_M_10to50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_ext", "DY", "DY_2l_M_50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
  /*
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("Total_noZpeak", "Total (no Z)", "total_noZpeak", -1, HistogramProperties((int) kGray, 1, 2));
  sampleList.emplace_back("Total", "Total", "total", -1, HistogramProperties((int) kGray, 1, 2));

  curdir->cd();

  std::vector<SampleSpecs> sampleProcessedList; sampleProcessedList.reserve(sampleList.size());
  std::vector<TChain*> treelist;
  for (auto const& sample:sampleList){
    if (sample.name.find("Total")!=std::string::npos) continue;
    TChain* tin = new TChain("SkimTree");
    TString cinput = cinput_core + '/' + sample.path.data() + "*.root";
    tin->Add(cinput);
    if (tin->GetEntries()==0){
      delete tin;
      continue;
    }
#define BRANCH_COMMAND(type, name) tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    treelist.push_back(tin);
    MELAout << "Extracted input file " << cinput << endl;

    sampleProcessedList.push_back(sample);
  }
  const size_t nsamples = sampleProcessedList.size();
  curdir->cd();

  constexpr int nbins_Nj=4;
  std::vector<TString> const channel_in_labels={ "ee", "mumu", "emu" };
  std::vector<TString> const channel_labels={ "ee", "mumu", "emu_ee_like", "emu_mumu_like" };
  std::vector<TString> const in_out_labels={ "in", "out" };
  std::vector<double> const bins_mll={ 40, 50, 60, 70, 76.2, 86.2, 96.2, 106.2, 115, 125, 135, 145, 155, 165, 175, 185, 210 };
  std::vector<double> const bins_pTmiss={ 70, 80, 90, 100, 110, 125, 140, 170, 200, 350 };
  std::vector<double> const bins_pTll={ 55, 65, 75, 85, 95, 115, 135, 190, 250 };
  std::vector<double> const bins_pTl={ 25, 30, 40, 50, 65, 95, 130, 250 };
  std::vector<double> const bins_pTll_over_mll={ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0 };

  TFile* finput_rewgt = TFile::Open(cinput_main + "/histograms_all.root", "read");
  finput_rewgt->cd();
  TH1F* h_rewgt_pTe[2][nbins_Nj]={ { 0 } }; // Channel (ee, emu), Njets, out-like selection
  TH1F* h_rewgt_pTmu[2][nbins_Nj]={ { 0 } }; // Channel (mumu, emu), Njets, out-like selection
  for (int ijet=0; ijet<nbins_Nj; ijet++){
    TString label_Nj;
    if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
    else label_Nj = Form("Nj_eq_%i", ijet);

    h_rewgt_pTe[0][ijet] = (TH1F*) finput_rewgt->Get(Form("%s/pTe_%s_%s_%s", sampleList.at(sampleList.size()-2).name.data(), channel_in_labels.at(0).Data(), label_Nj.Data(), in_out_labels.at(1).Data()));
    h_rewgt_pTmu[0][ijet] = (TH1F*) finput_rewgt->Get(Form("%s/pTmu_%s_%s_%s", sampleList.at(sampleList.size()-2).name.data(), channel_in_labels.at(1).Data(), label_Nj.Data(), in_out_labels.at(1).Data()));
    h_rewgt_pTe[1][ijet] = (TH1F*) finput_rewgt->Get(Form("%s/pTe_%s_%s_%s", sampleList.at(sampleList.size()-2).name.data(), channel_in_labels.at(2).Data(), label_Nj.Data(), in_out_labels.at(1).Data()));
    h_rewgt_pTmu[1][ijet] = (TH1F*) finput_rewgt->Get(Form("%s/pTmu_%s_%s_%s", sampleList.at(sampleList.size()-2).name.data(), channel_in_labels.at(2).Data(), label_Nj.Data(), in_out_labels.at(1).Data()));
  }
  curdir->cd();

  TFile* foutput = TFile::Open(coutput_main + "/histograms_all.root", "recreate");

  TDirectory* outdir_total = foutput->mkdir(sampleList.back().name.data());
  outdir_total->cd();
  TH1F* h_total_mll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_mTZZ[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_pTmiss[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_pTll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_pTe[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_pTmu[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);
      for (size_t io=0; io<in_out_labels.size(); io++){
        MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

        h_total_mll[ic][ijet][io] = new TH1F(
          Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data()
        ); h_total_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_mll[ic][ijet][io]);
        h_total_mTZZ[ic][ijet][io] = new TH1F(
          Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          8, 180, 500
        ); h_total_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_total_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_mTZZ[ic][ijet][io]);
        h_total_pTmiss[ic][ijet][io] = new TH1F(
          Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTmiss.size()-1, bins_pTmiss.data()
        ); h_total_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_total_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTmiss[ic][ijet][io]);
        h_total_pTll[ic][ijet][io] = new TH1F(
          Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTll.size()-1, bins_pTll.data()
        ); h_total_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTll[ic][ijet][io]);
        h_total_pTe[ic][ijet][io] = new TH1F(
          Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_total_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTe[ic][ijet][io]);
        h_total_pTmu[ic][ijet][io] = new TH1F(
          Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_total_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.back().setupHistogram(h_total_pTmu[ic][ijet][io]);
      }
    }
  }
  foutput->cd();

  TDirectory* outdir_total_noZpeak = foutput->mkdir(sampleList.at(sampleList.size()-2).name.data());
  outdir_total_noZpeak->cd();
  TH1F* h_total_noZpeak_mll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_noZpeak_mTZZ[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTmiss[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTe[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_total_noZpeak_pTmu[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);
      for (size_t io=0; io<in_out_labels.size(); io++){
        MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

        h_total_noZpeak_mll[ic][ijet][io] = new TH1F(
          Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_mll.size()-1, bins_mll.data()
        ); h_total_noZpeak_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_total_noZpeak_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_mll[ic][ijet][io]);
        h_total_noZpeak_mTZZ[ic][ijet][io] = new TH1F(
          Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          8, 180, 500
        ); h_total_noZpeak_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_total_noZpeak_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_mTZZ[ic][ijet][io]);
        h_total_noZpeak_pTmiss[ic][ijet][io] = new TH1F(
          Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTmiss.size()-1, bins_pTmiss.data()
        ); h_total_noZpeak_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_total_noZpeak_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTmiss[ic][ijet][io]);
        h_total_noZpeak_pTll[ic][ijet][io] = new TH1F(
          Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTll.size()-1, bins_pTll.data()
        ); h_total_noZpeak_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_total_noZpeak_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTll[ic][ijet][io]);
        h_total_noZpeak_pTe[ic][ijet][io] = new TH1F(
          Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_noZpeak_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_total_noZpeak_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTe[ic][ijet][io]);
        h_total_noZpeak_pTmu[ic][ijet][io] = new TH1F(
          Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
          bins_pTl.size()-1, bins_pTl.data()
        ); h_total_noZpeak_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_total_noZpeak_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleList.at(sampleList.size()-2).setupHistogram(h_total_noZpeak_pTmu[ic][ijet][io]);
      }
    }
  }
  foutput->cd();

  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << isample << " (" << sampleProcessedList.at(isample).name << ")" << endl;
    auto const& tin = treelist.at(isample);

    foutput->cd();
    TDirectory* outdir = foutput->mkdir(sampleProcessedList.at(isample).name.data());

    TH1F* h_mll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    TH1F* h_mTZZ[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    TH1F* h_pTmiss[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    TH1F* h_pTll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    TH1F* h_pTe[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    TH1F* h_pTmu[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
    outdir->cd();
    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        TString label_Nj;
        if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
        else label_Nj = Form("Nj_eq_%i", ijet);
        for (size_t io=0; io<in_out_labels.size(); io++){
          MELAout << "Creating histograms for " << channel_labels.at(ic) << ", " << label_Nj << ", " << in_out_labels.at(io) << endl;

          h_mll[ic][ijet][io] = new TH1F(
            Form("mll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_mll.size()-1, bins_mll.data()
          ); h_mll[ic][ijet][io]->GetXaxis()->SetTitle("m_{ll} (GeV)"); h_mll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_mll[ic][ijet][io]);
          h_mTZZ[ic][ijet][io] = new TH1F(
            Form("mTZZ_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            8, 180, 500
          ); h_mTZZ[ic][ijet][io]->GetXaxis()->SetTitle("m^{ZZ}_{T} (GeV)"); h_mTZZ[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_mTZZ[ic][ijet][io]);
          h_pTmiss[ic][ijet][io] = new TH1F(
            Form("pTmiss_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTmiss.size()-1, bins_pTmiss.data()
          ); h_pTmiss[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)"); h_pTmiss[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTmiss[ic][ijet][io]);
          h_pTll[ic][ijet][io] = new TH1F(
            Form("pTll_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTll.size()-1, bins_pTll.data()
          ); h_pTll[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{ll} (GeV)"); h_pTll[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTll[ic][ijet][io]);
          h_pTe[ic][ijet][io] = new TH1F(
            Form("pTe_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTl.size()-1, bins_pTl.data()
          ); h_pTe[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{e} (GeV)"); h_pTe[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTe[ic][ijet][io]);
          h_pTmu[ic][ijet][io] = new TH1F(
            Form("pTmu_%s_%s_%s", channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()), "",
            bins_pTl.size()-1, bins_pTl.data()
          ); h_pTmu[ic][ijet][io]->GetXaxis()->SetTitle("p_{T}^{#mu} (GeV)"); h_pTmu[ic][ijet][io]->GetYaxis()->SetTitle("N_{events}"); sampleProcessedList.at(isample).setupHistogram(h_pTmu[ic][ijet][io]);
        }
      }
    }
    foutput->cd();

    int nEntries = tin->GetEntries();
    MELAout << "Looping over " << nEntries << " entries now..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);

      if (pt_l1<25.f || pt_l2<25.f) continue;
      if (pt_ll<55.f) continue;
      int index_channel = 0*is_ee + 1*is_mumu + 2*is_emu;

      if (index_channel<2) event_wgt *= event_wgt_OSSF_triggers;
      else event_wgt *= event_wgt_OSDF_triggers;

      if (mass_ll<40.f || mass_ll>=200.f) continue;

      int index_Njets = std::min(nbins_Nj-1, static_cast<int>(event_Njets));
      int index_INorOUT=-1;
      if (std::abs(mass_ll-91.2f)<15.f){
        index_INorOUT = 0;
        if (pTmiss<125.f) continue;
        if (event_Njets_btagged!=0) continue;
      }
      else{
        index_INorOUT = 1;
        if (pTmiss<70.f) continue;
        if ((event_Njets>0 && event_Njets_btagged==0) || (event_Njets==0 && event_Njets20_btagged!=0)) continue;
      }

      float event_wgt_extra=1;
      if (index_channel==2){
        event_wgt_extra=1;
        float const& tmp_pTe = (std::abs(id_l1)==11 ? pt_l1 : pt_l2);
        float const& tmp_pTmu = (std::abs(id_l1)==13 ? pt_l1 : pt_l2);
        float wgt_e_eelike = h_rewgt_pTe[0][index_Njets]->GetBinContent(h_rewgt_pTe[0][index_Njets]->GetXaxis()->FindBin(tmp_pTe));
        float wgt_e_emulike = h_rewgt_pTe[1][index_Njets]->GetBinContent(h_rewgt_pTe[1][index_Njets]->GetXaxis()->FindBin(tmp_pTe));
        if (wgt_e_emulike!=0.f) event_wgt_extra *= wgt_e_eelike/wgt_e_emulike;
        else event_wgt_extra=0;
        float wgt_mu_eelike = h_rewgt_pTe[0][index_Njets]->GetBinContent(h_rewgt_pTe[0][index_Njets]->GetXaxis()->FindBin(tmp_pTmu));
        float wgt_mu_emulike = h_rewgt_pTmu[1][index_Njets]->GetBinContent(h_rewgt_pTmu[1][index_Njets]->GetXaxis()->FindBin(tmp_pTmu));
        if (wgt_mu_emulike!=0.f) event_wgt_extra *= wgt_mu_eelike/wgt_mu_emulike;
        else event_wgt_extra=0;
      }
      h_mll[index_channel][index_Njets][index_INorOUT]->Fill(mass_ll, event_wgt*event_wgt_extra);
      h_mTZZ[index_channel][index_Njets][index_INorOUT]->Fill(mTZZ, event_wgt*event_wgt_extra);
      h_pTmiss[index_channel][index_Njets][index_INorOUT]->Fill(pTmiss, event_wgt*event_wgt_extra);
      h_pTll[index_channel][index_Njets][index_INorOUT]->Fill(pt_ll, event_wgt*event_wgt_extra);
      if (std::abs(id_l1)==11) h_pTe[index_channel][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt*event_wgt_extra);
      else h_pTmu[index_channel][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt*event_wgt_extra);
      if (std::abs(id_l2)==11) h_pTe[index_channel][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt*event_wgt_extra);
      else h_pTmu[index_channel][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt*event_wgt_extra);

      if (index_channel==2){
        event_wgt_extra=1;
        float const& tmp_pTe = (std::abs(id_l1)==11 ? pt_l1 : pt_l2);
        float const& tmp_pTmu = (std::abs(id_l1)==13 ? pt_l1 : pt_l2);
        float wgt_e_mumulike = h_rewgt_pTmu[0][index_Njets]->GetBinContent(h_rewgt_pTmu[0][index_Njets]->GetXaxis()->FindBin(tmp_pTe));
        float wgt_e_emulike = h_rewgt_pTe[1][index_Njets]->GetBinContent(h_rewgt_pTe[1][index_Njets]->GetXaxis()->FindBin(tmp_pTe));
        if (wgt_e_emulike!=0.f) event_wgt_extra *= wgt_e_mumulike/wgt_e_emulike;
        else event_wgt_extra=0;
        float wgt_mu_mumulike = h_rewgt_pTmu[0][index_Njets]->GetBinContent(h_rewgt_pTmu[0][index_Njets]->GetXaxis()->FindBin(tmp_pTmu));
        float wgt_mu_emulike = h_rewgt_pTmu[1][index_Njets]->GetBinContent(h_rewgt_pTmu[1][index_Njets]->GetXaxis()->FindBin(tmp_pTmu));
        if (wgt_mu_emulike!=0.f) event_wgt_extra *= wgt_mu_mumulike/wgt_mu_emulike;
        else event_wgt_extra=0;

        h_mll[index_channel+1][index_Njets][index_INorOUT]->Fill(mass_ll, event_wgt*event_wgt_extra);
        h_mTZZ[index_channel+1][index_Njets][index_INorOUT]->Fill(mTZZ, event_wgt*event_wgt_extra);
        h_pTmiss[index_channel+1][index_Njets][index_INorOUT]->Fill(pTmiss, event_wgt*event_wgt_extra);
        h_pTll[index_channel+1][index_Njets][index_INorOUT]->Fill(pt_ll, event_wgt*event_wgt_extra);
        if (std::abs(id_l1)==11) h_pTe[index_channel+1][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt*event_wgt_extra);
        else h_pTmu[index_channel+1][index_Njets][index_INorOUT]->Fill(pt_l1, event_wgt*event_wgt_extra);
        if (std::abs(id_l2)==11) h_pTe[index_channel+1][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt*event_wgt_extra);
        else h_pTmu[index_channel+1][index_Njets][index_INorOUT]->Fill(pt_l2, event_wgt*event_wgt_extra);
      }
    }

    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=0; ijet<nbins_Nj; ijet++){
        for (size_t io=0; io<in_out_labels.size(); io++){
          if (
            sampleProcessedList.at(isample).name.find("DY")==std::string::npos
            &&
            sampleProcessedList.at(isample).name.find("ZZ_2l2nu")==std::string::npos
            &&
            sampleProcessedList.at(isample).name.find("WZ_3l")==std::string::npos
            ){
            h_total_noZpeak_mll[ic][ijet][io]->Add(h_mll[ic][ijet][io], 1.);
            h_total_noZpeak_mTZZ[ic][ijet][io]->Add(h_mTZZ[ic][ijet][io], 1.);
            h_total_noZpeak_pTmiss[ic][ijet][io]->Add(h_pTmiss[ic][ijet][io], 1.);
            h_total_noZpeak_pTll[ic][ijet][io]->Add(h_pTll[ic][ijet][io], 1.);
            h_total_noZpeak_pTe[ic][ijet][io]->Add(h_pTe[ic][ijet][io], 1.);
            h_total_noZpeak_pTmu[ic][ijet][io]->Add(h_pTmu[ic][ijet][io], 1.);
          }

          outdir->WriteTObject(h_mll[ic][ijet][io]); h_total_mll[ic][ijet][io]->Add(h_mll[ic][ijet][io], 1.); delete h_mll[ic][ijet][io];
          outdir->WriteTObject(h_mTZZ[ic][ijet][io]); h_total_mTZZ[ic][ijet][io]->Add(h_mTZZ[ic][ijet][io], 1.); delete h_mTZZ[ic][ijet][io];
          outdir->WriteTObject(h_pTmiss[ic][ijet][io]); h_total_pTmiss[ic][ijet][io]->Add(h_pTmiss[ic][ijet][io], 1.); delete h_pTmiss[ic][ijet][io];
          outdir->WriteTObject(h_pTll[ic][ijet][io]); h_total_pTll[ic][ijet][io]->Add(h_pTll[ic][ijet][io], 1.); delete h_pTll[ic][ijet][io];
          outdir->WriteTObject(h_pTe[ic][ijet][io]); h_total_pTe[ic][ijet][io]->Add(h_pTe[ic][ijet][io], 1.); delete h_pTe[ic][ijet][io];
          outdir->WriteTObject(h_pTmu[ic][ijet][io]); h_total_pTmu[ic][ijet][io]->Add(h_pTmu[ic][ijet][io], 1.); delete h_pTmu[ic][ijet][io];
        }
      }
    }
    outdir->Close();
  }

  for (size_t ic=0; ic<channel_labels.size(); ic++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      for (size_t io=0; io<in_out_labels.size(); io++){
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_mll[ic][ijet][io]); delete h_total_noZpeak_mll[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_mTZZ[ic][ijet][io]); delete h_total_noZpeak_mTZZ[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTmiss[ic][ijet][io]); delete h_total_noZpeak_pTmiss[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTll[ic][ijet][io]); delete h_total_noZpeak_pTll[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTe[ic][ijet][io]); delete h_total_noZpeak_pTe[ic][ijet][io];
        outdir_total_noZpeak->WriteTObject(h_total_noZpeak_pTmu[ic][ijet][io]); delete h_total_noZpeak_pTmu[ic][ijet][io];

        outdir_total->WriteTObject(h_total_mll[ic][ijet][io]); delete h_total_mll[ic][ijet][io];
        outdir_total->WriteTObject(h_total_mTZZ[ic][ijet][io]); delete h_total_mTZZ[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTmiss[ic][ijet][io]); delete h_total_pTmiss[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTll[ic][ijet][io]); delete h_total_pTll[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTe[ic][ijet][io]); delete h_total_pTe[ic][ijet][io];
        outdir_total->WriteTObject(h_total_pTmu[ic][ijet][io]); delete h_total_pTmu[ic][ijet][io];
      }
    }
  }
  outdir_total_noZpeak->Close();
  outdir_total->Close();

  foutput->Close();
  finput_rewgt->Close();

  for (auto& tin:treelist) delete tin;
#undef BRANCH_COMMANDS
}


void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool adjustYLow,
  float factorYHigh
){
  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  bool hasData = false;

  std::vector<bool> hHasErrors;

  double ymin = 0;
  if (adjustYLow) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (hname.Contains("Data")) hasData = true;
    bool hasErrors=false;
    for (int ix=1; ix<=hist->GetNbinsX(); ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (adjustYLow) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //MELAout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : 1.5);
  else ymax /= (factorYHigh>0.f ? factorYHigh : 1.5);
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

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
  if (!hasData) text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  else text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
  text = pt->AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  bool firstHist = true;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString const& hlabel = hlabels.at(is);

    hist->SetTitle("");
    if (hlabel!="") legend->AddEntry(hist, hlabel, "f");

    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    if (firstHist){
      hist->Draw(stropt);
      firstHist = false;
    }
    else{
      hist->Draw(stropt+"same");
    }
  }

  // Re-draw data
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (!hname.Contains("Data")) continue;
    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    hist->Draw(stropt+"same");
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


void makePlots(TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200313");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString const cinput_core = "output/TTBarClosure/SkimTrees/" + strdate;
  TString cinput_step1 = cinput_core + "/Step1";
  TString cinput_step2 = cinput_core + "/Step2";
  TString coutput_main = cinput_core + "/Plots";
  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;
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

  TString const proc_total_noZpeak = "Total_noZpeak";
  constexpr int nbins_Nj=4;
  std::vector<TString> const channel_in_labels={ "ee", "mumu", "emu" };
  std::vector<TString> const channel_labels={ "ee", "mumu", "emu_ee_like", "emu_mumu_like" };
  std::vector<TString> const in_out_labels={ "in", "out" };

  std::vector<TString> const channel_in_hlabels={ "ee", "#mu#mu", "e#mu" };
  std::vector<TString> const channel_hlabels={ "ee", "#mu#mu", "e#mu (ee-like)", "e#mu (#mu#mu-like)" };
  std::vector<TString> const in_out_hlabels={ "'in'", "'out'" };

  TFile* finput_step1 = TFile::Open(cinput_step1 + "/histograms_all.root", "read");
  TFile* finput_step2 = TFile::Open(cinput_step2 + "/histograms_all.root", "read");

  TH1F* h_step1_mll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_step1_mTZZ[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_step1_pTmiss[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_step1_pTll[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_step1_pTe[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
  TH1F* h_step1_pTmu[3][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection

  TH1F* h_step2_mll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_step2_mTZZ[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_step2_pTmiss[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_step2_pTll[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_step2_pTe[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection
  TH1F* h_step2_pTmu[4][nbins_Nj][2]={ { { 0 } } }; // Channel (ee, mumu, emu [ee], emu [mumu]), Njets, in/out-like selection

  for (size_t io=0; io<in_out_labels.size(); io++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);
      for (size_t ic=0; ic<channel_in_labels.size(); ic++){
        h_step1_mll[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/mll_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step1_mTZZ[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/mTZZ_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step1_pTll[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/pTll_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step1_pTmiss[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/pTmiss_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step1_pTe[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/pTe_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step1_pTmu[ic][ijet][io] = (TH1F*) finput_step1->Get(Form("%s/pTmu_%s_%s_%s", proc_total_noZpeak.Data(), channel_in_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
      }
      for (size_t ic=0; ic<channel_labels.size(); ic++){
        h_step2_mll[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/mll_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step2_mTZZ[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/mTZZ_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step2_pTmiss[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/pTmiss_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step2_pTll[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/pTll_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step2_pTe[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/pTe_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
        h_step2_pTmu[ic][ijet][io] = (TH1F*) finput_step2->Get(Form("%s/pTmu_%s_%s_%s", proc_total_noZpeak.Data(), channel_labels.at(ic).Data(), label_Nj.Data(), in_out_labels.at(io).Data()));
      }
    }
  }

  // Plot raw distributions from Step1
  for (size_t io=0; io<in_out_labels.size(); io++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);

      TString hlabel_Nj;
      if (ijet==nbins_Nj-1) hlabel_Nj = Form("N_{j}#geq%i", ijet);
      else hlabel_Nj = Form("N_{j}=%i", ijet);

      TString selectionLabels = hlabel_Nj + ", " + in_out_hlabels.at(io) + " (step 1)";
      std::vector<TH1F*> hlist;
      std::vector<TString> hlabels;

      // mTZZ
      hlist.push_back((TH1F*) h_step1_mTZZ[0][ijet][io]->Clone(Form("copy_%s", h_step1_mTZZ[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_mTZZ[1][ijet][io]->Clone(Form("copy_%s", h_step1_mTZZ[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_mTZZ[2][ijet][io]->Clone(Form("copy_%s", h_step1_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0){
        hlist.push_back((TH1F*) h_step1_mTZZ[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_mTZZ[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step1_mTZZ_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTmiss
      hlist.push_back((TH1F*) h_step1_pTmiss[0][ijet][io]->Clone(Form("copy_%s", h_step1_pTmiss[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTmiss[1][ijet][io]->Clone(Form("copy_%s", h_step1_pTmiss[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTmiss[2][ijet][io]->Clone(Form("copy_%s", h_step1_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTmiss[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTmiss[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step1_pTmiss_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTll
      hlist.push_back((TH1F*) h_step1_pTll[0][ijet][io]->Clone(Form("copy_%s", h_step1_pTll[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTll[1][ijet][io]->Clone(Form("copy_%s", h_step1_pTll[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTll[2][ijet][io]->Clone(Form("copy_%s", h_step1_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTll[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTll[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step1_pTll_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTe
      hlist.push_back((TH1F*) h_step1_pTe[0][ijet][io]->Clone(Form("copy_%s", h_step1_pTe[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      //hlist.push_back((TH1F*) h_step1_pTe[1][ijet][io]->Clone(Form("copy_%s", h_step1_pTe[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTe[2][ijet][io]->Clone(Form("copy_%s", h_step1_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTe[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(2.*h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        //hlist.push_back((TH1F*) h_step1_pTe[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        //hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step1_pTe_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTmu
      //hlist.push_back((TH1F*) h_step1_pTmu[0][ijet][io]->Clone(Form("copy_%s", h_step1_pTmu[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTmu[1][ijet][io]->Clone(Form("copy_%s", h_step1_pTmu[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step1_pTmu[2][ijet][io]->Clone(Form("copy_%s", h_step1_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0){
        //hlist.push_back((TH1F*) h_step1_pTmu[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        //hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTmu[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(2.*h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step1_pTmu_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      if (io==1){
        // mll
        hlist.push_back((TH1F*) h_step1_mll[0][ijet][io]->Clone(Form("copy_%s", h_step1_mll[0][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.push_back((TH1F*) h_step1_mll[1][ijet][io]->Clone(Form("copy_%s", h_step1_mll[1][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.push_back((TH1F*) h_step1_mll[2][ijet][io]->Clone(Form("copy_%s", h_step1_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_in_hlabels.at(2)); hlist.back()->SetLineColor(kViolet); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        if (io==0){
          hlist.push_back((TH1F*) h_step1_mll[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
          hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
          hlist.push_back((TH1F*) h_step1_mll[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
          hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        }
        makePlot(
          coutput_main, lumi,
          Form("c_step1_mll_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
          hlist, hlabels,
          selectionLabels
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }
    }
  }

  // Plot raw distributions from Step2
  for (size_t io=0; io<in_out_labels.size(); io++){
    for (int ijet=0; ijet<nbins_Nj; ijet++){
      TString label_Nj;
      if (ijet==nbins_Nj-1) label_Nj = Form("Nj_geq_%i", ijet);
      else label_Nj = Form("Nj_eq_%i", ijet);

      TString hlabel_Nj;
      if (ijet==nbins_Nj-1) hlabel_Nj = Form("N_{j}#geq%i", ijet);
      else hlabel_Nj = Form("N_{j}=%i", ijet);

      TString selectionLabels = hlabel_Nj + ", " + in_out_hlabels.at(io) + " (step 2)";
      std::vector<TH1F*> hlist;
      std::vector<TString> hlabels;

      // mTZZ
      hlist.push_back((TH1F*) h_step2_mTZZ[0][ijet][io]->Clone(Form("copy_%s", h_step2_mTZZ[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_mTZZ[1][ijet][io]->Clone(Form("copy_%s", h_step2_mTZZ[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_mTZZ[2][ijet][io]->Clone(Form("copy_%s", h_step2_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
      hlist.push_back((TH1F*) h_step2_mTZZ[3][ijet][io]->Clone(Form("copy_%s", h_step2_mTZZ[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
      if (io==0){
        hlist.push_back((TH1F*) h_step1_mTZZ[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_mTZZ[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_mTZZ[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step2_mTZZ_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTmiss
      hlist.push_back((TH1F*) h_step2_pTmiss[0][ijet][io]->Clone(Form("copy_%s", h_step2_pTmiss[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTmiss[1][ijet][io]->Clone(Form("copy_%s", h_step2_pTmiss[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTmiss[2][ijet][io]->Clone(Form("copy_%s", h_step2_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
      hlist.push_back((TH1F*) h_step2_pTmiss[3][ijet][io]->Clone(Form("copy_%s", h_step2_pTmiss[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTmiss[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTmiss[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTmiss[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step2_pTmiss_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTll
      hlist.push_back((TH1F*) h_step2_pTll[0][ijet][io]->Clone(Form("copy_%s", h_step2_pTll[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTll[1][ijet][io]->Clone(Form("copy_%s", h_step2_pTll[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTll[2][ijet][io]->Clone(Form("copy_%s", h_step2_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
      hlist.push_back((TH1F*) h_step2_pTll[3][ijet][io]->Clone(Form("copy_%s", h_step2_pTll[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTll[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTll[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step2_pTll_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTe
      hlist.push_back((TH1F*) h_step2_pTe[0][ijet][io]->Clone(Form("copy_%s", h_step2_pTe[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTe[1][ijet][io]->Clone(Form("copy_%s", h_step2_pTe[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTe[2][ijet][io]->Clone(Form("copy_%s", h_step2_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
      hlist.push_back((TH1F*) h_step2_pTe[3][ijet][io]->Clone(Form("copy_%s", h_step2_pTe[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTe[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTe[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTe[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step2_pTe_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      // pTmu
      hlist.push_back((TH1F*) h_step2_pTmu[0][ijet][io]->Clone(Form("copy_%s", h_step2_pTmu[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTmu[1][ijet][io]->Clone(Form("copy_%s", h_step2_pTmu[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      hlist.push_back((TH1F*) h_step2_pTmu[2][ijet][io]->Clone(Form("copy_%s", h_step2_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
      hlist.push_back((TH1F*) h_step2_pTmu[3][ijet][io]->Clone(Form("copy_%s", h_step2_pTmu[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
      if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
      if (io==0){
        hlist.push_back((TH1F*) h_step1_pTmu[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step1_pTmu[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_pTmu[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
      }
      makePlot(
        coutput_main, lumi,
        Form("c_step2_pTmu_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
        hlist, hlabels,
        selectionLabels
      );
      for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();

      if (io==1){
        // mll
        hlist.push_back((TH1F*) h_step2_mll[0][ijet][io]->Clone(Form("copy_%s", h_step2_mll[0][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(0)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.push_back((TH1F*) h_step2_mll[1][ijet][io]->Clone(Form("copy_%s", h_step2_mll[1][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(1)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(1); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        hlist.push_back((TH1F*) h_step2_mll[2][ijet][io]->Clone(Form("copy_%s", h_step2_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2)); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        if (io==0) hlist.back()->Scale(h_step2_mll[0][ijet][1-io]->Integral(1, h_step2_mll[0][ijet][1-io]->GetNbinsX()) / h_step2_mll[2][ijet][1-io]->Integral(1, h_step2_mll[2][ijet][1-io]->GetNbinsX()));
        hlist.push_back((TH1F*) h_step2_mll[3][ijet][io]->Clone(Form("copy_%s", h_step2_mll[3][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3)); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(7); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
        if (io==0) hlist.back()->Scale(h_step2_mll[1][ijet][1-io]->Integral(1, h_step2_mll[1][ijet][1-io]->GetNbinsX()) / h_step2_mll[3][ijet][1-io]->Integral(1, h_step2_mll[3][ijet][1-io]->GetNbinsX()));
        if (io==0){
          hlist.push_back((TH1F*) h_step1_mll[2][ijet][io]->Clone(Form("copy_ee_like_%s", h_step1_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(2) + " (step 1)"); hlist.back()->SetLineColor(kRed); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
          hlist.back()->Scale(h_step1_mll[0][ijet][1-io]->Integral(1, h_step1_mll[0][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
          hlist.push_back((TH1F*) h_step1_mll[2][ijet][io]->Clone(Form("copy_mumu_like_%s", h_step1_mll[2][ijet][io]->GetName()))); hlabels.push_back(channel_hlabels.at(3) + " (step 1)"); hlist.back()->SetLineColor(kBlue); hlist.back()->SetLineStyle(2); hlist.back()->SetLineWidth(2); hlist.back()->SetMarkerColor(hlist.back()->GetLineColor());
          hlist.back()->Scale(h_step1_mll[1][ijet][1-io]->Integral(1, h_step1_mll[1][ijet][1-io]->GetNbinsX()) / h_step1_mll[2][ijet][1-io]->Integral(1, h_step1_mll[2][ijet][1-io]->GetNbinsX()));
        }
        makePlot(
          coutput_main, lumi,
          Form("c_step2_mll_%s_%s", label_Nj.Data(), in_out_labels.at(io).Data()),
          hlist, hlabels,
          selectionLabels
        );
        for (auto& hh:hlist) delete hh; hlist.clear(); hlabels.clear();
      }
    }
  }


  finput_step2->Close();
  finput_step1->Close();
}
