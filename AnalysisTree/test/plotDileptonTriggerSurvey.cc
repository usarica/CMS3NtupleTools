#include "common_includes.h"
#include "OffshellCutflow.h"
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
  TString const& canvasname,
  TH2F* hist,
  TString selectionLabels
);


using namespace SystematicsHelpers;
void plot(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int id_ll=0
){
  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;
  if (id_ll!=-13*13 && id_ll!=-13*11 && id_ll!=-11*11) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  bool const periodIsLikeData = SampleHelpers::testDataPeriodIsLikeData();

  // Jet ID options
  bool applyPUIdToAK4Jets=true, applyTightLeptonVetoIdToAK4Jets=false, useJetOverlapStripping = false;
  // MET options
  bool use_MET_Puppi=false;
  bool use_MET_XYCorr=true, use_MET_JERCorr=false, use_MET_ParticleMomCorr=true;
  bool use_MET_p4Preservation=false;
  bool use_MET_corrections=false;
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers;
  if (id_ll==-13*11) requiredTriggers = std::vector<TriggerHelpers::TriggerType>{
    TriggerHelpers::kMuEle, TriggerHelpers::kMuEle_Extra,
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt, TriggerHelpers::kSingleEle_Prescaled
  };
  else if (id_ll==-11*11) requiredTriggers = std::vector<TriggerHelpers::TriggerType>{
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt, TriggerHelpers::kSingleEle_L1EG,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt, TriggerHelpers::kSingleEle_Prescaled
  };
  else requiredTriggers = std::vector<TriggerHelpers::TriggerType>{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt
  };
  if (SampleHelpers::getDataYear()==2016 && id_ll==-13*13) requiredTriggers.push_back(TriggerHelpers::kDoubleMu_Prescaled);
  auto triggerCheckList_all = TriggerHelpers::getHLTMenus(requiredTriggers);
  std::vector<std::string> triggerCheckList_baseline;
  switch (SampleHelpers::getDataYear()){
  case 2016:
  {
    if (id_ll==-13*11) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele25_eta2p1_WPTight_Gsf_v*", "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"
    };
    else if (id_ll==-13*13) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
    };
    else triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele25_eta2p1_WPTight_Gsf_v*", "HLT_Ele27_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*", "HLT_DoubleEle33_CaloIdL_MW_v*", "HLT_DoublePhoton60_v*"
    };
    break;
  }
  case 2017:
  {
    if (id_ll==-13*11) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_IsoMu27_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"
    };
    else if (id_ll==-13*13) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_IsoMu27_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"
    };
    else triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele35_WPTight_Gsf_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_DoubleEle33_CaloIdL_MW_v*", "HLT_DoublePhoton70_v*"
    };
    break;
  }
  case 2018:
  {
    if (id_ll==-13*11) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele32_WPTight_Gsf_v*", "HLT_Photon200_v*",
      "HLT_IsoMu24_v*",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"
    };
    else if (id_ll==-13*13) triggerCheckList_baseline=std::vector<std::string>{
      "HLT_IsoMu24_v*",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"
    };
    else triggerCheckList_baseline=std::vector<std::string>{
      "HLT_Ele32_WPTight_Gsf_v*", "HLT_Photon200_v*",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", "HLT_DoubleEle25_CaloIdL_MW_v*", "HLT_DoublePhoton70_v*"
    };
    break;
  }
  default:
    break;
  }
  for (auto& strig:triggerCheckList_all) HelperFunctions::replaceString(strig, "_v*", "_v");
  for (auto& strig:triggerCheckList_baseline) HelperFunctions::replaceString(strig, "_v*", "_v");
  std::vector<std::string> triggerCheckList_extras;
  for (auto const& strig:triggerCheckList_all){ if (!HelperFunctions::checkListVariable(triggerCheckList_baseline, strig)) triggerCheckList_extras.push_back(strig); }
  MELAout << "Baseline triggers: " << triggerCheckList_baseline << endl;
  MELAout << "Triggers to check: " << triggerCheckList_extras << endl;

  TDirectory* curdir = gDirectory;

  bool isData = false;
  {
    std::vector<TString> sampledirs;
    SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
    isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  }
  if (isData && theGlobalSyst!=sNominal) return;

  TString cinput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(cinput, "_MINIAOD", "");
  bool containsPeriodString = false;
  for (auto const& pp:validDataPeriods){
    if (cinput.Contains(pp)){ containsPeriodString = true; break; }
  }
  bool isGGSig = cinput.Contains("/GluGluH");
  bool isVBFSig = cinput.Contains("/VBF");
  bool isTT = cinput.Contains("/TTTo2L2Nu");
  bool isQQZZ = cinput.Contains("/ZZTo2L2Nu");
  bool isDY = cinput.Contains("/DY");
  float xsec_scale = 1;
  if (isGGSig) xsec_scale = 0.00813704*0.541;
  if (isVBFSig) xsec_scale = 0.00846586*0.541;
  bool apply_pTmiss = (isGGSig || isVBFSig || isTT || isQQZZ);
  bool apply_dPhi = (isTT || isQQZZ);
  bool apply_btagging = (isGGSig || isVBFSig || isDY || isQQZZ);
  bool apply_mll = (isGGSig || isVBFSig || isDY || isQQZZ);
  bool apply_mllwide = (isTT);
  TString strproc, strproctitle;
  if (isGGSig) strproc = Form("gg#rightarrowH#rightarrow2%s2#nu", (id_ll==-121 ? "e" : "#mu"));
  if (isVBFSig) strproc = Form("VBF H#rightarrow2%s2#nu", (id_ll==-121 ? "e" : "#mu"));
  if (isTT) strproc = Form("t#bar{t}#rightarrow%s", (id_ll==-121 ? "2e2#nu" : (id_ll==-169 ? "2#mu2#nu" : "e#mu2#nu")));
  if (isQQZZ) strproc = Form("q#bar{q}#rightarrow2%s2#nu", (id_ll==-121 ? "e" : "#mu"));
  if (isDY) strproc = Form("Z#rightarrow2%s", (id_ll==-121 ? "e" : "#mu"));
  if (isGGSig) strproctitle = "ggH";
  if (isVBFSig) strproctitle = "VBFH";
  if (isTT) strproctitle = "tt2l2nu";
  if (isQQZZ) strproctitle = "qqZZ";
  if (isDY) strproctitle = "DY";

  TString cinput_base_dir;
  /*if (!SampleHelpers::checkRunOnCondor()) cinput_base_dir = "output/";
  else*/ cinput_base_dir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/";

  // Set input directory
  TChain* tin = nullptr;
  {
    TString cinput_main =
      cinput_base_dir + "DileptonTriggerSurvey/SkimTrees/" + strdate
      + "/AK4Jets"
      + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
      + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
      + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
      + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
      + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
    if (!isData) cinput_main = cinput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
    cinput_main = cinput_main
      + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
      + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
    if (!isData) cinput_main = cinput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");

    tin = new TChain("SkimTree");
    if (isData || containsPeriodString){
      for (auto const& pp:validDataPeriods){
        if (periodIsLikeData && pp!=SampleHelpers::getDataPeriod()) continue;

        TString cinput_tmp = cinput;
        if (containsPeriodString) HelperFunctions::replaceString(cinput_tmp, SampleHelpers::getDataPeriod(), pp);
        else if (isData && cinput_tmp==Form("Run%i", SampleHelpers::getDataYear())) cinput_tmp = Form("Run%s", pp.Data());

        TString strinput = Form("%s/%s/%s", cinput_main.Data(), pp.Data(), cinput_tmp.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
        strinput += ".root";

        MELAout << "Adding input " << strinput << endl;
        tin->Add(strinput);
      }
    }
    else{
      TString strinput = Form("%s/%s/%s", cinput_main.Data(), period.Data(), cinput.Data());
      strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
      strinput += ".root";

      MELAout << "Adding input " << strinput << endl;
      tin->Add(strinput);
    }

    MELAout << "Tree has a total of " << tin->GetEntries() << " entries..." << endl;
    curdir->cd();
  }

  TString coutput_main = 
    TString("output/DileptonTriggerSurveyChecks/") + strdate
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  coutput_main = coutput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  coutput_main = coutput_main + "/" + period;
  if (id_ll==-11*11) coutput_main = coutput_main + "/ee";
  else if (id_ll==-13*13) coutput_main = coutput_main + "/mumu";
  else coutput_main = coutput_main + "/mue";
  coutput_main = coutput_main + "/" + cinput;

  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + "/histograms.root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

#define BRANCH_OPTIONAL_COMMANDS \
  BRANCH_COMMAND(float, KFactor_QCD_NNLO_ggZZ_Sig_Nominal) \
  BRANCH_COMMAND(float, p_Gen_CPStoBWPropRewgt) \
  BRANCH_COMMAND(float, p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM) \
  BRANCH_COMMAND(float, p_Gen_JJEW_SIG_ghv1_1_MCFM)
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

  std::unordered_map<std::string, float> strig_val_map;
  for (auto const& strig:triggerCheckList_all) strig_val_map[strig]=0;
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 1; bool exists_##NAME = false;
  BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  {
#define BRANCH_COMMAND(TYPE, NAME) exists_##NAME = tin->GetBranchStatus(#NAME)==1;
    BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) if (exists_##NAME){ MELAout << "Booking " << #NAME << "..." << endl; tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME); }
    BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) MELAout << "Booking " << #NAME << "..." << endl; tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (auto& it:strig_val_map){
      auto const& strig = it.first;
      MELAout << "Booking " << strig << "..." << endl;
      tin->SetBranchStatus(Form("event_wgt_%s", strig.data()), 1);
      tin->SetBranchAddress(Form("event_wgt_%s", strig.data()), &(it.second));
    }
  }

  std::vector<float> binning_pt{
    5, 15, 25, 30, 35, 40, 50, 75, 100, 130, 170, 200, 500
  };
  
  std::unordered_map<std::string, TH2F*> hist_trig_map;
  {
    hist_trig_map["baseline"] = new TH2F(
      "h_triggers_baseline", "Baseline",
      binning_pt.size()-1, binning_pt.data(),
      binning_pt.size()-1, binning_pt.data()
    );
    hist_trig_map["all_extras"] = new TH2F(
      "h_triggers_all_extras", "All extras",
      binning_pt.size()-1, binning_pt.data(),
      binning_pt.size()-1, binning_pt.data()
    );
    for (auto const& strig:triggerCheckList_extras){
      hist_trig_map[strig] = new TH2F(
        Form("h_%s", strig.data()), Form("%s%s", strig.data(), "*"),
        binning_pt.size()-1, binning_pt.data(),
        binning_pt.size()-1, binning_pt.data()
      );
    }
  }
  {
    TString xlabel = (id_ll==-13*13 ? "Leading-p_{T}^{#mu}" : (id_ll==-13*11 ? "p_{T}^{#mu}" : "Leading-p_{T}^{e}"));
    TString ylabel = (id_ll==-13*13 ? "Subleading-p_{T}^{#mu}" : (id_ll==-13*11 ? "p_{T}^{e}" : "Subleading-p_{T}^{e}"));
    for (auto& it:hist_trig_map){
      auto& hh=it.second;
      hh->Sumw2();
      hh->GetXaxis()->SetTitle(xlabel);
      hh->GetYaxis()->SetTitle(ylabel);
      hh->SetOption("colz");
    }
  }

  {
    using namespace OffshellCutflow;

    foutput->cd();

    constexpr unsigned int ncuts = 7;
    std::vector<TString> const cutlabels{
      "Id cut",
      "pT, eta, mass cuts on leptons or boson",
      "dPhi ll, MET",
      "dPhi ll+jets, MET",
      "dPhi jets, MET",
      "b-tagging",
      "pTmiss and pT boson"
    };
    assert(cutlabels.size()==ncuts);

    double sum_wgts_cuts[ncuts]={ 0 };
    double sum_wgts_total = 0;
    const int nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEvent(ev);
      HelperFunctions::progressbar(ev, nEntries);

      float wgt = event_wgt*event_wgt_triggers*event_wgt_SFs*xsec_scale;
#define BRANCH_COMMAND(TYPE, NAME) if (exists_##NAME) wgt *= NAME;
      BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
      sum_wgts_total += wgt;

      if (dilepton_id!=id_ll) continue;
      sum_wgts_cuts[0] += wgt;

      float pTx, pTy;
      if (id_ll==-13*11){
        pTx = (std::abs(leptons_id->front())==13 ? leptons_pt->front() : leptons_pt->back());
        pTy = (std::abs(leptons_id->front())==13 ? leptons_pt->back() : leptons_pt->front());
      }
      else{
        pTx = std::max(leptons_pt->front(), leptons_pt->back());
        pTy = std::min(leptons_pt->back(), leptons_pt->front());
      }

      if (apply_mll && !check_mll(dilepton_mass, (dilepton_id!=-11*13))) continue;
      if (apply_mllwide && !(dilepton_mass>=70. && dilepton_mass<200.)) continue;
      sum_wgts_cuts[1] += wgt;

      if (apply_dPhi && !check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      sum_wgts_cuts[2] += wgt;

      if (apply_dPhi && !check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      sum_wgts_cuts[3] += wgt;

      if (apply_dPhi && !check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss)) continue;
      sum_wgts_cuts[4] += wgt;

      if (apply_btagging && event_n_ak4jets_pt30_btagged_loose>0) continue;
      sum_wgts_cuts[5] += wgt;

      bool pass_boson_pt = check_pTboson(dilepton_pt);
      bool pass_pTmiss = check_pTmiss(event_pTmiss);
      if (isData && id_ll!=-13*11) pass_pTmiss = !pass_pTmiss; // Invert requirement
      if ((apply_pTmiss && !pass_pTmiss) || !pass_boson_pt) continue;
      sum_wgts_cuts[6] += wgt;

      wgt = std::abs(wgt);
      bool pass_baseline = false;
      bool pass_extras = false;
      for (auto const& strig:triggerCheckList_baseline){
        if (strig_val_map[strig]>0.f){
          pass_baseline=true;
          break;
        }
      }
      if (pass_baseline) hist_trig_map["baseline"]->Fill(pTx, pTy, wgt);
      else{
        for (auto const& strig:triggerCheckList_extras){
          if (strig_val_map[strig]==0.f) continue;
          hist_trig_map[strig]->Fill(pTx, pTy, wgt);
          pass_extras = true;
        }
      }
      if (pass_extras) hist_trig_map["all_extras"]->Fill(pTx, pTy, wgt);
    }
    MELAout << "Accumulated sum of weights:" << sum_wgts_total << endl;
    for (unsigned int icut=0; icut<ncuts; icut++){
      MELAout << "\t- " << cutlabels.at(icut) << ": " << sum_wgts_cuts[icut]
        << ", efficiency: " << sum_wgts_cuts[icut]/sum_wgts_cuts[0]
        << ",  recursive eff.: " << sum_wgts_cuts[icut]/sum_wgts_cuts[(icut==0 ? 0 : icut-1)]
        << endl;
    }

    delete tin;
  }

  foutput->cd();
  for (auto& it:hist_trig_map) HelperFunctions::wipeOverUnderFlows(it.second, false, true);
  for (auto& it:hist_trig_map){
    if (it.first!="baseline") it.second->Divide(hist_trig_map["baseline"]);
    foutput->WriteTObject(it.second);

    makePlot(
      coutput_main, lumi,
      Form("%s_%s", it.second->GetName(), strproctitle.Data()),
      it.second,
      strproc
    );
  }
  for (auto& it:hist_trig_map) delete it.second;
  foutput->Close();
}

void plotAll(TString period, TString prodVersion, TString strdate){
  std::vector<TString> procnames;
  if (period=="2016") procnames={
    "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUgenV698_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM",
    "/ZZTo2L2Nu_13TeV_powheg_pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM",
    "/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
  };
  else if (period=="2017") procnames={
    "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    "/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM"
  };
  else if (period=="2018") procnames={
    "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    "/VBF_HToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM",
    "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM",
    "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM"
  };
  for (auto const& procname:procnames){
    std::vector<int> channels{ -121,-169 };
    if (procname.Contains("/TT")) channels.push_back(-11*13);
    for (auto const& channel:channels) plot(procname, period, prodVersion, strdate, channel);
  }
}

void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString const& canvasname,
  TH2F* hist,
  TString selectionLabels
){
  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  TString strig = hist->GetTitle();
  hist->SetTitle("");
  selectionList.push_back(strig);

  TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 1000, 800);
  canvas->cd();
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
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
  text = pt->AddText(0.83, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  hist->GetXaxis()->SetRangeUser(25, hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX())-1e-5);
  hist->GetYaxis()->SetRangeUser(25, hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY())-1e-5);
  hist->Draw("colz");

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
  canvas->Close();
}
