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
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts="hist",
  bool adjustYLow=false,
  float factorYHigh=-1
);


float getAbsWeightThresholdByNeff(TTree* tree, std::vector<float*> const& vals, double thr_Neff, TVar::VerbosityLevel verbosity=TVar::ERROR){
  float res = -1;

  int nEntries = tree->GetEntries();
  thr_Neff = std::min(thr_Neff, double(nEntries)/3.*2.);
  unsigned int npos = 0;
  double Neff = 0;
  double sum_wgts[2]={ 0 }; // [0]: w, [1]: w^2
  std::vector<float> smallest_weights;
  if (verbosity>=TVar::ERROR) MELAout << "getAbsWeightThresholdByNeff: Determining the weight thresholds (Neff threshold = " << thr_Neff << ", number of events = " << nEntries << ")..." << endl;
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt_combined = 1;
    for (auto const& val:vals) wgt_combined *= std::abs(*val);

    HelperFunctions::addByLowest(smallest_weights, wgt_combined, false);
    if (ev%100000 == 0 || ev == nEntries-1){
      sum_wgts[0] = sum_wgts[1] = 0;
      npos = 0;
      for (auto const& wgt:smallest_weights){
        sum_wgts[0] += wgt;
        sum_wgts[1] += wgt*wgt;
        Neff = std::pow(sum_wgts[0], 2) / sum_wgts[1];
        npos++;
        if (Neff>=thr_Neff) break;
      }
      if (verbosity>=TVar::ERROR) MELAout << "\t- Current Neff = " << Neff << " over " << ev+1 << " events..." << endl;
    }
    if (Neff>=thr_Neff) break;
  }

  if (!smallest_weights.empty()){
    res = (sum_wgts[0] + std::sqrt(sum_wgts[1]*Neff)) / (Neff-1.);
    if (verbosity>=TVar::ERROR){
      unsigned int nVeto = 0;
      for (auto const& wgt:smallest_weights){
        if (wgt<res) continue;
        nVeto++;
      }

      MELAout
        << "\t- " << res
        << " is the default weight threshold calculated from sN=" << sum_wgts[0] << ", vN=" << sum_wgts[1] << ", nN=" << Neff
        << " (N=" << npos << " / " << smallest_weights.size() << ", wN=" << smallest_weights.at(npos-1) << ", wLast=" << smallest_weights.back()
        << "). Expected fraction of vetos: " << ((double) nVeto) / ((double) smallest_weights.size())
        << endl;
    }
  }
  else{
    if (verbosity>=TVar::INFO) MELAout << "\t- No weight threshold is found." << endl;
  }

  return res;
}


void getDataSampleDirs(std::vector<TString>& strsamples){
  SystematicsHelpers::SystematicVariationTypes const syst = SystematicsHelpers::sNominal;
  TString strSyst = SystematicsHelpers::getSystName(syst).data();

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  bool isDataLikePeriod = SampleHelpers::testDataPeriodIsLikeData();

  for (auto const& period:validDataPeriods){
    if (isDataLikePeriod && period!=SampleHelpers::getDataPeriod()) continue;
    strsamples.push_back(Form("%s/Run%s_%s%s", period.Data(), period.Data(), strSyst.Data(), ".root"));
  }
}
void getMCSampleDirs(std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> >& strsamples, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst){
  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString strPeriod = SampleHelpers::getDataPeriod();

  std::vector< std::pair< TString, std::vector<TString> > > sampleSpecs;
  switch (SampleHelpers::getDataYear()){
  case 2016:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_10to50" }
      },
      {
        "DY_2l",{ "DY_2l_M_50"}
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_inclusive", "WJets_lnu_inclusive_ext" }
      }
    };
    break;
  case 2017:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_10to50", "DY_2l_M_10to50_ext" }
      },
      {
        "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext", "DY_2l_M_50_ext2" }
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_0j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_1j", "WJets_lnu_1j_ext" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_2j", "WJets_lnu_2j_ext" }
      }
    };
    break;
  case 2018:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l", { "DY_2l_M_10to50", "DY_2l_M_10to50_ext" }
      },
      {
        "DY_2l", { "DY_2l_M_50", "DY_2l_M_50_ext" }
      },
      {
        "TT_2l2nu", { "TT_2l2nu" }
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_0j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_1j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_2j" }
      }
    };
    break;
  }
  for (auto const& s:sampleSpecs){
    std::vector<TString> sdirs;
    std::vector<std::pair<TString, TString>> sname_dir_pairs;
    for (auto const& strSampleSet:s.second) SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sdirs);
    sname_dir_pairs.reserve(sdirs.size());
    for (auto const& sname:sdirs){
      TString cinput = SampleHelpers::getSampleIdentifier(sname);
      HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
      HelperFunctions::replaceString(cinput, "_MINIAOD", "");
      cinput = cinput + "_*" + strSyst + ".root";
      cinput = strPeriod + "/" + cinput;
      sname_dir_pairs.emplace_back(sname, cinput);
    }
    strsamples.emplace_back(s.first, sname_dir_pairs);
  }
}

using namespace SystematicsHelpers;
void getDistributions(
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  unsigned int istep,
  bool useGenMatchedLeptons = false
){
  constexpr bool useJetOverlapStripping=false;
  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;
  // MET options
  constexpr bool use_MET_Puppi=false;
  constexpr bool use_MET_XYCorr=true;
  constexpr bool use_MET_JERCorr=true;
  constexpr bool use_MET_ParticleMomCorr=true;
  constexpr bool use_MET_p4Preservation=true;
  constexpr bool use_MET_corrections=true;

  if (istep>1) return;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "hadoop_skims", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> const strChannelNames{ "mumu", "ee", "mue", "mue_rewgt_mumu", "mue_rewgt_ee" };
  std::vector<TString> const strChannelTitles{ "#mu#mu", "ee", "#mue (un-rewgt.)", "#mue (rewgt. #mu#mu)", "#mue (rewgt. ee)" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> const strNjetsNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2", "Nj_geq_0" };
  std::vector<TString> const strNjetsTitles{ "N_{j}=0", "N_{j}=1", "N_{j}#geq2", "N_{j}#geq0" };
  const unsigned int nbins_njets = strNjetsNames.size();

  std::vector<TString> const strBtaggingRegionNames{ "Nbloose_eq_0", "Nbmed_geq_1" };
  std::vector<TString> const strBtaggingRegionTitles{ "N_{b}^{loose}=0", "N_{b}^{medium}#geq1" };
  const unsigned int nbins_nbtagged = strBtaggingRegionNames.size();

  constexpr float thr_corr_pTmiss = 100;

  TString const coutput_main = "output/NRBEstimates/" + strdate + "/Histograms/" + period;
  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  TriggerScaleFactorHandler triggerSFHandler;

  TString stroutput = coutput_main + Form("/histograms_%s", strSyst.Data());
  if (useGenMatchedLeptons) stroutput = stroutput + "_GenMatchedLeptons";
  stroutput = stroutput + Form("_Step%u", istep);
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();

  // Acquire corrections from previous steps
  std::vector<TGraph*> tg_incorr_list[2][2][nbins_njets]; // [Data, non-res MC][mumu, ee][Nj=0, 1, >=2, >=0]
  std::vector<TFile*> finputs_prevstep;
  if (istep>0){
    auto const& strBtaggingRegionName = strBtaggingRegionNames.at(1);
    for (unsigned int jstep=0; jstep<istep-1; jstep++){
      TString strinput_prevstep = stroutput;
      TString ownstep = Form("Step%u", istep);
      TString prevstep = Form("Step%u", jstep);
      HelperFunctions::replaceString(strinput_prevstep, ownstep, prevstep);
      TFile* ftmp = TFile::Open(strinput_prevstep, "read"); finputs_prevstep.push_back(ftmp);
      ftmp->cd();
      for (unsigned int ip=0; ip<2; ip++){
        TString strProcess = (ip==0 ? "Data" : "AllMC_NonRes");
        for (unsigned int ic=0; ic<2; ic++){
          auto const& strChannelName = strChannelNames.at(ic);
          for (unsigned int ij=0; ij<nbins_njets; ij++){
            auto const& strNjetsName = strNjetsNames.at(ij);

            TString strCutName = strChannelName + "_" + strNjetsName + "_" + strBtaggingRegionName;
            TString strname = Form("tg_corr_%s_%s_%s_Nominal", "pTmiss", strCutName.Data(), strProcess.Data());

            tg_incorr_list[ip][ic][ij].push_back((TGraph*) ftmp->Get(strname));
          }
        }
      }
      foutput->cd();
    }
  }

  TString const cinput_main_MC =
    "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + strdate
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY")
    + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER")
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default")
    + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  std::unordered_map<TChain*, double> norm_map;
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst);
  std::vector<std::pair<TString, TChain*>> samples_all;
  std::vector<TString> sgroups;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    auto const& sname_sfname_pairs = sgroup_sname_sfname_pair.second;
    if (!HelperFunctions::checkListVariable(sgroups, sgroup)) sgroups.push_back(sgroup);
    std::vector<TChain*> tins_collected;
    for (auto const& sname_sfname_pair:sname_sfname_pairs){
      auto const& sname = sname_sfname_pair.first;
      auto const& sfname = sname_sfname_pair.second;
      TString cinput = cinput_main_MC + "/" + sfname;
      foutput->cd();
      TChain* tin = new TChain("SkimTree");
      int nfiles = tin->Add(cinput);
      MELAout << "\t- Successfully added " << nfiles << " files for " << sname << " from " << cinput << "..." << endl;
      samples_all.emplace_back(sgroup, tin);
      tins_collected.push_back(tin);
      norm_map[tin] = 1;
      if (sname_sfname_pairs.size()>1){
        norm_map[tin] = 0;
        std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
        double sum_wgts = 0;
        bool hasCounters = true;
        {
          int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
          int bin_period = 1;
          for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
            if (validDataPeriods.at(iperiod)==SampleHelpers::getDataPeriod()){ bin_period += iperiod+1; break; }
          }
          for (auto const& fname:inputfilenames){
            TFile* ftmp = TFile::Open(fname, "read");
            TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
            if (!hCounters){
              hasCounters = false;
              sum_wgts = 0;
              break;
            }
            sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
            ftmp->Close();
            foutput->cd();
          }
          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
        }
        norm_map[tin] += sum_wgts;
      }
      foutput->cd();
    }
    {
      double sum_wgts_MC = 0;
      for (auto const& tin:tins_collected) sum_wgts_MC += norm_map[tin];
      for (auto const& tin:tins_collected) norm_map[tin] /= sum_wgts_MC;
    }
  }
  for (auto const& sgroup_tin_pair:samples_all) MELAout
    << "Relative normalization for sample in group " << sgroup_tin_pair.first << " = " << norm_map[sgroup_tin_pair.second]
    << endl;

  TString const cinput_main_data =
    "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + strdate
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY")
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  std::vector<TString> sfnames_data;
  getDataSampleDirs(sfnames_data);
  sgroups.push_back("Data");
  for (auto const& sfname:sfnames_data){
    TString cinput = cinput_main_data + "/" + sfname;
    foutput->cd();
    TChain* tin = new TChain("SkimTree");
    int nfiles = tin->Add(cinput);
    MELAout << "\t- Successfully added " << nfiles << " files for data from " << cinput << "..." << endl;
    samples_all.emplace_back("Data", tin);
    norm_map[tin] = 1;
    foutput->cd();
  }


  // Recorded variables
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(bool, leptons_is_TOmatched_SingleLepton) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(float, leptons_eff) \
  BRANCH_COMMAND(float, leptons_eff_DF) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched_fullCone) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_phi) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float> ME_Kfactor_values;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

    std::vector<TString> allbranchnames;
    {
      // Then check all leaves
      const TList* llist = (const TList*) tin->GetListOfLeaves();
      if (llist){
        for (int ib=0; ib<llist->GetSize(); ib++){
          auto const& bmem = llist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
         }
      }
      // Then check all branches
      const TList* blist = (const TList*) tin->GetListOfBranches();
      if (blist){
        for (int ib=0; ib<blist->GetSize(); ib++){
          auto const& bmem = blist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
        }
      }

      for (auto const& bname:allbranchnames){
        if (
          (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
          ||
          bname.BeginsWith("KFactor")
          ) ME_Kfactor_values[bname] = -1;
      }
    }

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (auto& it:ME_Kfactor_values){
      TString const& MEname = it.first;
      if (!HelperFunctions::checkListVariable(allbranchnames, MEname)) continue;
      float& MEval = it.second;
      tin->SetBranchStatus(MEname, 1); tin->SetBranchAddress(MEname, &MEval);
    }
  }

  // Build discriminants
  Discriminant* DjjVBF = DiscriminantClasses::constructKDFromType(
    DiscriminantClasses::kDjjVBF,
    ANALYSISTREEPKGDATAPATH+"RecoMEConstants/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth"
  );

  foutput->cd();

  // Create the histograms
  std::vector<TString> sgroups_allhists = sgroups;
  if (HelperFunctions::checkListVariable<TString>(sgroups, "DY_2l")){
    sgroups_allhists.push_back("DY_2l_genMatched");
    sgroups_allhists.push_back("DY_2l_failedMatches");
  }
  if (HelperFunctions::checkListVariable<TString>(sgroups, "qqZZ_2l2nu")){
    sgroups_allhists.push_back("qqZZ_2l2nu_genMatched");
    sgroups_allhists.push_back("qqZZ_2l2nu_failedMatches");
  }
  if (HelperFunctions::checkListVariable<TString>(sgroups, "qqWZ_3lnu")){
    sgroups_allhists.push_back("qqWZ_3lnu_genMatched");
    sgroups_allhists.push_back("qqWZ_3lnu_failedMatches");
  }
  sgroups_allhists.push_back("AllMC_NonRes");
  unsigned int const nprocesses = sgroups_allhists.size();

  std::vector<TH1F*> hlist;
  std::vector<ExtendedHistogram_1D*> hlist_SB;
  std::vector<ExtendedHistogram_1D*> hlist_SB_ttbar;
  for (auto const& sgroup:sgroups_allhists){
    int scolor = (int) kBlack;
    if (sgroup == "Data") scolor = (int) (kBlack);
    else if (sgroup == "AllMC_NonRes") scolor = (int) (kGray);
    else if (sgroup.Contains("DY_2l")) scolor = (int) (kCyan);
    else if (sgroup.Contains("qqZZ_2l2nu")) scolor = (int) (kYellow-3);
    else if (sgroup.Contains("qqWZ_3lnu")) scolor = (int) (kBlue);
    else if (sgroup == "TT_2l2nu") scolor = (int) (kOrange-3);
    else if (sgroup == "qqWW_2l2nu") scolor = (int) (kTeal-1);
    else if (sgroup == "WJets_lnu") scolor = (int) (kRed);
    else MELAerr << "Sample type " << sgroup << " is not recognized!" << endl;

    for (unsigned int ic=0; ic<nchannels; ic++){
      auto const& strChannelName = strChannelNames.at(ic);
      auto const& strChannelTitle = strChannelTitles.at(ic);
      for (unsigned int ij=0; ij<nbins_njets; ij++){
        auto const& strNjetsName = strNjetsNames.at(ij);
        auto const& strNjetsTitle = strNjetsTitles.at(ij);
        for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
          auto const& strBtaggingRegionName = strBtaggingRegionNames.at(ibt);
          auto const& strBtaggingRegionTitle = strBtaggingRegionTitles.at(ibt);

          TString strCutName = strChannelName + "_" + strNjetsName + "_" + strBtaggingRegionName;
          TString strCutTitle = strChannelTitle + "|" + strNjetsTitle + ", " + strBtaggingRegionTitle;

          ExtendedBinning binning_mTZZ(27, 150, 1500, "m_{T}^{ZZ} (GeV)");
          ExtendedBinning binning_mll(30, 91.2-15., 91.2+15., "m_{ll} (GeV)");
          ExtendedBinning binning_pTl1(19, 25, 500, "p_{T}^{l1} (GeV)");
          ExtendedBinning binning_pTl2(30, 25, 325, "p_{T}^{l2} (GeV)");
          ExtendedBinning binning_pTmiss(35, 125, 1000, "p_{T}^{miss} (GeV)");
          ExtendedBinning binning_pTmiss_vbin({ thr_corr_pTmiss, 110, 125, 150, 200, 300, 13000 }, "pTmiss", "p_{T}^{miss} (GeV)");
          ExtendedBinning binning_DjjVBF(10, 0, 1, "D_{2jet}^{VBF}");

          TH1F* htmp = nullptr;
#define HISTOGRAM_COMMAND(NAME, BINNING) \
          htmp = new TH1F(Form("h_%s_%s_%s", #NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING.getNbins(), BINNING.getBinning()); htmp->Sumw2(); \
          htmp->GetXaxis()->SetTitle(BINNING.getLabel()); htmp->GetYaxis()->SetTitle("Events / bin"); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist.push_back(htmp);

          if (istep>0){
            HISTOGRAM_COMMAND(mTZZ, binning_mTZZ);
            HISTOGRAM_COMMAND(mll, binning_mll);
            HISTOGRAM_COMMAND(pTl1, binning_pTl1);
            HISTOGRAM_COMMAND(pTl2, binning_pTl2);
            HISTOGRAM_COMMAND(DjjVBF, binning_DjjVBF);
          }
          HISTOGRAM_COMMAND(pTmiss, binning_pTmiss);

#undef HISTOGRAM_COMMAND

          ExtendedHistogram_1D* ehtmp = nullptr;
#define HISTOGRAM_COMMAND(NAME, BINNING) \
          ehtmp = new ExtendedHistogram_1D(Form("h_SB_%s_%s_%s", #NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING); \
          htmp = ehtmp->getHistogram(); \
          htmp->GetYaxis()->SetTitle("Events / bin"); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist_SB.push_back(ehtmp);

          HISTOGRAM_COMMAND(pTmiss, binning_pTmiss_vbin);

#undef HISTOGRAM_COMMAND

#define HISTOGRAM_COMMAND(NAME, BINNING) \
          ehtmp = new ExtendedHistogram_1D(Form("h_SB_ttbar_%s_%s_%s", #NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING); \
          htmp = ehtmp->getHistogram(); \
          htmp->GetYaxis()->SetTitle("Events / bin"); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist_SB_ttbar.push_back(ehtmp);

          HISTOGRAM_COMMAND(pTmiss, binning_pTmiss_vbin);

#undef HISTOGRAM_COMMAND
        }
      }
    }
  }

  using namespace OffshellCutflow;

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);
  for (auto const& spair:samples_all){
    auto const& sgroup = spair.first;
    auto const& tin = spair.second;

    // Reset ME and K factor values
    for (auto& it:ME_Kfactor_values) it.second = -1;
    bool is_qqVV = sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW");
    bool is_ggVV = sgroup.Contains("ggZZ") || sgroup.Contains("ggWW") || sgroup.Contains("GGH");
    bool isData = (sgroup == "Data");

    float* val_Kfactor_QCD = nullptr;
    float* val_Kfactor_EW = nullptr;
    if (is_qqVV){
      val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_qqVV_Bkg_Nominal")->second);
      switch (theGlobalSyst){
      case tEWDn:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWDn")->second);
        break;
      case tEWUp:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWUp")->second);
        break;
      default:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_Nominal")->second);
        break;
      }
    }
    if (is_ggVV){
      switch (theGlobalSyst){
      case tQCDScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn")->second);
        break;
      case tQCDScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp")->second);
        break;
      case tPDFScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn")->second);
        break;
      case tPDFScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp")->second);
        break;
      case tPDFReplicaDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn")->second);
        break;
      case tPDFReplicaUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp")->second);
        break;
      case tAsMZDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsDn")->second);
        break;
      case tAsMZUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsUp")->second);
        break;
      default:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second);
        break;
      }
    }

    int const nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;

      if (event_pTmiss<thr_corr_pTmiss) continue;
      bool pass_SR_pTmiss = check_pTmiss(event_pTmiss);

      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss)) continue;
      if (!check_pTl1(dilepton_pt)) continue;

      bool hasGenMatchedPair = isData || (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back());
      float pTl1 = std::max(leptons_pt->front(), leptons_pt->back());
      float pTl2 = std::min(leptons_pt->front(), leptons_pt->back());
      if (!check_pTl1(pTl1)) continue;
      if (!check_pTl2(pTl2)) continue;

      bool is_emu = (dilepton_id==-143);
      bool is_mumu = (dilepton_id==-169);
      bool is_ee = (dilepton_id==-121);
      if (std::abs(dilepton_mass-MZ_VAL_CUTS)>=30.f) continue;
      bool pass_SR_mll = check_mll(dilepton_mass, is_ee || is_mumu);

      event_wgt_SFs_PUJetId = std::min(event_wgt_SFs_PUJetId, 3.f);
      float wgt = event_wgt * event_wgt_SFs_muons * event_wgt_SFs_electrons * event_wgt_SFs_photons * event_wgt_SFs_PUJetId * event_wgt_SFs_btagging * norm_map[tin];
      if (val_Kfactor_QCD) wgt *= *val_Kfactor_QCD;
      if (val_Kfactor_EW) wgt *= *val_Kfactor_EW;

      float wgt_emu_rewgt_ee = 1;
      float wgt_emu_rewgt_mumu = 1;
      if (is_emu){
        if (std::abs(leptons_id->front())==13) wgt_emu_rewgt_ee *= leptons_eff_DF->front() / leptons_eff->front();
        else wgt_emu_rewgt_mumu *= leptons_eff_DF->front() / leptons_eff->front();
        if (std::abs(leptons_id->back())==13) wgt_emu_rewgt_ee *= leptons_eff_DF->back() / leptons_eff->back();
        else wgt_emu_rewgt_mumu *= leptons_eff_DF->back() / leptons_eff->back();
      }

      float SFself=1, effself=1;
      triggerSFHandler.getCombinedDileptonSFAndEff(
        theGlobalSyst,
        leptons_pt->front(), leptons_eta->front(), leptons_id->front(),
        leptons_pt->back(), leptons_eta->back(), leptons_id->back(),
        true,
        SFself, &effself
      );
      if (!isData) wgt *= SFself;
      if (is_emu){
        float SF_ee=1, eff_ee=1;
        float SF_mumu=1, eff_mumu=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          theGlobalSyst,
          leptons_pt->front(), leptons_eta->front(), 11,
          leptons_pt->back(), leptons_eta->back(), -11,
          true,
          SF_ee, &eff_ee
        );
        triggerSFHandler.getCombinedDileptonSFAndEff(
          theGlobalSyst,
          leptons_pt->front(), leptons_eta->front(), 13,
          leptons_pt->back(), leptons_eta->back(), -13,
          true,
          SF_mumu, &eff_mumu
        );
        wgt_emu_rewgt_ee *= eff_ee/effself;
        wgt_emu_rewgt_mumu *= eff_mumu/effself;

        for (auto const& tg:tg_incorr_list[!isData][0][1*(event_n_ak4jets_pt30==1) + 2*(event_n_ak4jets_pt30>=2)]) wgt_emu_rewgt_mumu *= tg->Eval(event_pTmiss);
        for (auto const& tg:tg_incorr_list[!isData][1][1*(event_n_ak4jets_pt30==1) + 2*(event_n_ak4jets_pt30>=2)]) wgt_emu_rewgt_ee *= tg->Eval(event_pTmiss);
      }

      bool const pass_sel_SR = pass_SR_pTmiss && pass_SR_mll;
      bool const pass_sel_SB = !pass_SR_mll;
      bool const pass_sel_SB_ttbar = dilepton_mass>=MZ_VAL_CUTS+15.f;

      // Fill histograms
      if (pass_sel_SR){
        DjjVBF->update({ ME_Kfactor_values["p_JJVBF_SIG_ghv1_1_JHUGen"], ME_Kfactor_values["p_JJQCD_SIG_ghg2_1_JHUGen"] }, event_mZZ);
        for (auto& hh:hlist){
          TString hname = hh->GetName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    sgroup == "DY_2l"
                    ||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("mTZZ")) hh->Fill(event_mTZZ, hwgt);
          if (hname.Contains("mll")) hh->Fill(dilepton_mass, hwgt);
          if (hname.Contains("pTl1")) hh->Fill(pTl1, hwgt);
          if (hname.Contains("pTl2")) hh->Fill(pTl2, hwgt);
          if (hname.Contains("pTmiss")) hh->Fill(event_pTmiss, hwgt);
          if (hname.Contains("DjjVBF") && event_n_ak4jets_pt30>=2) hh->Fill(*DjjVBF, hwgt);
        }
      }

      if (pass_sel_SB){
        for (auto& hh:hlist_SB){
          TString hname = hh->getName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    sgroup == "DY_2l"
                    ||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("pTmiss")) hh->fill(event_pTmiss, hwgt);
        }
      }

      if (pass_sel_SB_ttbar){
        for (auto& hh:hlist_SB_ttbar){
          TString hname = hh->getName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    sgroup == "DY_2l"
                    ||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("pTmiss")) hh->fill(event_pTmiss, hwgt);
        }
      }

    }
  }

  delete DjjVBF;

  for (auto& hh:hlist){
    foutput->WriteTObject(hh);
    delete hh;
  }

  {
    unsigned int const nhists_SB = hlist_SB.size()/(nprocesses * nchannels * nbins_njets * nbins_nbtagged);
    std::vector<unsigned int> idxs_process;
    {
      unsigned int idx=0;
      for (auto const& sgroup:sgroups_allhists){
        if (
          sgroup == "Data"
          ||
          sgroup == "AllMC_NonRes"
          ||
          sgroup == "TT_2l2nu"
          ) idxs_process.push_back(idx);
        idx++;
      }
    }

    for (auto const& idx_process:idxs_process){
      for (unsigned int ic=0; ic<2; ic++){
        for (unsigned int ij=0; ij<nbins_njets; ij++){
          for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
            for (unsigned int ih=0; ih<nhists_SB; ih++){
              auto const& eh_num = hlist_SB.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*(ic+nchannels*(idx_process)))));
              auto const& eh_den = hlist_SB.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*((ic+3)+nchannels*(idx_process)))));

              TString tgname = eh_num->getName();
              if (!tgname.Contains("pTmiss")) continue;
              HelperFunctions::replaceString<TString>(tgname, "h_SB", "tg_corr");

              TGraphErrors* tgnum = eh_num->getGraph(tgname+"_num");
              TGraphErrors* tgdenom = eh_den->getGraph(tgname+"_denom");

              unsigned int const npoints = tgnum->GetN();
              assert(npoints == tgdenom->GetN());

              std::vector<std::pair<double, double>> xy; xy.reserve(npoints);
              std::vector<std::pair<double, double>> xy_dn; xy_dn.reserve(npoints);
              std::vector<std::pair<double, double>> xy_up; xy_up.reserve(npoints);
              double* xxnum = tgnum->GetX();
              double* exnum = tgnum->GetEX();
              double* yynum = tgnum->GetY();
              double* eynum = tgnum->GetEY();
              double* xxdenom = tgdenom->GetX();
              double* exdenom = tgdenom->GetEX();
              double* yydenom = tgdenom->GetY();
              double* eydenom = tgdenom->GetEY();
              double xlast=0;
              for (unsigned int ip=0; ip<npoints; ip++){
                double cxx = (xxnum[ip]/std::pow(exnum[ip], 2) + xxdenom[ip]/std::pow(exdenom[ip], 2)) / (1./std::pow(exnum[ip], 2) + 1./std::pow(exdenom[ip], 2));

                double cyy = yynum[ip] / yydenom[ip];

                double cyy_dn = cyy - std::sqrt(std::pow((yynum[ip] - eynum[ip])/yydenom[ip] - cyy, 2) + std::pow(yynum[ip]/(yydenom[ip] + eydenom[ip]) - cyy, 2));
                double cyy_up = cyy + std::sqrt(std::pow((yynum[ip] + eynum[ip])/yydenom[ip] - cyy, 2) + std::pow(yynum[ip]/(yydenom[ip] - eydenom[ip]) - cyy, 2));

                xy.emplace_back(cxx, cyy);
                xy_dn.emplace_back(cxx, cyy_dn);
                xy_up.emplace_back(cxx, cyy_up);

                if (ip==npoints-1) xlast = cxx;
              }

              delete tgnum;
              delete tgdenom;

              TGraph* tg_tmp = nullptr;
              tg_tmp = HelperFunctions::makeGraphFromPair(xy, tgname+"_Nominal"); foutput->WriteTObject(tg_tmp); delete tg_tmp;
              tg_tmp = HelperFunctions::makeGraphFromPair(xy_dn, tgname+"_StatDn"); foutput->WriteTObject(tg_tmp); delete tg_tmp;
              tg_tmp = HelperFunctions::makeGraphFromPair(xy_up, tgname+"_StatUp"); foutput->WriteTObject(tg_tmp); delete tg_tmp;
            }
          }
        }
      }
    }
  }

  for (auto& hh:hlist_SB){
    foutput->WriteTObject(hh->getHistogram());
    delete hh;
  }
  for (auto& hh:hlist_SB_ttbar){
    foutput->WriteTObject(hh->getHistogram());
    delete hh;
  }

  for (auto& pp:samples_all) delete pp.second;

  for (auto& ftmp:finputs_prevstep) ftmp->Close();

  foutput->Close();
  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);

#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS
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
  sampleList.emplace_back("DY_2l_M_50_ext", "DY ll", "DY_2l_M_50_ext", -1, HistogramProperties((int) kCyan, 1, 2));
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
