#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "ParticleObjectHelpers.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "utils.hh"

using namespace reco;

inline void moveOverFlowToLastBin1D(TH1* hist) {
  int nbin = hist->GetNbinsX();
  if (hist->GetBinLowEdge(nbin+1) < 1499.9) return;  // only when the last bin is the infinity bin
  if (hist->GetBinContent(nbin+1) > 0) {
    // cout << "Moving the overflow for hist: " << hist->GetTitle() << " to its last bin!" << endl;
    double err = 0;
    hist->SetBinContent(nbin, hist->IntegralAndError(nbin, -1, err));
    hist->SetBinError(nbin, err);
    hist->SetBinContent(nbin+1, 0);
    hist->SetBinError(nbin+1, 0);
  }
}

template<typename... TArgs>
void plot1d(std::string name, double xval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH1D* currentHisto= new TH1D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    iter->second->Fill(xval, weight);
  }
}

inline bool isCloseObject(const float eta1, const float phi1, const float eta2, const float phi2, const float conesize = 0.4, float* deltaR = nullptr) {
  const float PI = TMath::Pi();
  float deltaEta = fabs(eta1 - eta2);
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(phi1 - phi2);
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;
  if (deltaR) *deltaR = sqrt(deltaR2);

  return true;
}

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

bool check_VBF_category(std::vector<AK4JetObject*> const& ak4jets_tight, ParticleObject::LorentzVector_t const& boson){
  if (ak4jets_tight.size()<2) return false;
  auto itFirstJet = ak4jets_tight.cbegin();
  AK4JetObject* leading_jet = *itFirstJet; itFirstJet++;
  AK4JetObject* subleading_jet = *itFirstJet; itFirstJet++;
  float leading_eta = leading_jet->eta();
  float subleading_eta = subleading_jet->eta();
  if (leading_eta<subleading_eta) std::swap(leading_eta, subleading_eta);
  if (std::abs(leading_eta - subleading_eta)<4.f) return false;
  if ((leading_jet->p4() + subleading_jet->p4()).M()<500.f) return false;
  float eta_cand = boson.Eta();
  if (eta_cand<=subleading_eta || eta_cand>=leading_eta) return false;
  for (auto it=itFirstJet; it!=ak4jets_tight.cend(); it++){
    float eta_jet = (*it)->eta();
    if (eta_jet>subleading_eta && eta_jet<leading_eta) return false;
  }
  return true;
}

const bool applyPUIdToAK4Jets = true;
const bool applyTightLeptonVetoIdToAK4Jets = false;

void getTrees(int procsel, int ichunk, int nchunks, TString strdate, TString year="2018") {
  if (procsel<0) return;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main = "output/DYestTest/SkimTrees/" + strdate + "_" + year;
  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure(year, "hadoop:200420_"+year);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  // std::vector<std::string> triggerCheckList_OSDF = OffshellTriggerHelpers::getHLTMenus({
  //     OffshellTriggerHelpers::kMuEle,
  //     OffshellTriggerHelpers::kSingleEle, OffshellTriggerHelpers::kSingleMu
  //   });
  // std::vector<std::string> triggerCheckList_OSSF = OffshellTriggerHelpers::getHLTMenus({
  //     OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kSingleEle,
  //     OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kSingleMu
  //   });

  // std::vector<std::string> triggerCheckList_muon = OffshellTriggerHelpers::getHLTMenus({
  //       OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kSingleMu
  //   });
  // std::vector<std::string> triggerCheckList_electron = OffshellTriggerHelpers::getHLTMenus({
  //       OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kSingleEle,
  //   });
  // std::vector<std::string> triggerCheckList_photon = OffshellTriggerHelpers::getHLTMenus({
  //     OffshellTriggerHelpers::kSinglePho,
  //   });

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_2l{
    // Main unprescaled triggers for signal regions
    TriggerHelpers::kDoubleMu, TriggerHelpers::kDoubleEle, TriggerHelpers::kMuEle,
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_el{
     TriggerHelpers::kDoubleEle, TriggerHelpers::kSingleEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_mu{
    TriggerHelpers::kSingleMu, TriggerHelpers::kDoubleMu
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_emu{ TriggerHelpers::kMuEle };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_ph{ TriggerHelpers::kSinglePho };

  std::vector<std::string> triggerCheckList_2l = TriggerHelpers::getHLTMenus(requiredTriggers_2l);
  std::vector<std::string> triggerCheckList_el = TriggerHelpers::getHLTMenus(requiredTriggers_el);
  std::vector<std::string> triggerCheckList_mu = TriggerHelpers::getHLTMenus(requiredTriggers_mu);
  std::vector<std::string> triggerCheckList_emu = TriggerHelpers::getHLTMenus(requiredTriggers_emu);
  // auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);
  auto triggerPropsCheckList_ph = TriggerHelpers::getHLTMenuProperties(requiredTriggers_ph);

  std::vector<SampleSpecs> sampleList;
  // ------------------------------- MC to be run over photon skims ---------------------------------------- //
  sampleList.emplace_back("GJets_HT_40-100",  "#gamma+jets (H_{T}-binned)", "GJets_HT_40-100",  -1, HistogramProperties((int) kCyan, 1, 2)); // 0
  sampleList.emplace_back("GJets_HT_100-200", "#gamma+jets (H_{T}-binned)", "GJets_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2)); // 1
  sampleList.emplace_back("GJets_HT_200-400", "#gamma+jets (H_{T}-binned)", "GJets_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2)); // 2
  sampleList.emplace_back("GJets_HT_400-600", "#gamma+jets (H_{T}-binned)", "GJets_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2)); // 3
  sampleList.emplace_back("GJets_HT_600-inf", "#gamma+jets (H_{T}-binned)", "GJets_HT_600-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 4
  sampleList.emplace_back("QCD_HT", "QCD HT binned", "QCD_HT", -1, HistogramProperties((int) kCyan, 1, 2)); // 5
  sampleList.emplace_back("ZJets_nunu", "ZJets Znunu HT binned", "ZJets_nunu", -1, HistogramProperties((int) kCyan, 1, 2)); // 6
  sampleList.emplace_back("ZGJets_nunu_pTG_40-130", "ZGJets Znunu", "ZGJets_nunu_pTG_40-130", -1, HistogramProperties((int) kCyan, 1, 2)); // 7
  sampleList.emplace_back("ZGJets_nunu_pTG_130-inf", "ZGJets Znunu", "ZGJets_nunu_pTG_130-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 8
  sampleList.emplace_back("ZGJets_nunu_nlo_incl", "ZGJets Znunu", "ZGJets_nunu_nlo_inclusive", -1, HistogramProperties((int) kCyan, 1, 2)); // 9
  sampleList.emplace_back("ZGJets_nunu_nlo_pTG_130-inf", "ZGJets Znunu", "ZGJets_nunu_nlo_pTG_130-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 10
  sampleList.emplace_back("qqWG_lnu", "WGToLNuG", "qqWG_lnu", -1, HistogramProperties((int) kCyan, 1, 2)); // 11
  sampleList.emplace_back("WZG", "WZG", "WZG", -1, HistogramProperties((int) kCyan, 1, 2)); // 12
  sampleList.emplace_back("TTGJets", "TTGJets", "TTGJets", -1, HistogramProperties((int) kCyan, 1, 2)); // 13
  sampleList.emplace_back("TGJets", "TTGJets", "TGJets", -1, HistogramProperties((int) kCyan, 1, 2)); // 14
  int phskim_size = sampleList.size();    /// ----------- line between gjets skims and 2l skims by default ------------------- //
  sampleList.emplace_back("ZGJets_ll_pTG_40-130", "ZGJets Zll", "ZGJets_ll_pTG_40-130", -1, HistogramProperties((int) kCyan, 1, 2)); // 15
  sampleList.emplace_back("ZGJets_ll_pTG_130-inf", "ZGJets Zll", "ZGJets_ll_pTG_130-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 16
  sampleList.emplace_back("ZGJets_ll_nlo_incl", "ZGJets Zll", "ZGJets_ll_nlo_inclusive", -1, HistogramProperties((int) kCyan, 1, 2)); // 17
  sampleList.emplace_back("ZGJets_ll_nlo_pTG_130-inf", "ZGJets Zll", "ZGJets_ll_nlo_pTG_130-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 18
  sampleList.emplace_back("DY_2l_M_50_ext", "DY ll (NLO)", "DY_2l_M_50_ext", -1, HistogramProperties((int) kCyan, 1, 2)); // 19
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 20
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 21
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 22
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2)); // 23
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 24
  sampleList.emplace_back("TTZ_2l2nu", "TTZ_2l2nu", "TTZ_2l2nu", -1, HistogramProperties((int) kCyan, 1, 2)); // 25
  sampleList.emplace_back("TTW_lnu", "TTW_lnu", "TTW_lnu", -1, HistogramProperties((int) kCyan, 1, 2)); // 26
  sampleList.emplace_back("WJets_lnu_HT", "WJets_lnu_HT binned", "WJets_lnu_HT", -1, HistogramProperties((int) kCyan, 1, 2)); // 27
  sampleList.emplace_back("ST_t-channel_top_5f", "ST_t-channel_top_5f", "ST_t-channel_top_5f", -1, HistogramProperties((int) kCyan, 1, 2)); // 28
  sampleList.emplace_back("ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f", "ST_t-channel_antitop_5f", -1, HistogramProperties((int) kCyan, 1, 2)); // 29
  sampleList.emplace_back("ST_s-channel_top_leptonDecays", "ST_s-channel_top_leptonDecays", "ST_s-channel_top_leptonDecays", -1, HistogramProperties((int) kCyan, 1, 2)); // 30
  sampleList.emplace_back("ST_s-channel_antitop_leptonDecays", "ST_s-channel_antitop_leptonDecays", "ST_s-channel_antitop_leptonDecays", -1, HistogramProperties((int) kCyan, 1, 2)); // 31
  sampleList.emplace_back("ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays", -1, HistogramProperties((int) kCyan, 1, 2)); // 32
  sampleList.emplace_back("ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_antitop_5f_NoFullyHadronicDecays", -1, HistogramProperties((int) kCyan, 1, 2)); // 33
  sampleList.emplace_back("DY_2l_M_50_HT_70-100", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_70-100", -1, HistogramProperties((int) kCyan, 1, 2)); // 34
  sampleList.emplace_back("DY_2l_M_50_HT_100-200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2)); // 35
  sampleList.emplace_back("DY_2l_M_50_HT_200-400", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2)); // 36
  sampleList.emplace_back("DY_2l_M_50_HT_400-600", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2)); // 37
  sampleList.emplace_back("DY_2l_M_50_HT_600-800", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_600-800", -1, HistogramProperties((int) kCyan, 1, 2)); // 38
  sampleList.emplace_back("DY_2l_M_50_HT_800-1200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_800-1200", -1, HistogramProperties((int) kCyan, 1, 2)); // 39
  sampleList.emplace_back("DY_2l_M_50_HT_1200-2500", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_1200-2500", -1, HistogramProperties((int) kCyan, 1, 2)); // 40
  sampleList.emplace_back("DY_2l_M_50_HT_2500-inf", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_2500-inf", -1, HistogramProperties((int) kCyan, 1, 2)); // 41
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 42
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 43
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 44
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 45
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 46
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 47
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 48
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 49
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 50
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 51
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 52
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2)); // 53

  // ------------------------------- Data photon skim ---------------------------------------- //
  sampleList.emplace_back("2018A_phskim", "Observed (2018A)", "Run2018A", -1, HistogramProperties((int) kCyan, 1, 2)); // 54
  sampleList.emplace_back("2018B_phskim", "Observed (2018B)", "Run2018B", -1, HistogramProperties((int) kCyan, 1, 2)); // 55
  sampleList.emplace_back("2018C_phskim", "Observed (2018C)", "Run2018C", -1, HistogramProperties((int) kCyan, 1, 2)); // 56
  sampleList.emplace_back("2018D_phskim", "Observed (2018D)", "Run2018D", -1, HistogramProperties((int) kCyan, 1, 2)); // 57

  sampleList.emplace_back("2017B_phskim", "Observed (2017B)", "Run2017B", -1, HistogramProperties((int) kCyan, 1, 2)); // 58
  sampleList.emplace_back("2017C_phskim", "Observed (2017C)", "Run2017C", -1, HistogramProperties((int) kCyan, 1, 2)); // 59
  sampleList.emplace_back("2017D_phskim", "Observed (2017D)", "Run2017D", -1, HistogramProperties((int) kCyan, 1, 2)); // 60
  sampleList.emplace_back("2017E_phskim", "Observed (2017E)", "Run2017E", -1, HistogramProperties((int) kCyan, 1, 2)); // 61
  sampleList.emplace_back("2017F_phskim", "Observed (2017F)", "Run2017F", -1, HistogramProperties((int) kCyan, 1, 2)); // 62

  sampleList.emplace_back("2016B_phskim", "Observed (2016B)", "Run2016B", -1, HistogramProperties((int) kCyan, 1, 2)); // 63
  sampleList.emplace_back("2016C_phskim", "Observed (2016C)", "Run2016C", -1, HistogramProperties((int) kCyan, 1, 2)); // 64
  sampleList.emplace_back("2016D_phskim", "Observed (2016D)", "Run2016D", -1, HistogramProperties((int) kCyan, 1, 2)); // 65
  sampleList.emplace_back("2016E_phskim", "Observed (2016E)", "Run2016E", -1, HistogramProperties((int) kCyan, 1, 2)); // 66
  sampleList.emplace_back("2016F_phskim", "Observed (2016F)", "Run2016F", -1, HistogramProperties((int) kCyan, 1, 2)); // 67
  sampleList.emplace_back("2016G_phskim", "Observed (2016G)", "Run2016G", -1, HistogramProperties((int) kCyan, 1, 2)); // 68
  sampleList.emplace_back("2016H_phskim", "Observed (2016H)", "Run2016H", -1, HistogramProperties((int) kCyan, 1, 2)); // 69
  // ------------------------------- Data dilep skim ---------------------------------------- //
  sampleList.emplace_back("2018A_llskim", "Observed (2018A)", "Run2018A", -1, HistogramProperties((int) kCyan, 1, 2)); // 70
  sampleList.emplace_back("2018B_llskim", "Observed (2018B)", "Run2018B", -1, HistogramProperties((int) kCyan, 1, 2)); // 71
  sampleList.emplace_back("2018C_llskim", "Observed (2018C)", "Run2018C", -1, HistogramProperties((int) kCyan, 1, 2)); // 72
  sampleList.emplace_back("2018D_llskim", "Observed (2018D)", "Run2018D", -1, HistogramProperties((int) kCyan, 1, 2)); // 73

  sampleList.emplace_back("2017B_llskim", "Observed (2017B)", "Run2017B", -1, HistogramProperties((int) kCyan, 1, 2)); // 74
  sampleList.emplace_back("2017C_llskim", "Observed (2017C)", "Run2017C", -1, HistogramProperties((int) kCyan, 1, 2)); // 75
  sampleList.emplace_back("2017D_llskim", "Observed (2017D)", "Run2017D", -1, HistogramProperties((int) kCyan, 1, 2)); // 76
  sampleList.emplace_back("2017E_llskim", "Observed (2017E)", "Run2017E", -1, HistogramProperties((int) kCyan, 1, 2)); // 77
  sampleList.emplace_back("2017F_llskim", "Observed (2017F)", "Run2017F", -1, HistogramProperties((int) kCyan, 1, 2)); // 78

  sampleList.emplace_back("2016B_llskim", "Observed (2016B)", "Run2016B", -1, HistogramProperties((int) kCyan, 1, 2)); // 79
  sampleList.emplace_back("2016C_llskim", "Observed (2016C)", "Run2016C", -1, HistogramProperties((int) kCyan, 1, 2)); // 80
  sampleList.emplace_back("2016D_llskim", "Observed (2016D)", "Run2016D", -1, HistogramProperties((int) kCyan, 1, 2)); // 81
  sampleList.emplace_back("2016E_llskim", "Observed (2016E)", "Run2016E", -1, HistogramProperties((int) kCyan, 1, 2)); // 82
  sampleList.emplace_back("2016F_llskim", "Observed (2016F)", "Run2016F", -1, HistogramProperties((int) kCyan, 1, 2)); // 83
  sampleList.emplace_back("2016G_llskim", "Observed (2016G)", "Run2016G", -1, HistogramProperties((int) kCyan, 1, 2)); // 84
  sampleList.emplace_back("2016H_llskim", "Observed (2016H)", "Run2016H", -1, HistogramProperties((int) kCyan, 1, 2)); // 85

  sampleList.emplace_back("singleph_2016E_phskim", "Observed (2016C)", "SinglePhoton_2016E", -1, HistogramProperties((int) kCyan, 1, 2)); // 85


  //sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  if ((unsigned int) procsel>=sampleList.size()) return;
  auto& sample = sampleList.at(procsel);
  TString treename = "cms3ntuple/Events";
  if (
    sample.name.find("Sig") != std::string::npos ||
    sample.name.find("Bkg") != std::string::npos ||
    sample.name.find("BSI") != std::string::npos || 
    strdate.Contains("fullsamp")
   ) {
    SampleHelpers::configure(year, "hadoop:200420_"+year);
    SampleHelpers::setInputDirectory("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Production/");
  } else {
    SampleHelpers::setInputDirectory("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/");
    if (strdate.Contains("llskim") || strdate.Contains("skim2l"))
      treename = "cms3ntuple/Dilepton";
    else if (strdate.Contains("phskim") || strdate.Contains("skimph"))
      treename = "cms3ntuple/SinglePhoton";
    else if (procsel < phskim_size || sample.name.find("_phskim") != std::string::npos)
      treename = "cms3ntuple/SinglePhoton";
    else 
      treename = "cms3ntuple/Dilepton";
    if (year == "2016") {
      if (sample.name == "DY_2l_M_50_ext")
        sample.name = sample.path = "DY_2l_M_50";
    }
  }
  MELAout << "Running over tree " << treename << endl;

  map<string,TH1*> hvec;  // test histograms

  // Get sample specifications
  TString const& strSampleSet = sample.path;
  bool const isData = SampleHelpers::checkSampleIsData(strSampleSet);
  if (isData && theGlobalSyst!=SystematicsHelpers::sNominal) return;
  // if (isData) return;

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
  PhotonScaleFactorHandler photonSFHandler;
  BtagScaleFactorHandler btagSFHandler;

  METCorrectionHandler metCorrectionHandler;

  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(true);

  bool applyMETResCorr = true;
  // if (applyMETResCorr) metCorrector.setup();

  TString stroutput = Form("%s/%s_%s_%s", coutput_main.Data(), sample.name.data(), year.Data(), SystematicsHelpers::getSystName(theGlobalSyst).data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
  stroutput += ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  TTree* tout = new TTree("SkimTree", "");
#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
  // Event variables
  BRANCH_COMMAND(float, event_wgt);
  BRANCH_COMMAND(float, event_wgt_SFs);
  BRANCH_COMMAND(float, pTmiss);
  BRANCH_COMMAND(float, phimiss);
  BRANCH_COMMAND(float, mTZZ);
  BRANCH_COMMAND(float, mZZ);
  BRANCH_COMMAND(unsigned, event_Njets);
  BRANCH_COMMAND(unsigned, event_Njets_btagged);
  BRANCH_COMMAND(unsigned, event_Njets_btagged_medium);
  BRANCH_COMMAND(unsigned, event_Njets20);
  BRANCH_COMMAND(unsigned, event_Njets20_btagged_loose);
  BRANCH_COMMAND(unsigned, event_Njets20_btagged_medium);
  BRANCH_COMMAND(unsigned, event_nvtxs_good);
  BRANCH_COMMAND(unsigned, event_Nphotons);
  BRANCH_COMMAND(unsigned, event_HT);
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers);
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers);
  BRANCH_COMMAND(float, event_wgt_trig_muon);
  BRANCH_COMMAND(float, event_wgt_trig_electron);
  BRANCH_COMMAND(float, event_wgt_trig_photon);
  BRANCH_COMMAND(float, event_DjjVBF);
  BRANCH_COMMAND(float, event_DjjVBF_rl);
  BRANCH_COMMAND(float, met_uncorr_pt);
  BRANCH_COMMAND(float, met_uncorr_phi);
  BRANCH_COMMAND(float, genmet_pt);
  BRANCH_COMMAND(float, genmet_phi);

  BRANCH_COMMAND(bool, is_ee);
  BRANCH_COMMAND(bool, is_mumu);
  BRANCH_COMMAND(bool, is_emu);
  BRANCH_COMMAND(bool, is_gamma);
  BRANCH_COMMAND(bool, is_VBFcat);
  BRANCH_COMMAND(bool, has_electrons_inHEM1516);
  BRANCH_COMMAND(bool, has_photons_inHEM1516);
  BRANCH_COMMAND(bool, has_ak4jets_inHEM1516);

  BRANCH_COMMAND(float, pt_boson);
  BRANCH_COMMAND(float, eta_boson);
  BRANCH_COMMAND(float, phi_boson);
  BRANCH_COMMAND(float, mass_boson);

  BRANCH_COMMAND(bool, pass_lepveto);
  BRANCH_COMMAND(bool, pass_trackveto);
  BRANCH_COMMAND(int, id_track);
  BRANCH_COMMAND(float, pt_track);
  BRANCH_COMMAND(float, eta_track);
  BRANCH_COMMAND(float, phi_track);
  // L1, L2
  BRANCH_COMMAND(int, id_l1);
  BRANCH_COMMAND(float, pt_l1);
  BRANCH_COMMAND(float, eta_l1);
  BRANCH_COMMAND(float, phi_l1);
  BRANCH_COMMAND(int, id_l2);
  BRANCH_COMMAND(float, pt_l2);
  BRANCH_COMMAND(float, eta_l2);
  BRANCH_COMMAND(float, phi_l2);

  BRANCH_COMMAND(bool, convveto_photon);
  BRANCH_COMMAND(bool, id_photon_Hgg);
  BRANCH_COMMAND(bool, isEB_photon);
  BRANCH_COMMAND(float, pt_photon);
  BRANCH_COMMAND(float, eta_photon);
  BRANCH_COMMAND(float, phi_photon);

  BRANCH_COMMAND(float, pt_jet1);
  BRANCH_COMMAND(float, eta_jet1);
  BRANCH_COMMAND(float, phi_jet1);
  BRANCH_COMMAND(float, pt_jet2);
  BRANCH_COMMAND(float, eta_jet2);
  BRANCH_COMMAND(float, phi_jet2);

  BRANCH_COMMAND(float, mindphi_jet_met);
  BRANCH_COMMAND(float, dphi_boson_met);
  BRANCH_COMMAND(float, dphi_lljets_met);
  BRANCH_COMMAND(float, dphi_jet20_met);
  BRANCH_COMMAND(float, dphi_lljets20_met);
#undef BRANCH_COMMAND


  bool isFirstInputFile=true;
  for (auto const& sname : sampledirs) {
    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), treename, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);

    // Set data tracking options
    eventFilter.setTrackDataEvents(isData);
    eventFilter.setCheckUniqueDataEvent(isData && !isFirstInputFile);

    bool is_QCD = (sample.name.find("QCD_HT") == 0);
    bool is_ZGJets = (sample.name.find("ZGJets_nunu") == 0);
    bool is_ZGJets_LO = (sample.name.find("ZGJets_nunu") == 0 && sample.name.find("nlo") == string::npos);
    bool is_ZllGJets = (sample.name.find("ZGJets_ll") == 0);
    bool is_ZJets = (sample.name.find("ZJets_nunu") == 0);
    if (is_QCD || is_ZGJets || is_ZJets) genInfoHandler.setAcquireGenParticles(true);

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
    MELAout << "getHistograms: DjjVBF is built!" << endl;

    const int nEntries = sample_tree.getSelectedNEvents();

    double sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;

    // float sum_evts = 0.f;  // sum of nevents, weighted by -1 for negative weighted ones
    // float sum_ewgt = 0.f;  // sum of event weight, as a final validation

    unsigned* ptr_nvtx_good=nullptr;
    float* ptr_genmet_pt=nullptr;
    float* ptr_genmet_phi=nullptr;

    sample_tree.bookBranch<unsigned>("vtxs_nvtxs_good", 0);
    sample_tree.getValRef("vtxs_nvtxs_good", ptr_nvtx_good);

    if (!isData) {

      sample_tree.bookBranch<float>("xsec", 0.f);
      sample_tree.bookBranch<float>("genmet_met", 0);
      sample_tree.bookBranch<float>("genmet_metPhi", 0);

      sample_tree.getValRef("genmet_met", ptr_genmet_pt);
      sample_tree.getValRef("genmet_metPhi", ptr_genmet_phi);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      TString firstFile = SampleHelpers::getDatasetFirstFileName(sname);
      TFile* fFirstFile = TFile::Open(firstFile, "read");
      TH2D* hCounters = (TH2D*) fFirstFile->Get("cms3ntuple/Counters");
      if (hCounters) {
        // sum_wgts = hCounters->GetBinContent(1, 1);
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            MELAout << "\t- Failed to find the counters histogram in " << fname << endl;
            sum_wgts = getSumWeights(year, sname);
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          ftmp->Close();
        }
        // sum_wgts = getSumWeights(year, sname);
        MELAout << "Determined the sum of weights for this sample to be " << sum_wgts << endl;
      } else {
        if (true) {
          sum_wgts = getSumWeights(year, sname);
          MELAout << "Determined the sum of weights for this sample to be " << sum_wgts << endl;
        } else {
          MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
          for (int ev=0; ev<nEntries; ev++) {
            HelperFunctions::progressbar(ev, nEntries);
            sample_tree.getSelectedEvent(ev);

            genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
            auto const& genInfo = genInfoHandler.getGenInfo();
            double genwgt = genInfo->getGenWeight(true);

            simEventHandler.constructSimEvent(theGlobalSyst);
            double puwgt = simEventHandler.getPileUpWeight();
            sum_wgts += genwgt * puwgt;
          }
          MELAout << "Determined the sum of weights for this sample to be " << sum_wgts << " from the list." << endl;
        }
      }
      fFirstFile->Close();

    }
    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;

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
    if (nchunks>0) {
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Running over sample " << sample.name << endl;
    MELAout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    unsigned n_evts_acc=0;
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++) {
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);
      // if (ev>1'000'000) break;
      // if (ev>10) break;

      if (!isData && firstEvent) {
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        xsec *= 1000.;
      }
      if (firstEvent) firstEvent=false;

      event_wgt = (isData ? 1.f : xsec * lumi) / sum_wgts;
      
      eventFilter.constructFilters();

      int istep=0;
      float weight = event_wgt;

      string lepcat = "lepcat";

      auto fill_passedsteps = [&](string s="", bool extra=false) {
        if (extra) {
          plot1d("h_passed_ossteps_"+lepcat+s, istep , weight, hvec, ";OS selections" , 20,  0, 20);
        }
        plot1d("h_passed_ossteps"+s, istep++ , weight, hvec, ";OS selections" , 20,  0, 20);
      };
      fill_passedsteps();

      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;
      fill_passedsteps();

      event_nvtxs_good = *ptr_nvtx_good;
        
      if (!isData) {
        genmet_pt = *ptr_genmet_pt;
        genmet_phi = *ptr_genmet_phi;

        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();
        double genwgt = genInfo->getGenWeight(true);
        if (genwgt == -1.0 || genwgt == 0.) continue;
        event_wgt *= genwgt;
        fill_passedsteps();

        if (sample.name == "ggZZ_2l2nu_BSI" || sample.name == "ggWW_2l2nu_BSI") { event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (sample.name == "ggZZ_2l2nu_Sig" || sample.name == "ggWW_2l2nu_Sig") { event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (sample.name == "ggZZ_2l2nu_BSI_g4" || sample.name == "ggWW_2l2nu_BSI_g4") { event_wgt *= (genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz4_1_MCFM"]*2.5502 + genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*(-2.55052+std::pow(2.55052, 2)) + genInfo->extras.LHE_ME_weights["p_Gen_GG_BKG_MCFM"]*(-2.5502+1.))*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (sample.name == "ggZZ_2l2nu_Sig_g4" || sample.name == "ggWW_2l2nu_Sig_g4") { event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*std::pow(2.55052, 2)*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
        else if (sample.name == "VBFZZ_2l2nu_BSI" || sample.name == "VBFWW_2l2nu_BSI") { event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
        else if (sample.name == "VBFZZ_2l2nu_Sig" || sample.name == "VBFWW_2l2nu_Sig") { event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
        else if (sample.name == "VBFZZ_2l2nu_Sig_g4" || sample.name == "VBFWW_2l2nu_Sig_g4") { event_wgt *= (genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*std::pow(2.55052, 2) + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv4_1_MCFM"]*(-std::pow(2.55052, 2)+std::pow(2.55052, 4)) + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BKG_MCFM"]*(-std::pow(2.55052, 2)+1.))*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
        if (event_wgt==0.f) continue;
        fill_passedsteps();

        simEventHandler.constructSimEvent(theGlobalSyst);
        event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();
        if (event_wgt==0.f) continue;
        fill_passedsteps();
      }

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      float SF_muons = 1;
      float SF_electrons = 1;
      float SF_photons = 1;
      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      auto const& photons = photonHandler.getProducts();
      auto selphotons = photons;
      selphotons.clear();

      // for event passed histo
      event_wgt_trig_photon = eventFilter.getTriggerWeight(triggerPropsCheckList_ph, nullptr, nullptr, &photons, nullptr, nullptr, nullptr);
      if (is_gamma && isData) weight *= event_wgt_trig_photon;
      fill_passedsteps();

      event_Nphotons = 0;
      int nTightPhoton = 0;
      PhotonObject* thePhoton = nullptr;
      pt_photon = eta_photon = phi_photon = -999.;
      for (auto ph : photons) {
        // if (!ParticleSelectionHelpers::isTightParticle(ph)) continue; // tight id
        event_Nphotons++;
        pt_photon = ph->pt();
        eta_photon = ph->eta();
        phi_photon = ph->phi();
        convveto_photon = ph->testSelectionBit(PhotonSelectionHelpers::kConversionSafe);
        id_photon_Hgg = ph->extras.id_cutBased_HGG_Bits;
        isEB_photon = ph->isEB();
        if (!ph->isEB()) continue;  // photons with EB only
        float SF_photon = 1;
        if (!isData) photonSFHandler.getIdIsoSFAndEff(theGlobalSyst, ph, SF_photon, nullptr);
        if (SF_photon == 0.f) continue;
        SF_photons *= SF_photon;

        nTightPhoton++;
        thePhoton = ph;
        selphotons.push_back(ph);
      };
      is_gamma = (nTightPhoton == 1);
      string possuf = "";
      if (is_gamma) {
        pt_photon = thePhoton->pt();
        eta_photon = thePhoton->eta();
        phi_photon = thePhoton->phi();
        
        possuf = (thePhoton->isEB())? "_barrel" : (thePhoton->isEE())? "_endcap" : "_gap";
      }

      auto fill_passedsteps2 = [&](string s="", bool extra=false) {
        if (is_gamma) {
          plot1d("h_photon_pt_step"+to_string(istep)+s+possuf, pt_photon, weight, hvec, ";pt_{ph} [GeV]" , 160,  0, 800);
          plot1d("h_photon_eta_step"+to_string(istep)+s, eta_photon, weight, hvec, ";eta_{ph} [GeV]" , 96, -2.4, 2.4);
          plot1d("h_photon_phi_step"+to_string(istep)+s+possuf, phi_photon, weight, hvec, ";phi_{ph} [GeV]" , 128, -3.2, 3.2);
        }
        fill_passedsteps("", extra);
      };
      fill_passedsteps2("_tw");

      // Special veto on any prompt (gen) photon with pT > 25 GeV for QCD
      if (is_QCD) {
        bool hasPromptPhoton = false;
        for (auto genpart : genInfoHandler.getGenParticles()) {
          if (!genpart->extras.isPromptFinalState) continue;
          if (genpart->pdgId() != 22) continue;
          plot1d("h_GenPhoton_pt", genpart->pt(), weight, hvec, ";p_{T}(gen-#gamma)", 60,  0, 300);
          if (genpart->pt() < 25) continue;
          hasPromptPhoton = true;
          break;
        }
        if (hasPromptPhoton) continue;
      } else if (is_ZGJets_LO) {
        // ParticleObject::LorentzVector_t gennunu_p4;
        // int ngennu = 0;
        // for (auto genpart : genInfoHandler.getGenParticles()) {
        //   if (!genpart->extras.isPromptFinalState) continue;
        //   int id = abs(genpart->pdgId());
        //   if (id != 12 && id != 14 && id != 16) continue;
        //   // No way to check that they come from Z...
        //   gennunu_p4 += genpart->p4();
        //   ngennu++;
        // }
        // float gennunu_pt = gennunu_p4.pt();
        // plot1d("h_nGenNu", ngennu, weight, hvec, ";N_{#gamma}(gen)" , 6,  0, 6);
        // plot1d("h_GenNuNu_pt", gennunu_pt, weight, hvec, ";p_{T}(gen#nu#nu)" , 120,  0, 600);
        // plot1d("h_diff_genmet_GenNuNu_pt", fabs(genmet_pt-gennunu_pt), weight, hvec, ";#Delta p_{T}(genmet, gen#nu#nu)", 80,  0, 400);

        const vector<float> ptcats               = {0., 100., 120., 160., 200., 240., 280., 320., 360., 400., 440., 800., 880., 920., 960., 7500.};
        const vector<float> scalesNLO_gennunu_pt = {1., 1.32, 1.50, 1.75, 1.90, 1.96, 2.03, 2.07, 2.13, 2.16, 2.20, 2.16, 2.13, 2.10, 2.05,};
        int icat = std::upper_bound(ptcats.begin(), ptcats.end(), std::min(genmet_pt, 960.f)) - ptcats.begin() - 1;
        event_wgt *= scalesNLO_gennunu_pt.at(icat);

      } else if (is_ZllGJets) {
        // to add ll into the met
      
      }

      if (is_ZJets || is_ZGJets || is_ZllGJets) {
        bool hasPromptPhoton = false;
        for (auto genpart : genInfoHandler.getGenParticles()) {
          if (genpart->pdgId() != 22) continue;
          plot1d("h_GenPhoton_all_pt", genpart->pt(), weight, hvec, ";p_{T}(gen-#gamma)", 60,  0, 300);
          if (!genpart->extras.isPromptFinalState) continue;
          plot1d("h_GenPhoton_prompt_pt", genpart->pt(), weight, hvec, ";p_{T}(gen-#gamma)", 60,  0, 300);
          hasPromptPhoton = true;
        }
      }

      plot1d("h_genmet_pt", genmet_pt, weight, hvec, ";genmet [GeV]" , 120,  0, 600);

      isotrackHandler.constructIsotracks(&muons, &electrons);
      pass_trackveto = true;
      for (auto const& isotrack:isotrackHandler.getProducts()) {
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)) {
          pass_trackveto = false;
          id_track = isotrack->pdgId();
          pt_track = isotrack->pt();
          eta_track = isotrack->eta();
          phi_track = isotrack->phi();
          break;
        }
      }
      if (!pass_trackveto) continue;
      fill_passedsteps2("_tv");

      int nLeptons = 0;
      for (auto mu : muons) {
        if (ParticleSelectionHelpers::isTightParticle(mu)) nLeptons++;
      }
      for (auto el : electrons) {
        if (ParticleSelectionHelpers::isTightParticle(el)) nLeptons++;
      }

      if (nLeptons != 0 && nLeptons != 2) continue;
      fill_passedsteps2("_lv");

      event_wgt_SFs = SF_muons*SF_electrons*SF_photons;
      // MELAout << "Pass line " << __LINE__ << endl;

      jetHandler.constructJetMET(theGlobalSyst, &muons, &electrons, &selphotons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      met_uncorr_pt = pfmet->pt(false, false, false);
      met_uncorr_phi = pfmet->phi(false, false, false);
      if (applyMETResCorr && !isData)
        metCorrectionHandler.applyCorrections(simEventHandler.getChosenDataPeriod(), genmet_pt, genmet_phi, pfmet, true);

      pTmiss = pfmet->pt(true, true, true);
      phimiss = pfmet->phi(true, true, true);

      has_electrons_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, &electrons, nullptr, nullptr, nullptr);
      has_photons_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, nullptr, &photons, nullptr, nullptr);
      has_ak4jets_inHEM1516 = !eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, nullptr);

      event_HT = 0;
      float SF_btag = 1;
      ParticleObject::LorentzVector_t ak4jets_sump4;
      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      event_Njets_btagged = event_Njets_btagged_medium = event_Njets20_btagged_loose = event_Njets20_btagged_medium = event_Njets = event_Njets20 = 0;
      ParticleObject::LorentzVector_t ak4jets20_sump4;
      for (auto* jet:ak4jets) {
        // section for jet > 30 GeV
        if (ParticleSelectionHelpers::isTightJet(jet)) {
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) event_Njets_btagged++;
          if (jet->getBtagValue()>=btag_medium_thr) event_Njets_btagged_medium++;
          ak4jets_sump4 += jet->p4();
          event_HT += jet->pt();

          float btagSF = 1;
          if (!isData) btagSFHandler.getSFAndEff(theGlobalSyst, jet, btagSF, nullptr);
          if (btagSF != 0.f) SF_btag *= btagSF;
        }
        // block for jet 20
        if (jet->testSelectionBit(AK4JetSelectionHelpers::kTightId) && jet->pt()>=20.f && fabs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight &&
            (!applyPUIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kPUJetId)) &&
            (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
            ) {
          event_Njets20++;
          if (jet->getBtagValue()>=btagvalue_thr) event_Njets20_btagged_loose++;
          if (jet->getBtagValue()>=btag_medium_thr) event_Njets20_btagged_medium++;
          ak4jets20_sump4 += jet->p4();
        }

      }
      event_Njets = ak4jets_tight.size();
      pt_jet1  = (event_Njets > 0)? ak4jets_tight[0]->pt()  : -9;
      eta_jet1 = (event_Njets > 0)? ak4jets_tight[0]->eta() : -9;
      phi_jet1 = (event_Njets > 0)? ak4jets_tight[0]->phi() : -9;
      pt_jet2  = (event_Njets > 1)? ak4jets_tight[1]->pt()  : -9;
      eta_jet2 = (event_Njets > 1)? ak4jets_tight[1]->eta() : -9;
      phi_jet2 = (event_Njets > 1)? ak4jets_tight[1]->phi() : -9;

      if (thePhoton) {
        for (auto* jet : ak4jets) {
          if (isCloseObject(jet->eta(), jet->phi(), eta_photon, phi_photon, 0.3)) {
            cout << __LINE__ << ": event_nJets= " << event_Njets << ", jet->eta()= " << jet->eta() << ", jet->phi()= " << jet->phi() << ", eta_photon= " << eta_photon
                 << ", phi_photon= " << phi_photon << ", selphotons.size()= " << selphotons.size() << ", tightphoton= " << ParticleSelectionHelpers::isLooseParticle(thePhoton) << endl;
          }
        }
      }

      // To deal with the case where the endcap photon may be used, njets 
      if (event_Nphotons >= 1 && thePhoton == nullptr && event_Njets > 2) {
        if (isCloseObject(eta_photon, phi_photon, eta_jet1, phi_jet1, 0.3)) {
          pt_jet1  = ak4jets_tight.at(2)->pt();
          eta_jet1 = ak4jets_tight.at(2)->eta();
          phi_jet1 = ak4jets_tight.at(2)->phi();
          event_Njets--;
        } else if (isCloseObject(eta_photon, phi_photon, eta_jet2, phi_jet2, 0.3)) {
          pt_jet2  = ak4jets_tight.at(2)->pt();
          eta_jet2 = ak4jets_tight.at(2)->eta();
          phi_jet2 = ak4jets_tight.at(2)->phi();
          event_Njets--;
        }
      }

      event_wgt_SFs *= SF_btag;

      mindphi_jet_met = 4.0;
      for (auto const& jet:ak4jets_tight) {
        float dphi_tmp = 4.0;
        HelperFunctions::deltaPhi(float(jet->phi()), phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
        if (dphi_tmp < mindphi_jet_met) mindphi_jet_met = dphi_tmp;
      }
      bool pass_dPhi_jet_met = (mindphi_jet_met > 0.2);
      if (!pass_dPhi_jet_met) continue;
      fill_passedsteps2("_djm");

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptons) {
        if (dilepton->isValid() && dilepton->isOS() && dilepton->nTightDaughters()==2) {
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }
      is_gamma = (nTightDilep < 1) && is_gamma;

      if ((nTightDilep != 1) && event_Nphotons < 1) continue;  // to allow endcap photons to be stored

      is_ee = is_mumu = is_emu = false;
      pass_lepveto = false;
      fill_passedsteps2("_nV");

      ParticleObject::LorentzVector_t boson_p4;
      if (is_gamma) {
        mass_boson = 91.2;  // to be replaced with randomized M
        pt_boson  = thePhoton->pt();
        eta_boson = thePhoton->eta();
        phi_boson = thePhoton->phi();
        boson_p4 = thePhoton->p4();
        boson_p4.SetE(sqrt(pow(thePhoton->p4().P(), 2) + 91.2*91.2));

        pass_lepveto = (electrons.size() == 0) && (muons.size() == 0);
      } else if (nTightDilep == 1) {
        std::vector<ParticleObject*> dileptonDaughters = theChosenDilepton->getDaughters();
        ParticleObjectHelpers::sortByGreaterPt(dileptonDaughters);
        ParticleObject* leadingLepton = dileptonDaughters.front();
        ParticleObject* subleadingLepton = dileptonDaughters.back();

        if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -121) is_ee=true;
        else if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -143) is_emu=true;
        else if (leadingLepton->pdgId() * subleadingLepton->pdgId() == -169) is_mumu=true;

        mass_boson = theChosenDilepton->m();
        pt_boson = theChosenDilepton->pt();
        eta_boson = theChosenDilepton->eta();
        phi_boson = theChosenDilepton->phi();
        boson_p4 = theChosenDilepton->p4();

        id_l1 = leadingLepton->pdgId();
        pt_l1 = leadingLepton->pt();
        eta_l1 = leadingLepton->eta();
        phi_l1 = leadingLepton->phi();
        id_l2 = subleadingLepton->pdgId();
        pt_l2 = subleadingLepton->pt();
        eta_l2 = subleadingLepton->eta();
        phi_l2 = subleadingLepton->phi();

        // if (std::abs(id_l1)==11) electronSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (ElectronObject const*) leadingLepton, isData, false);
        // else muonSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (MuonObject const*) leadingLepton, isData, false);
        // if (std::abs(id_l2)==11) electronSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (ElectronObject const*) subleadingLepton, isData, false);
        // else muonSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (MuonObject const*) subleadingLepton, isData, false);
        if (is_mumu) pass_lepveto = (electrons.size() == 0) && (muons.size() == 2);
        else if (is_ee) pass_lepveto = (electrons.size() == 2) && (muons.size() == 0);
        else if (is_emu) pass_lepveto = (electrons.size() == 1) && (muons.size() == 1);

      } else if (event_Nphotons == 1) {
        mass_boson = 91.2;  // to be replaced with randomized M
        pt_boson  = pt_photon;
        eta_boson = eta_photon;
        phi_boson = phi_photon;
        boson_p4 = ParticleObject::PolarLorentzVector_t(pt_boson, eta_boson, phi_boson, 91.2);
        boson_p4.SetE(sqrt(pow(boson_p4.P(), 2) + 91.2*91.2));

        pass_lepveto = (electrons.size() == 0) && (muons.size() == 0);
      }

      bool is_Zllg = (is_ee || is_mumu) && (event_Nphotons == 1);  // used to lift cuts for llg events

      bool pass_pt_boson = (pt_boson > 50.0);
      if (!is_Zllg && !pass_pt_boson) continue;
      else if (is_Zllg && pt_photon < 50.0) continue;
      fill_passedsteps2("_ptV");

      dphi_boson_met = 4.0;
      HelperFunctions::deltaPhi(phi_boson, phimiss, dphi_boson_met);
      dphi_boson_met = std::abs(dphi_boson_met);
      bool pass_dPhi_boson_met = (dphi_boson_met > 0.5);
      if (!is_Zllg && !pass_dPhi_boson_met) continue;
      fill_passedsteps2("_dbm");
      
      dphi_lljets_met = 4.0;
      ParticleObject::LorentzVector_t p4_lljets = ak4jets_sump4 + boson_p4;
      HelperFunctions::deltaPhi(float(p4_lljets.phi()), phimiss, dphi_lljets_met);
      dphi_lljets_met = std::abs(dphi_lljets_met);
      dphi_lljets20_met = 4.0;
      dphi_jet20_met = 4.0;
      if (event_Njets == 0 && event_Njets20 > 0) {
        HelperFunctions::deltaPhi(float(ak4jets20_sump4.phi()), phimiss, dphi_lljets20_met);
        dphi_lljets20_met = std::abs(dphi_lljets20_met);
        HelperFunctions::deltaPhi(float(ak4jets20_sump4.phi()), phimiss, dphi_jet20_met);
        dphi_jet20_met = std::abs(dphi_jet20_met);
      }

      mTZZ = sqrt(pow(sqrt(pow(pt_boson, 2) + pow(mass_boson, 2)) + sqrt(pow(pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((boson_p4 + pfmet->p4(false, false, false)).Pt(), 2));

      float etamiss_approx = eta_boson;
      ParticleObject::LorentzVector_t pfmet_p4_approx; pfmet_p4_approx = ParticleObject::PolarLorentzVector_t(pTmiss, etamiss_approx, phimiss, PDGHelpers::Zmass);
      ParticleObject::LorentzVector_t pfmet_ZZ_p4_approx = pfmet_p4_approx + boson_p4;
      mZZ = pfmet_ZZ_p4_approx.M();

      event_wgt_trig_muon = eventFilter.getTriggerWeight(triggerCheckList_mu);
      event_wgt_trig_electron = eventFilter.getTriggerWeight(triggerCheckList_el);
      // event_wgt_trig_photon = eventFilter.getTriggerWeight(triggerPropsCheckList_ph, &muons, &electrons, &photons, &ak4jets, &ak8jets, pfmet);;
      // event_wgt_OSSF_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_OSSF, &muons, &electrons, nullptr, nullptr, nullptr, nullptr);
      event_wgt_OSDF_triggers = eventFilter.getTriggerWeight(triggerCheckList_emu);

      // float trigwgt = eventFilter.getTriggerWeight(triggerPropsCheckList, &muons, &electrons, &photons, &ak4jets, &ak8jets, pfmet);
      if (event_wgt_trig_muon==0.f && event_wgt_trig_electron==0.f && event_wgt_OSDF_triggers==0.f && event_wgt_trig_photon==0.f) continue;
      fill_passedsteps2("_trig");

      event_DjjVBF = -1;
      is_VBFcat = false;
      if (ak4jets_tight.size() >= 2){
        is_VBFcat = check_VBF_category(ak4jets_tight, boson_p4);

        std::unordered_map<std::string, float> ME_values;
        SimpleParticleCollection_t daughters;
        daughters.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(pfmet_ZZ_p4_approx)));
        SimpleParticleCollection_t associated;
        associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(ak4jets_tight.at(0)->p4())));
        associated.push_back(SimpleParticle_t(0, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(ak4jets_tight.at(1)->p4())));

        //if (firstValidEvent) CMS3MELAHelpers::melaHandle->setVerbosity(TVar::DEBUG);
        CMS3MELAHelpers::melaHandle->setCandidateDecayMode(TVar::CandidateDecay_Stable);
        CMS3MELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
        MEblock.computeMELABranches();
        MEblock.pushMELABranches();
        MEblock.getBranchValues(ME_values); // Record the MEs into the EDProducer product
        CMS3MELAHelpers::melaHandle->resetInputEvent();

        DjjVBF->update({ ME_values["p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"], ME_values["p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"] }, mZZ);
        event_DjjVBF = float(*DjjVBF);

        if (is_Zllg) {
          ParticleObject::LorentzVector_t rlmet_gamma_p4_approx;
          rlmet_gamma_p4_approx = ParticleObject::PolarLorentzVector_t(pTmiss, eta_photon, phimiss, PDGHelpers::Zmass)
                                + ParticleObject::PolarLorentzVector_t(pt_photon, eta_photon, phi_photon, PDGHelpers::Zmass); 
          daughters.clear();
          daughters.push_back(SimpleParticle_t(25, ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(rlmet_gamma_p4_approx)));

          CMS3MELAHelpers::melaHandle->setInputEvent(&daughters, &associated, nullptr, false);
          MEblock.computeMELABranches();
          MEblock.pushMELABranches();
          MEblock.getBranchValues(ME_values); // Record the MEs into the EDProducer product
          CMS3MELAHelpers::melaHandle->resetInputEvent();

          DjjVBF->update({ ME_values["p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"], ME_values["p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"] }, rlmet_gamma_p4_approx.M());
          event_DjjVBF_rl = float(*DjjVBF);
        }
      }

      // auto fillhists = [&](string s) {
      //   // Z quantities
      //   plot1d("h_mll"+s,   boson.M() , weight, hvec, ";M_{ll} [GeV]" , 125,  0, 500);
      //   plot1d("h_boson_pt"+s,  boson.Pt() , weight, hvec, ";p_{T}(ll) [GeV]" , 200,  0, 800);
      //   plot1d("h_boson_eta"+s, boson.Eta(), weight, hvec, ";#eta(ll) [GeV]"     , 50,  -5, 5);
      //   plot1d("h_mtZZ"+s, mtZZ  , weight, hvec, ";M_{T}(ZZ) [GeV]"   , 200,  0, 2000);
      //   plot1d("h_mZZ"+s, mZZ_approx  , weight, hvec, ";M^{ZZ}(approx) [GeV]"   , 200,  0, 2000);
      // };

      tout->Fill();
      n_evts_acc++;
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc << " / " << (ev_end - ev_start) << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;
  } // End loop over samples list

  foutput->WriteTObject(tout);

  // plot1d("h_sum_wgts", 0, nEntries, hvec, ";Bin;Sum of weights" , 5, 0, 5);
  // plot1d("h_sum_wgts", 1, sum_evts, hvec, ";Bin;Sum of weights" , 5, 0, 5);
  // plot1d("h_sum_wgts", 2, sum_wgts, hvec, ";Bin;Sum of weights" , 5, 0, 5);

  for (auto& h : hvec) {
    if (h.first.find("hnum") == 0 || h.first.find("hden") == 0) continue;
    if (h.first.find("phi") == string::npos && h.first.find("2d") == string::npos)
      moveOverFlowToLastBin1D(h.second);
    string dirname = "OffShell";
    vector<string> jetsufs = {"_eq0j", "_eq1j", "_eq2j", };
    vector<string> lepsufs = {"_gamma", "_ee", "_mumu", "_ll"};
    vector<string> dirsufs;
    for (string jsuf : jetsufs) {
      for (string lsuf : lepsufs) {
        dirsufs.push_back(jsuf+lsuf);
      }
    }
    dirsufs.insert(dirsufs.end(), jetsufs.begin(), jetsufs.end());
    dirsufs.insert(dirsufs.end(), lepsufs.begin(), lepsufs.end());

    for (string dsuf : dirsufs) {
      if (TString(h.first).EndsWith(dsuf.c_str())) {
        dirname += dsuf;
        break;
      }
    }
    TDirectory* dir = (TDirectory*) foutput->Get(dirname.c_str());
    if (dir == nullptr) dir = foutput->mkdir(dirname.c_str());
    dir->cd();
    h.second->Write();
  }

  foutput->Close();

  SampleHelpers::addToCondorTransferList(stroutput);
}

