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


void getTrees(int procsel, int ichunk, int nchunks, TString strdate){
  if (procsel<0) return;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString const coutput_main = "output/TTBarClosure/SkimTrees/" + strdate;
  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", "hadoop:200326");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<std::string> triggerCheckList_OSDF = TriggerHelpers::getHLTMenus(
    {
      TriggerHelpers::kMuEle,
      TriggerHelpers::kSingleEle, TriggerHelpers::kSingleMu
    }
  );
  std::vector<std::string> triggerCheckList_OSSF = TriggerHelpers::getHLTMenus(
    {
      TriggerHelpers::kDoubleEle, TriggerHelpers::kSingleEle,
      TriggerHelpers::kDoubleMu, TriggerHelpers::kSingleMu
    }
  );

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_70-100", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_70-100", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_100-200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_100-200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_200-400", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_200-400", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_400-600", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_400-600", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_600-800", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_600-800", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_800-1200", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_800-1200", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_1200-2500", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_1200-2500", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT_2500-inf", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT_2500-inf", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "GGH_WW2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1500_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));
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
  if (isData) return;

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

  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(true);

  BaseTree sample_tree(SampleHelpers::getDatasetFileName(sampledirs.front()), EVENTS_TREE_NAME, "", "");
  sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());

  const int nEntries = sample_tree.getSelectedNEvents();

  float sum_wgts = (isData ? 1.f : 0.f);
  float xsec = 1;
  if (!isData){
    sample_tree.bookBranch<float>("xsec", 0.f);

    simEventHandler.bookBranches(&sample_tree);
    simEventHandler.wrapTree(&sample_tree);

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    TString firstFile = SampleHelpers::getDatasetFirstFileName(sampledirs.front());
    TFile* fFirstFile = TFile::Open(firstFile, "read");
    TH2F* hCounters = (TH2F*) fFirstFile->Get("cms3ntuple/Counters");
    if (hCounters) sum_wgts = hCounters->GetBinContent(1, 1);
    else{
      MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
      for (int ev=0; ev<nEntries; ev++){
        HelperFunctions::progressbar(ev, nEntries);
        sample_tree.getSelectedEvent(ev);

        genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
        auto const& genInfo = genInfoHandler.getGenInfo();
        float genwgt = genInfo->getGenWeight(true);

        simEventHandler.constructSimEvent(SystematicsHelpers::sNominal);
        float puwgt = simEventHandler.getPileUpWeight();
        sum_wgts += genwgt * puwgt;
      }
    }
    fFirstFile->Close();
  }

  TString stroutput = Form("%s/%s", coutput_main.Data(), sample.name.data());
  if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
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
  BRANCH_COMMAND(unsigned int, event_Njets);
  BRANCH_COMMAND(unsigned int, event_Njets_btagged);
  BRANCH_COMMAND(unsigned int, event_Njets20);
  BRANCH_COMMAND(unsigned int, event_Njets20_btagged);
  BRANCH_COMMAND(float, event_wgt_OSSF_triggers);
  BRANCH_COMMAND(float, event_wgt_OSDF_triggers);
  BRANCH_COMMAND(bool, is_ee);
  BRANCH_COMMAND(bool, is_mumu);
  BRANCH_COMMAND(bool, is_emu);
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
#undef BRANCH_COMMAND

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
  MELAout << "Looping over " << nEntries << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;
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

    eventFilter.constructFilters();
    if (isData && !eventFilter.isUniqueDataEvent()) continue;
    if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

    event_wgt_OSSF_triggers = eventFilter.getTriggerWeight(triggerCheckList_OSSF);
    event_wgt_OSDF_triggers = eventFilter.getTriggerWeight(triggerCheckList_OSDF);
    if ((event_wgt_OSSF_triggers + event_wgt_OSDF_triggers) == 0.f) continue;

    genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
    auto const& genInfo = genInfoHandler.getGenInfo();
    event_wgt *= genInfo->getGenWeight(true);

    if (event_wgt==0.f) continue;

    if (sample.name == "ggZZ_2l2nu_BSI" || sample.name == "ggWW_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
    else if (sample.name == "ggZZ_2l2nu_Sig" || sample.name == "ggWW_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
    else if (sample.name == "ggZZ_2l2nu_BSI_g4" || sample.name == "ggWW_2l2nu_BSI_g4"){ event_wgt *= (genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz4_1_MCFM"]*2.5502 + genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*(-2.55052+std::pow(2.55052, 2)) + genInfo->extras.LHE_ME_weights["p_Gen_GG_BKG_MCFM"]*(-2.5502+1.))*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
    else if (sample.name == "ggZZ_2l2nu_Sig_g4" || sample.name == "ggWW_2l2nu_Sig_g4"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM"]*std::pow(2.55052, 2)*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
    else if (sample.name == "VBFZZ_2l2nu_BSI" || sample.name == "VBFWW_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
    else if (sample.name == "VBFZZ_2l2nu_Sig" || sample.name == "VBFWW_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
    else if (sample.name == "VBFZZ_2l2nu_Sig_g4" || sample.name == "VBFWW_2l2nu_Sig_g4"){ event_wgt *= (genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*std::pow(2.55052, 2) + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv4_1_MCFM"]*(-std::pow(2.55052, 2)+std::pow(2.55052, 4)) + genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BKG_MCFM"]*(-std::pow(2.55052, 2)+1.))*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
    if (event_wgt==0.f) continue;

    simEventHandler.constructSimEvent(SystematicsHelpers::sNominal);
    event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();
    if (event_wgt==0.f) continue;
    //MELAout << "Pass line " << __LINE__ << endl;

    muonHandler.constructMuons(SystematicsHelpers::sNominal);
    electronHandler.constructElectrons(SystematicsHelpers::sNominal);
    photonHandler.constructPhotons(SystematicsHelpers::sNominal);
    particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

    auto const& muons = muonHandler.getProducts();
    auto const& electrons = electronHandler.getProducts();
    auto const& photons = photonHandler.getProducts();

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

    if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons, &photons, &ak4jets, &ak8jets)) continue;

    ParticleObject::LorentzVector_t ak4jets_sump4;
    std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
    std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
    std::vector<AK4JetObject*> ak4jets_tight_pt20; ak4jets_tight_pt20.reserve(ak4jets.size());
    std::vector<AK4JetObject*> ak4jets_tight_pt20_btagged; ak4jets_tight_pt20_btagged.reserve(ak4jets.size());
    for (auto* jet:ak4jets){
      if (ParticleSelectionHelpers::isTightJet(jet)){
        ak4jets_tight.push_back(jet);
        if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);
        ak4jets_sump4 += jet->p4();
      }
      if (
        jet->testSelectionBit(AK4JetSelectionHelpers::kTightId) && jet->testSelectionBit(AK4JetSelectionHelpers::kPUJetId)
        &&
        jet->pt()>=20.f && fabs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
        ){
        ak4jets_tight_pt20.push_back(jet);
        if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_pt20_btagged.push_back(jet);
      }
    }
    event_Njets = ak4jets_tight.size();
    event_Njets_btagged = ak4jets_tight_btagged.size();
    event_Njets20 = ak4jets_tight_pt20.size();
    event_Njets20_btagged = ak4jets_tight_pt20_btagged.size();

    bool pass_dPhi_jet_met = true;
    for (auto const& jet:ak4jets_tight){
      float dphi_tmp;
      HelperFunctions::deltaPhi(float(jet->phi()), phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
      if (dphi_tmp<=0.2){ pass_dPhi_jet_met=false; break; }
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

    if (std::abs(id_l1)==11) electronSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (ElectronObject const*) leadingLepton, isData, false);
    else muonSFHandler.getIdIsoEffAndError(eff_l1, efferr_l1, (MuonObject const*) leadingLepton, isData, false);
    if (std::abs(id_l2)==11) electronSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (ElectronObject const*) subleadingLepton, isData, false);
    else muonSFHandler.getIdIsoEffAndError(eff_l2, efferr_l2, (MuonObject const*) subleadingLepton, isData, false);

    tout->Fill();
  } // End loop over events

  foutput->WriteTObject(tout);
  foutput->Close();

  SampleHelpers::addToCondorTransferList(stroutput);
}

void count(TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure("2018", "hadoop:200203");
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString const cinput_main = "output/TTBarClosure/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Counts";
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = Form("%s/Integrals.txt", coutput_main.Data());
  MELAout.open(stroutput_txt.Data());

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;

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
  /*
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kYellow-3, 1, 2));
  /*
  sampleList.emplace_back("qqWW_2l2nu", "WW#rightarrow2l2#nu", "qqWW_2l2nu", -1, HistogramProperties((int) kTeal-1, 1, 2));
  sampleList.emplace_back("DY_2l_M_50_HT", "DY ll (H_{T}-binned)", "DY_2l_M_50_HT", -1, HistogramProperties((int) kCyan, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_Sig", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig_g4", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig_g4", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_Sig_g4", "gg#rightarrowWW sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_Sig_g4", "EW WW+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_Sig_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_BSI_g4", "gg#rightarrowZZ total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_BSI_g4", "EW ZZ+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggWW_2l2nu_BSI_g4", "gg#rightarrowWW total (#Gamma_{H}=#Gamma_{H}^{SM})", "ggWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("VBFWW_2l2nu_BSI_g4", "EW WW+jj total (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFWW_2l2nu_BSI_g4", -1, HistogramProperties((int) kYellow-3, 1, 2));
  */

  TDirectory* curdir = gDirectory;
  curdir->cd();

  std::vector<SampleSpecs> sampleProcessedList; sampleProcessedList.reserve(sampleList.size());
  std::vector<TFile*> finputlist;
  std::vector<TTree*> treelist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample.path.data() + ".root";
    if (!HostHelpers::FileReadable(cinput)) continue;
    TFile* finput = TFile::Open(cinput, "read");
    TTree* tin = (TTree*) finput->Get("SkimTree");
    if (!tin){
      finput->Close();
      continue;
    }
#define BRANCH_COMMAND(type, name) tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    finputlist.push_back(finput);
    treelist.push_back(tin);
    MELAout << "Extracted input file " << cinput << endl;
    curdir->cd();
    sampleProcessedList.push_back(sample);
  }
  const size_t nsamples = sampleProcessedList.size();

  std::vector<TString> channel_labels={ "ee", "mumu", "emu" };
  std::vector<TString> in_out_labels={ "in", "out" };

  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << isample << " (" << sampleProcessedList.at(isample).name << ")" << endl;
    std::vector<HistogramObject>& hspecs = sampleProcessedList.at(isample).hlist_1D;
    TFile*& finput = finputlist.at(isample);
    finput->cd();
    TTree* const& tin = treelist.at(isample);

    float sum_wgts[3][3][2]={ { { 0 } } }; // Channel (ee, mumu, emu), Njets, in/out-like selection
    float sum_wgts_emu_ee_like[3][2]={ { 0 } }; // Njets, in/out-like selection
    float sum_wgts_emu_mumu_like[3][2]={ { 0 } }; // Njets, in/out-like selection

    int nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);

      if (pt_ll<55.f || pt_l1<25.f || pt_l2<25.f) continue;
      int index_channel = 0*is_ee + 1*is_mumu + 2*is_emu;

      if (index_channel<2) event_wgt *= event_wgt_OSSF_triggers;
      else event_wgt *= event_wgt_OSDF_triggers;

      if (mass_ll<40.f || mass_ll>=200.f) continue;

      float eff_weight_ee_like=1;
      float eff_weight_mumu_like=1;
      if (index_channel==2){
        eff_weight_ee_like = eff_weight_mumu_like = 1.f/(eff_l1*eff_l2);

        float dummy_err=0;
        float eff_e1=1, eff_e2=1, eff_m1=1, eff_m2=1;

        electronSFHandler.getIdIsoEffAndError(eff_e1, dummy_err, pt_l1, eta_l1, false, false);
        electronSFHandler.getIdIsoEffAndError(eff_e2, dummy_err, pt_l2, eta_l2, false, false);

        muonSFHandler.getIdIsoEffAndError(eff_m1, dummy_err, pt_l1, eta_l1, false, false);
        muonSFHandler.getIdIsoEffAndError(eff_m2, dummy_err, pt_l2, eta_l2, false, false);

        eff_weight_ee_like *= eff_e1*eff_e2;
        eff_weight_mumu_like *= eff_m1*eff_m2;
      }

      int index_Njets = static_cast<int>(event_Njets)-1;
      int index_INorOUT=-1;
      if (std::abs(mass_ll-91.2f)<15.f){
        index_INorOUT = 0;
        if (pTmiss<125.f) continue;
        if (event_Njets==0 || event_Njets_btagged!=0) continue;
      }
      else{
        index_INorOUT = 1;
        if (pTmiss<70.f) continue;
        if (event_Njets_btagged==0) continue;
      }

      sum_wgts[index_channel][index_Njets][index_INorOUT] += event_wgt;
      if (index_channel==2){
        sum_wgts_emu_ee_like[index_Njets][index_INorOUT] += event_wgt*eff_weight_ee_like;
        sum_wgts_emu_mumu_like[index_Njets][index_INorOUT] += event_wgt*eff_weight_mumu_like;
      }
    }

    for (size_t ic=0; ic<channel_labels.size(); ic++){
      for (int ijet=1; ijet<=3; ijet++){
        for (size_t io=0; io<in_out_labels.size(); io++){
          MELAout << channel_labels.at(ic) << ", ";
          MELAout << "Njets = " << ijet << ", ";
          MELAout << in_out_labels.at(io) << ": " << sum_wgts[ic][ijet-1][io] << endl;
        }
      }
    }
    MELAout << "Weighted emu contributions:" << endl;
    for (int ijet=1; ijet<=3; ijet++){
      for (size_t io=0; io<in_out_labels.size(); io++){
        MELAout << "ee-like, ";
        MELAout << "Njets = " << ijet << ", ";
        MELAout << in_out_labels.at(io) << ": " << sum_wgts_emu_ee_like[ijet-1][io] << endl;

        MELAout << "mumu-like, ";
        MELAout << "Njets = " << ijet << ", ";
        MELAout << in_out_labels.at(io) << ": " << sum_wgts_emu_mumu_like[ijet-1][io] << endl;
      }
    }

    curdir->cd();
  }

  for (TFile*& finput:finputlist) finput->Close();
  MELAout.close();
}
