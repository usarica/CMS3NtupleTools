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


void count(int procsel, int ichunk, int nchunks, TString strdate){
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", "hadoop:200203");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<std::string> triggerCheckList = OffshellTriggerHelpers::getHLTMenus(
    {
      OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kSingleEle, OffshellTriggerHelpers::kMuEle,
      OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kSingleMu
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
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBF_ZZ2L2Nu_M1000_POWHEG", -1, HistogramProperties((int) kYellow-3, 1, 2));

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

  genInfoHandler.setAcquireLHEMEWeights(true);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);

  eventFilter.setTrackDataEvents(true);
  eventFilter.setCheckUniqueDataEvent(true);

  for (size_t isample=0; isample<sampleList.size(); isample++){
    if (procsel>=0 && static_cast<size_t>(procsel)!=isample) continue;

    auto& sample = sampleList.at(isample);

    if (
      sample.name == "ggZZ_2l2nu_BSI"
      ||
      sample.name == "ggZZ_2l2nu_Sig"
      ||
      sample.name == "VBFZZ_2l2nu_BSI"
      ||
      sample.name == "VBFZZ_2l2nu_Sig"
      ) SampleHelpers::configure("2018", "hadoop:200101");
    else{
      SampleHelpers::configure("2018", "hadoop:200203");
      SampleHelpers::setInputDirectory("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims/DileptonSkims/2018");
    }

    bool const isData = SampleHelpers::checkSampleIsData(sample.path);
    if (isData) continue;

    std::vector<TString> sampledirs;
    SampleHelpers::constructSamplesList(sample.path, theGlobalSyst, sampledirs);

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

#define BRANCH_COMMAND(type, name) type name = 0;
    // Event variables
    BRANCH_COMMAND(bool, event_pass_ZZselection);
    BRANCH_COMMAND(bool, event_has_loose_photons);
    BRANCH_COMMAND(bool, event_has_loose_extraLeptons);
    BRANCH_COMMAND(bool, event_has_tight_extraLeptons);
    BRANCH_COMMAND(float, event_wgt);
    BRANCH_COMMAND(float, event_pfmet_pTmiss);
    BRANCH_COMMAND(float, event_pfmet_phimiss);
    BRANCH_COMMAND(float, event_pfmet_mTZZ);
    BRANCH_COMMAND(unsigned int, event_Njets);
    BRANCH_COMMAND(unsigned int, event_Nisotracks);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_fromPV);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_pfcand);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_lost);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_highPurity);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_tight);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_lepId);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_hadId);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_unId);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_lepIdAndCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_hadIdAndCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_unIdAndCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_lepIdAndCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_hadIdAndCuts);
    BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_unIdAndCuts);
    // LL
    BRANCH_COMMAND(float, pt_ll);
    BRANCH_COMMAND(float, eta_ll);
    BRANCH_COMMAND(float, phi_ll);
    BRANCH_COMMAND(float, mass_ll);
    // L1, L2
    BRANCH_COMMAND(int, id_l1);
    BRANCH_COMMAND(float, pt_l1);
    BRANCH_COMMAND(int, id_l2);
    BRANCH_COMMAND(float, pt_l2);
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

    float Nevents_noLepVeto[3][2]={ 0 };
    float Nevents_wLepVeto[3][2]={ 0 };

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
      //MELAout << "Pass line " << __LINE__ << endl;

      float triggerWeight = eventFilter.getTriggerWeight(triggerCheckList);
      event_wgt *= triggerWeight;

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();

      event_wgt *= genInfo->getGenWeight(true);
      if (event_wgt==0.f) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

      if (sample.name == "ggZZ_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
      else if (sample.name == "ggZZ_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
      else if (sample.name == "VBFZZ_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      else if (sample.name == "VBFZZ_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      if (event_wgt==0.f) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

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

      // Tight electrons, muons
      std::vector<MuonObject*> muons_tight;
      std::vector<ElectronObject*> electrons_tight;
      for (auto& muon:muons){ if (ParticleSelectionHelpers::isTightParticle(muon)) muons_tight.push_back(muon); }
      for (auto& electron:electrons){ if (ParticleSelectionHelpers::isTightParticle(electron)) electrons_tight.push_back(electron); }

      dileptonHandler.constructDileptons(&muons, &electrons);
      auto const& dileptons = dileptonHandler.getProducts();
      DileptonObject* theChosenDilepton = nullptr;
      size_t nTightDilep = 0;
      for (auto const& dilepton:dileptons){
        if (dilepton->isValid() && dilepton->isOS() && dilepton->isSF() && dilepton->nTightDaughters()==2){
          if (!theChosenDilepton) theChosenDilepton = dilepton;
          nTightDilep++;
        }
      }
      if (!theChosenDilepton || nTightDilep>1) continue;
      if (theChosenDilepton->pt()<55.f || std::abs(theChosenDilepton->m()-91.2f)>=15.f) continue;

      std::vector<ParticleObject*> dileptonDaughters = theChosenDilepton->getDaughters();
      ParticleObjectHelpers::sortByGreaterPt(dileptonDaughters);
      ParticleObject* leadingLepton = dileptonDaughters.front();
      ParticleObject* subleadingLepton = dileptonDaughters.back();
      if (leadingLepton->pt()<25.f || subleadingLepton->pt()<25.f) continue;

      std::vector<MuonObject*> muons_jetclean;
      std::vector<ElectronObject*> electrons_jetclean;
      for (auto& muon:muons_tight){
        if (theChosenDilepton->hasDaughter(muon)){
          muons_jetclean.push_back(muon);
        }
      }
      for (auto& electron:electrons_tight){
        if (theChosenDilepton->hasDaughter(electron)){
          electrons_jetclean.push_back(electron);
        }
      }

      jetHandler.constructJetMET(theGlobalSyst, &muons_jetclean, &electrons_jetclean, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      event_pfmet_pTmiss = pfmet->pt();
      event_pfmet_phimiss = pfmet->phi();
      if (event_pfmet_pTmiss<125.f) continue;

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_wLepVeto; ak4jets_tight_wLepVeto.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged_wLepVeto; ak4jets_tight_btagged_wLepVeto.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);
          if (jet->extras.pass_leptonVetoId){
            ak4jets_tight_wLepVeto.push_back(jet);
            if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged_wLepVeto.push_back(jet);
          }
        }
      }
      //if (!ak4jets_tight_btagged.empty()) continue;

      if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons_tight, &photons, &ak4jets, &ak8jets)) continue;

      isotrackHandler.constructIsotracks(&muons_tight, &electrons_tight);
      auto const& isotracks = isotrackHandler.getProducts();
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotracks){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;

      mass_ll = theChosenDilepton->m();
      pt_ll = theChosenDilepton->pt();
      eta_ll = theChosenDilepton->eta();
      phi_ll = theChosenDilepton->phi();

      event_pfmet_mTZZ = sqrt(pow(sqrt(pow(pt_ll, 2) + pow(mass_ll, 2)) + sqrt(pow(event_pfmet_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfmet->p4()).Pt(), 2));

      id_l1 = leadingLepton->pdgId();
      pt_l1 = leadingLepton->pt();
      id_l2 = subleadingLepton->pdgId();
      pt_l2 = subleadingLepton->pt();

      float abs_dPhi_ll_pfmet = theChosenDilepton->deltaPhi(event_pfmet_phimiss); abs_dPhi_ll_pfmet = std::abs(abs_dPhi_ll_pfmet);
      if (abs_dPhi_ll_pfmet<=0.5) continue;

      if (ak4jets_tight_btagged.empty()){
        bool doCount = true;
        for (auto const& jet:ak4jets_tight){
          float dphi_tmp;
          HelperFunctions::deltaPhi(float(jet->phi()), event_pfmet_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          if (dphi_tmp<0.5){ doCount=false; break; }
        }
        if (doCount){
          int Njets = ak4jets_tight.size();
          int ix=-1, iy=-1;
          if (Njets==0) ix=0;
          else if (Njets==1) ix=1;
          else if (Njets==2) ix=2;
          iy = (event_pfmet_mTZZ>=350.f ? 1 : 0);
          if (ix>=0) Nevents_noLepVeto[ix][iy] += event_wgt;
        }
      }
      if (ak4jets_tight_btagged_wLepVeto.empty()){
        bool doCount = true;
        for (auto const& jet:ak4jets_tight_wLepVeto){
          float dphi_tmp;
          HelperFunctions::deltaPhi(float(jet->phi()), event_pfmet_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
          if (dphi_tmp<0.5){ doCount=false; break; }
        }
        if (doCount){
          int Njets = ak4jets_tight_wLepVeto.size();
          int ix=-1, iy=-1;
          if (Njets==0) ix=0;
          else if (Njets==1) ix=1;
          else if (Njets==2) ix=2;
          iy = (event_pfmet_mTZZ>350.f ? 1 : 0);
          if (ix>=0) Nevents_wLepVeto[ix][iy] += event_wgt;
        }
      }

    } // End loop over events

    MELAout << "Event counts before lepton veto for sample " << sample_tree.sampleIdentifier << endl;
    for (unsigned int ix=0; ix<3; ix++){
      for (unsigned int iy=0; iy<2; iy++){
        MELAout << "Nj = " << ix << ", ";
        MELAout << (iy==0 ? "mTZZ<350" : "mTZZ>=350") << ": " << Nevents_noLepVeto[ix][iy];
        MELAout << endl;
      }
    }
    MELAout << "Event counts after lepton veto for sample " << sample_tree.sampleIdentifier << endl;
    for (unsigned int ix=0; ix<3; ix++){
      for (unsigned int iy=0; iy<2; iy++){
        MELAout << "Nj = " << ix << ", ";
        MELAout << (iy==0 ? "mTZZ<350" : "mTZZ>=350") << ": " << Nevents_wLepVeto[ix][iy];
        MELAout << endl;
      }
    }

  } // End loop over samples
}
