#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include <DataFormats/MuonReco/interface/Muon.h>


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
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString const coutput_main = "output/IsotracksFromData/SkimTrees/" + strdate;
  gSystem->mkdir(coutput_main, true);

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
  genInfoHandler.setAcquireGenParticles(true);

  eventFilter.setTrackDataEvents(false);
  eventFilter.setCheckUniqueDataEvent(false);

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

    TString stroutput = Form("%s/%s", coutput_main.Data(), sample.name.data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    TTree* tout = new TTree("SkimTree", "");
#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
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
    // L3
    BRANCH_COMMAND(bool, has_l3);
    BRANCH_COMMAND(int, id_l3);
    BRANCH_COMMAND(bool, isPrompt_genmatch_l3);
    BRANCH_COMMAND(int, genmatch_id_l3);
    BRANCH_COMMAND(float, dRmatch_genpart_l3);
    BRANCH_COMMAND(float, pt_l3);
    BRANCH_COMMAND(float, eta_l3);
    BRANCH_COMMAND(float, phi_l3);
    BRANCH_COMMAND(float, mass_l3);
    BRANCH_COMMAND(bool, passLooseSelection_l3);
    BRANCH_COMMAND(bool, passTightSelection_l3);
    // Isotrack closest to L3
    BRANCH_COMMAND(bool, has_isotrack_closest_l3);
    BRANCH_COMMAND(bool, fromPV_isotrack_closest_l3);
    BRANCH_COMMAND(bool, is_pfCand_isotrack_closest_l3);
    BRANCH_COMMAND(bool, is_lostTrack_isotrack_closest_l3);
    BRANCH_COMMAND(bool, is_highPurityTrack_isotrack_closest_l3);
    BRANCH_COMMAND(bool, is_tightTrack_isotrack_closest_l3);
    BRANCH_COMMAND(int, id_isotrack_closest_l3);
    BRANCH_COMMAND(float, dRmatch_isotrack_closest_l3);
    BRANCH_COMMAND(float, pt_isotrack_closest_l3);
    BRANCH_COMMAND(float, eta_isotrack_closest_l3);
    BRANCH_COMMAND(float, phi_isotrack_closest_l3);
    BRANCH_COMMAND(float, mass_isotrack_closest_l3);
    BRANCH_COMMAND(float, pfIso03_ch_isotrack_closest_l3);
    BRANCH_COMMAND(float, dxy_isotrack_closest_l3);
    BRANCH_COMMAND(float, dz_isotrack_closest_l3);
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

      event_has_loose_photons = false;
      event_has_loose_extraLeptons = false;
      event_has_tight_extraLeptons = false;

      id_l1 = 0;
      pt_l1 = -1;
      id_l2 = 0;
      pt_l2 = -1;

      has_l3 = false;
      id_l3 = 0;
      isPrompt_genmatch_l3 = false;
      genmatch_id_l3 = 0;
      dRmatch_genpart_l3 = -1;
      pt_l3 = -1;
      eta_l3 = 0;
      phi_l3 = 0;
      mass_l3 = 0;
      passLooseSelection_l3 = false;
      passTightSelection_l3 = false;

      has_isotrack_closest_l3 = false;
      fromPV_isotrack_closest_l3 = false;
      is_pfCand_isotrack_closest_l3 = false;
      is_lostTrack_isotrack_closest_l3 = false;
      is_highPurityTrack_isotrack_closest_l3 = false;
      is_tightTrack_isotrack_closest_l3 = false;
      id_isotrack_closest_l3 = 0;
      dRmatch_isotrack_closest_l3 = -1;
      pt_isotrack_closest_l3 = -1;
      eta_isotrack_closest_l3 = 0;
      phi_isotrack_closest_l3 = 0;
      mass_isotrack_closest_l3 = 0;
      pfIso03_ch_isotrack_closest_l3 = 0;
      dxy_isotrack_closest_l3 = 0;
      dz_isotrack_closest_l3 = 0;

      eventFilter.constructFilters();
      //if (isData && !eventFilter.isUniqueDataEvent()) continue;
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

      float triggerWeight = eventFilter.getTriggerWeight(triggerCheckList);
      event_wgt *= triggerWeight;

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();
      auto const& genparticles = genInfoHandler.getGenParticles();

      event_wgt *= genInfo->getGenWeight(true);
      if (event_wgt==0.f) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

      if (sample.name == "ggZZ_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
      else if (sample.name == "ggZZ_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]*genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"]; }
      else if (sample.name == "VBFZZ_2l2nu_BSI"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      else if (sample.name == "VBFZZ_2l2nu_Sig"){ event_wgt *= genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"]*genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"]; }
      if (event_wgt==0.f){
        /*
        if (sample.name == "ggZZ_2l2nu_BSI"){
          MELAout << "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM = " << genInfo->extras.LHE_ME_weights["p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM"] << endl;
          MELAout << "p_Gen_CPStoBWPropRewgt = " << genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"] << endl;
          MELAout << "KFactor_QCD_NNLO_ggZZ_Sig_Nominal = " << genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"] << endl;
        }
        else if (sample.name == "ggZZ_2l2nu_Sig"){
          MELAout << "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM = " << genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"] << endl;
          MELAout << "p_Gen_CPStoBWPropRewgt = " << genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"] << endl;
          MELAout << "KFactor_QCD_NNLO_ggZZ_Sig_Nominal = " << genInfo->extras.Kfactors["KFactor_QCD_NNLO_ggZZ_Sig_Nominal"] << endl;
        }
        else if (sample.name == "VBFZZ_2l2nu_BSI"){
          MELAout << "p_Gen_JJEW_BSI_ghv1_1_MCFM = " << genInfo->extras.LHE_ME_weights["p_Gen_JJEW_BSI_ghv1_1_MCFM"] << endl;
          MELAout << "p_Gen_CPStoBWPropRewgt = " << genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"] << endl;
        }
        else if (sample.name == "VBFZZ_2l2nu_Sig"){
          MELAout << "p_Gen_JJEW_SIG_ghv1_1_MCFM = " << genInfo->extras.LHE_ME_weights["p_Gen_JJEW_SIG_ghv1_1_MCFM"] << endl;
          MELAout << "p_Gen_CPStoBWPropRewgt = " << genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"] << endl;
        }
        */
        continue;
      }
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
      for (auto const& photon:photons){
        if (ParticleSelectionHelpers::isLooseParticle(photon)){
          event_has_loose_photons = true;
          break;
        }
      }

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
          float distance_Zpeak = std::abs(dilepton->m()-91.2);
          if (!theChosenDilepton || distance_Zpeak<std::abs(theChosenDilepton->m()-91.2)) theChosenDilepton = dilepton;
        }
      }
      if (!theChosenDilepton) continue;
      //MELAout << "Pass Z presence" << endl;

      std::vector<ParticleObject*> dileptonDaughters = theChosenDilepton->getDaughters();
      ParticleObjectHelpers::sortByGreaterPt(dileptonDaughters);
      ParticleObject* leadingLepton = dileptonDaughters.front();
      ParticleObject* subleadingLepton = dileptonDaughters.back();
      if (ev%10000 == 0){
        MELAout << "leadingLepton p4 = " << leadingLepton->p4() << endl;
        MELAout << "subleadingLepton p4 = " << subleadingLepton->p4() << endl;
      }

      std::vector<MuonObject*> muons_jetclean;
      std::vector<ElectronObject*> electrons_jetclean;
      for (auto& muon:muons_tight){
        if (theChosenDilepton->hasDaughter(muon)){
          muons_jetclean.push_back(muon);
          if (ev%10000 == 0) MELAout << "Muon with p4 = " << muon->p4() << " can be used for cleaning." << endl;
        }
      }
      for (auto& electron:electrons_tight){
        if (theChosenDilepton->hasDaughter(electron)){
          electrons_jetclean.push_back(electron);
          if (ev%10000 == 0) MELAout << "Electron with p4 = " << electron->p4() << " can be used for cleaning." << endl;
        }
      }

      jetHandler.constructJetMET(theGlobalSyst, &muons_jetclean, &electrons_jetclean, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      event_pfmet_pTmiss = pfmet->pt();
      event_pfmet_phimiss = pfmet->phi();

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);
        }
      }
      if (!ak4jets_tight_btagged.empty()) continue;
      event_Njets = ak4jets_tight.size();
      //MELAout << "Pass line " << __LINE__ << endl;

      if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons_tight, &photons, &ak4jets, &ak8jets)) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

      isotrackHandler.constructIsotracks(nullptr, nullptr);
      auto const& isotracks = isotrackHandler.getProducts();
      event_Nisotracks = isotracks.size();
      if (event_Nisotracks==0) continue;
      //MELAout << "Pass line " << __LINE__ << endl;

      mass_ll = theChosenDilepton->m();
      pt_ll = theChosenDilepton->pt();
      eta_ll = theChosenDilepton->eta();
      phi_ll = theChosenDilepton->phi();

      event_pfmet_mTZZ = sqrt(pow(sqrt(pow(pt_ll, 2) + pow(mass_ll, 2)) + sqrt(pow(event_pfmet_pTmiss, 2) + pow(PDGHelpers::Zmass, 2)), 2) - pow((theChosenDilepton->p4() + pfmet->p4()).Pt(), 2));

      id_l1 = leadingLepton->pdgId();
      pt_l1 = leadingLepton->pt();
      id_l2 = subleadingLepton->pdgId();
      pt_l2 = subleadingLepton->pt();

      event_pass_ZZselection = (
        std::abs(mass_ll-91.2f)<15.f
        &&
        pt_ll>=55.f && pt_l1>=25.f && pt_l2>=25.f
        &&
        event_pfmet_pTmiss>=125.f
        );
      for (auto const& jet:ak4jets_tight){
        float dphi_tmp;
        HelperFunctions::deltaPhi(float(jet->phi()), event_pfmet_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
        if (dphi_tmp<0.5){ event_pass_ZZselection=false; break; }
      }
      float abs_dPhi_ll_pfmet = theChosenDilepton->deltaPhi(event_pfmet_phimiss); abs_dPhi_ll_pfmet = std::abs(abs_dPhi_ll_pfmet);
      event_pass_ZZselection &= abs_dPhi_ll_pfmet>0.5;

      std::vector<IsotrackObject*> isotracks_cleaned;
      for (auto& isotrack:isotracks){
        if (reco::deltaR(isotrack->p4(), leadingLepton->p4())<0.3) continue;
        else if (reco::deltaR(isotrack->p4(), subleadingLepton->p4())<0.3) continue;
        isotracks_cleaned.push_back(isotrack);
      }

      std::vector<ParticleObject*> extraLeptons;
      for (auto& part:muons){
        if (theChosenDilepton->hasDaughter(part)) continue;
        extraLeptons.push_back(part);
        if (ParticleSelectionHelpers::isTightParticle(part)) event_has_tight_extraLeptons = true;
        // Manual selection for loose and soft muons based on older selection requirements
        auto const& extras = part->extras;
        bool is_loose = (
          MuonSelectionHelpers::relPFIso_DR0p4(*part)<0.15
          &&
          (extras.POG_selector_bits & Muon::CutBasedIdLoose) == Muon::CutBasedIdLoose
          &&
          part->testSelectionBit(MuonSelectionHelpers::kLooseKin)
          );
        if (is_loose) event_has_loose_extraLeptons = true;
      }
      for (auto& part:electrons){
        if (theChosenDilepton->hasDaughter(part)) continue;
        extraLeptons.push_back(part);
        if (ParticleSelectionHelpers::isTightParticle(part)) event_has_tight_extraLeptons = true;
        // Manual selection for loose electrons based on older selection requirements
        auto const& extras = part->extras;
        bool is_loose = (
          ElectronSelectionHelpers::relPFIso_DR0p3(*part)<0.2
          &&
          extras.id_cutBased_Fall17V2_Loose_Bits == 1023
          &&
          part->testSelectionBit(ElectronSelectionHelpers::kLooseKin)
          );
        if (is_loose) event_has_loose_extraLeptons = true;
      }
      ParticleObjectHelpers::sortByGreaterPt(extraLeptons);

      // Find isotrack closest to L3
      auto const* lepton_l3 = (extraLeptons.empty() ? nullptr : extraLeptons.front());
      IsotrackObject* isotrack_closest_l3 = nullptr;
      if (lepton_l3){
        has_l3 = true;
        id_l3 = lepton_l3->pdgId();
        pt_l3 = lepton_l3->pt();
        eta_l3 = lepton_l3->eta();
        phi_l3 = lepton_l3->phi();
        mass_l3 = lepton_l3->m();

        GenParticleObject* genparticle_matched_l3 = nullptr;
        for (auto const& genparticle:genparticles){
          if (genparticle->status()!=1) continue;
          if (PDGHelpers::isANeutrino(genparticle->pdgId())) continue;
          float tmpDR = reco::deltaR(genparticle->p4(), lepton_l3->p4());
          if (dRmatch_genpart_l3<0.f || tmpDR<dRmatch_genpart_l3){
            dRmatch_genpart_l3 = tmpDR;
            genparticle_matched_l3 = genparticle;
          }
        }
        if (genparticle_matched_l3){
          isPrompt_genmatch_l3 = genparticle_matched_l3->extras.isPromptFinalState;
          genmatch_id_l3 = genparticle_matched_l3->pdgId();
        }

        ElectronObject const* electron_l3 = dynamic_cast<ElectronObject const*>(lepton_l3);
        MuonObject const* muon_l3 = dynamic_cast<MuonObject const*>(lepton_l3);

        passTightSelection_l3 = ParticleSelectionHelpers::isTightParticle(lepton_l3);
        passLooseSelection_l3 = (
          (
            electron_l3
            &&
            ElectronSelectionHelpers::relPFIso_DR0p3(*electron_l3)<0.2
            &&
            electron_l3->extras.id_cutBased_Fall17V2_Loose_Bits == 1023
            &&
            electron_l3->testSelectionBit(ElectronSelectionHelpers::kLooseKin)
            )
          ||
          (
            muon_l3
            &&
            MuonSelectionHelpers::relPFIso_DR0p4(*muon_l3)<0.15
            &&
            (muon_l3->extras.POG_selector_bits & Muon::CutBasedIdLoose) == Muon::CutBasedIdLoose
            &&
            muon_l3->testSelectionBit(MuonSelectionHelpers::kLooseKin)
            )
          );

        for (auto& isotrack:isotracks_cleaned){
          float tmpDR = reco::deltaR(isotrack->p4(), lepton_l3->p4());
          if (!isotrack_closest_l3 || dRmatch_isotrack_closest_l3>tmpDR){
            isotrack_closest_l3 = isotrack;
            dRmatch_isotrack_closest_l3 = tmpDR;
          }
        }
      }

      event_Nisotracks_fromPV = 0;
      event_Nisotracks_pfcand = 0;
      event_Nisotracks_lost = 0;
      event_Nisotracks_highPurity = 0;
      event_Nisotracks_tight = 0;
      event_Nisotracks_lepId = 0;
      event_Nisotracks_hadId = 0;
      event_Nisotracks_unId = 0;
      event_Nisotracks_veto_lepIdAndCuts = 0;
      event_Nisotracks_veto_hadIdAndCuts = 0;
      event_Nisotracks_veto_unIdAndCuts = 0;
      event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts = 0;
      event_Nisotracks_stdveto_lepIdAndCuts = 0;
      event_Nisotracks_stdveto_hadIdAndCuts = 0;
      event_Nisotracks_stdveto_unIdAndCuts = 0;
      for (auto& isotrack:isotracks_cleaned){
        float pt = isotrack->pt();
        auto const& extras = isotrack->extras;

        if (extras.fromPV) event_Nisotracks_fromPV++;
        if (extras.is_pfCand) event_Nisotracks_pfcand++;
        if (extras.is_lostTrack) event_Nisotracks_lost++;
        if (extras.is_highPurityTrack) event_Nisotracks_highPurity++;
        if (extras.is_tightTrack) event_Nisotracks_tight++;
        if (std::abs(isotrack->pdgId())==11 || std::abs(isotrack->pdgId())==13){
          event_Nisotracks_lepId++;
          if (extras.fromPV && (extras.is_highPurityTrack || std::abs(isotrack->pdgId())==11) && std::abs(extras.dxy)<0.1 && std::abs(extras.dz)<0.1 && extras.pfIso03_ch<5.){
            event_Nisotracks_veto_lepIdAndCuts++;
          }
          if (extras.fromPV && extras.is_pfCand/* && std::abs(extras.dxy)<0.1*/ && std::abs(extras.dz)<0.1 && extras.pfIso03_ch<std::min(5., 0.1*isotrack->pt())){
            event_Nisotracks_stdveto_lepIdAndCuts++;
          }
        }
        if (std::abs(isotrack->pdgId())>100){
          event_Nisotracks_hadId++;
          if (extras.fromPV && std::abs(extras.dxy)<0.05 && std::abs(extras.dz)<0.1 && extras.pfIso03_ch<5.){
            event_Nisotracks_veto_hadIdAndCuts++;
            for (auto const& jet:ak4jets_tight){
              float dR_isotrack_jet = reco::deltaR(isotrack->p4(), jet->p4());
              if (
                dR_isotrack_jet<0.1
                &&
                extras.pfIso03_ch<1.f
                ){
                event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts++;
                break;
              }
            }
          }
          if (extras.fromPV && extras.is_pfCand /*&& std::abs(extras.dxy)<0.1*/ && std::abs(extras.dz)<0.1 && isotrack->pt()>=10.f && extras.pfIso03_ch<std::min(5., 0.1*isotrack->pt())){
            event_Nisotracks_stdveto_hadIdAndCuts++;
          }
        }
        if (std::abs(isotrack->pdgId())==0){
          event_Nisotracks_unId++;
          if (extras.fromPV && extras.is_tightTrack && std::abs(extras.dxy)<0.05 && std::abs(extras.dz)<0.1 && extras.pfIso03_ch<5.f){
            event_Nisotracks_veto_unIdAndCuts++;
          }
          if (extras.fromPV && extras.is_pfCand/* && std::abs(extras.dxy)<0.1*/ && std::abs(extras.dz)<0.1 && isotrack->pt()>=10.f && extras.pfIso03_ch<std::min(5., 0.1*isotrack->pt())){
            event_Nisotracks_stdveto_unIdAndCuts++;
          }
        }
      }

      if (isotrack_closest_l3){
        auto const& extras_isotrack_closest_l3 = isotrack_closest_l3->extras;

        has_isotrack_closest_l3 = true;
        fromPV_isotrack_closest_l3 = extras_isotrack_closest_l3.fromPV;
        is_pfCand_isotrack_closest_l3 = extras_isotrack_closest_l3.is_pfCand;
        is_lostTrack_isotrack_closest_l3 = extras_isotrack_closest_l3.is_lostTrack;
        is_highPurityTrack_isotrack_closest_l3 = extras_isotrack_closest_l3.is_highPurityTrack;
        is_tightTrack_isotrack_closest_l3 = extras_isotrack_closest_l3.is_tightTrack;
        id_isotrack_closest_l3 = isotrack_closest_l3->pdgId();
        //dRmatch_isotrack_closest_l3 = reco::deltaR(isotrack_closest_l3->p4(), lepton_l3->p4());
        pt_isotrack_closest_l3 = isotrack_closest_l3->pt();
        eta_isotrack_closest_l3 = isotrack_closest_l3->eta();
        phi_isotrack_closest_l3 = isotrack_closest_l3->phi();
        mass_isotrack_closest_l3 = isotrack_closest_l3->m();
        pfIso03_ch_isotrack_closest_l3 = extras_isotrack_closest_l3.pfIso03_ch;
        dxy_isotrack_closest_l3 = extras_isotrack_closest_l3.dxy;
        dz_isotrack_closest_l3 = extras_isotrack_closest_l3.dz;
      }

      tout->Fill();
    } // End loop over events

    foutput->WriteTObject(tout);
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples
}

void makePlots(TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();
  const float lumi = SampleHelpers::getIntegratedLuminosity("2018");

  TString const cinput_main = "output/IsotracksFromData/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Plots" + "/Nevents";
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = Form("%s/Integrals.txt", coutput_main.Data());
  MELAout.open(stroutput_txt.Data());

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
  BRANCH_COMMAND(bool, event_pass_ZZselection) \
  BRANCH_COMMAND(bool, event_has_loose_photons) \
  BRANCH_COMMAND(bool, event_has_loose_extraLeptons) \
  BRANCH_COMMAND(bool, event_has_tight_extraLeptons) \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_pfmet_pTmiss) \
  BRANCH_COMMAND(float, event_pfmet_phimiss) \
  BRANCH_COMMAND(float, event_pfmet_mTZZ) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_fromPV) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_pfcand) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_lost) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_highPurity) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_tight) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_lepId) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_hadId) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_unId) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_lepIdAndCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_hadIdAndCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_unIdAndCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_lepIdAndCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_hadIdAndCuts) \
  BRANCH_COMMAND(unsigned int, event_Nisotracks_stdveto_unIdAndCuts) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(int, id_l1) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(int, id_l2) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(bool, has_l3) \
  BRANCH_COMMAND(int, id_l3) \
  BRANCH_COMMAND(bool, isPrompt_genmatch_l3) \
  BRANCH_COMMAND(int, genmatch_id_l3) \
  BRANCH_COMMAND(float, dRmatch_genpart_l3) \
  BRANCH_COMMAND(float, pt_l3) \
  BRANCH_COMMAND(float, eta_l3) \
  BRANCH_COMMAND(float, phi_l3) \
  BRANCH_COMMAND(float, mass_l3) \
  BRANCH_COMMAND(bool, passLooseSelection_l3) \
  BRANCH_COMMAND(bool, passTightSelection_l3) \
  BRANCH_COMMAND(bool, has_isotrack_closest_l3) \
  BRANCH_COMMAND(bool, fromPV_isotrack_closest_l3) \
  BRANCH_COMMAND(bool, is_pfCand_isotrack_closest_l3) \
  BRANCH_COMMAND(bool, is_lostTrack_isotrack_closest_l3) \
  BRANCH_COMMAND(bool, is_highPurityTrack_isotrack_closest_l3) \
  BRANCH_COMMAND(bool, is_tightTrack_isotrack_closest_l3) \
  BRANCH_COMMAND(int, id_isotrack_closest_l3) \
  BRANCH_COMMAND(float, dRmatch_isotrack_closest_l3) \
  BRANCH_COMMAND(float, pt_isotrack_closest_l3) \
  BRANCH_COMMAND(float, eta_isotrack_closest_l3) \
  BRANCH_COMMAND(float, phi_isotrack_closest_l3) \
  BRANCH_COMMAND(float, mass_isotrack_closest_l3) \
  BRANCH_COMMAND(float, pfIso03_ch_isotrack_closest_l3) \
  BRANCH_COMMAND(float, dxy_isotrack_closest_l3) \
  BRANCH_COMMAND(float, dz_isotrack_closest_l3)
#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("TT_2l2nu", "t#bar{t}#rightarrow2l2#nu", "TT_2l2nu", -1, HistogramProperties((int) kOrange-3, 1, 2));
  sampleList.emplace_back("qqZZ_2l2nu", "ZZ#rightarrow2l2#nu", "qqZZ_2l2nu", -1, HistogramProperties((int) kRed, 1, 2));
  sampleList.emplace_back("qqWZ_3lnu_MG", "WZ#rightarrow3l1#nu", "qqWZ_3lnu_MG", -1, HistogramProperties((int) kYellow-3, 1, 2));
  //sampleList.emplace_back("qqWZ_lnu2q", "WZ#rightarrowl#nuq", "qqWZ_lnu2q", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("ggZZ_2l2nu_Sig", "gg#rightarrowZZ sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "ggZZ_2l2nu_Sig", -1, HistogramProperties((int) kBlue, 2, 2));
  sampleList.emplace_back("VBFZZ_2l2nu_Sig", "EW ZZ+jj sig. (#Gamma_{H}=#Gamma_{H}^{SM})", "VBFZZ_2l2nu_Sig", -1, HistogramProperties((int) kViolet, 2, 2));

  TDirectory* curdir = gDirectory;
  curdir->cd();

  const size_t nsamples = sampleList.size();
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
  }
  curdir->cd();

  TString strChannel, strChannelLabel;
  getChannelTitleLabel(3, strChannel, strChannelLabel);

  std::vector< std::vector<CutSpecs> > cutsets; getCutSets(cutsets);

  std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;
  std::vector<TString> hnames;
  for (size_t isample=0; isample<nsamples; isample++){
    MELAout << "Processing sample " << sampleList.at(isample).name << endl;
    std::vector<HistogramObject>& hspecs = sampleList.at(isample).hlist_1D;
    TFile*& finput = finputlist.at(isample);
    finput->cd();
    TTree* const& tin = treelist.at(isample);

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
      MELAout << "\t- Creating histograms for cut set " << cuttitle << endl;

      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pfmet_pTmiss", cuttitle.Data()),
        Form("%s|%s|%s", sampleList.at(isample).label.data(), strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{miss,PF} (GeV)", "Events / 25 GeV", 35, 125, 1000, finput
      );
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pfmet_mTZZ", cuttitle.Data()),
        Form("%s|%s|%s", sampleList.at(isample).label.data(), strChannelLabel.Data(), cutlabel.Data()),
        "m_{T}^{ZZ,PF} (GeV)", "Events / 25 GeV", 80, 200, 2200, finput
      );
    }

    sampleList.at(isample).setup();
    std::vector<TH1F*> hlist; hlist.reserve(hspecs.size());
    for (auto& hspec:hspecs){
      TH1F* htmp = &(hspec.hist); hlist.push_back(htmp);
      htmp->SetLineWidth(2);
      htmp->SetFillStyle(0);
      if (isample==0) hnames.push_back(htmp->GetName());
    }

    for (int ev=0; ev<tin->GetEntries(); ev++){
      tin->GetEntry(ev);

      if (!event_pass_ZZselection) continue;

      auto it_hist = hlist.begin();
      size_t icutset=0;
      for (auto const& cutset:cutsets){
        bool doFill=true;
        for (auto const& cut:cutset){
          TString const& cutvar = cut.cutvar;
          float cutval=0;
          if (cutvar == "event_has_loose_extraLeptons_flag") cutval = event_has_loose_extraLeptons;
          else if (cutvar == "event_has_isotracks_stdveto_lepId_flag") cutval = event_Nisotracks_stdveto_lepIdAndCuts>0;
          else if (cutvar == "event_has_isotracks_stdveto_anyId_flag") cutval = event_Nisotracks_stdveto_lepIdAndCuts>0 || event_Nisotracks_stdveto_hadIdAndCuts>0;
          else if (cutvar == "event_has_isotracks_stdveto_anyId_or_tightLeptons_flag") cutval = event_Nisotracks_stdveto_lepIdAndCuts>0 || event_Nisotracks_stdveto_hadIdAndCuts>0 || event_has_tight_extraLeptons;
          else if (cutvar == "event_has_isotracks_veto_lepId_flag") cutval = event_Nisotracks_veto_lepIdAndCuts>0;
          else if (cutvar == "event_has_isotracks_veto_lepId_hadId_flag") cutval = event_Nisotracks_veto_lepIdAndCuts>0 || event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts>0;
          else if (cutvar == "event_has_isotracks_veto_lepId_hadId_or_tightLeptons_flag") cutval = event_Nisotracks_veto_lepIdAndCuts>0 || event_Nisotracks_veto_hadIdAndCuts_jetMatchingCuts>0 || event_has_tight_extraLeptons;
          else if (cutvar == "Nj") cutval = event_Njets;
          doFill &= cut.testCut(cutval);

          //if (ev<10) MELAout << "\t\t-Tested cut " << cutvar << " | " << cut.cutlow << " | " << cut.cuthigh << ": " << doFill << endl;
        }

        {
          if (doFill) (*it_hist)->Fill(event_pfmet_pTmiss, event_wgt); it_hist++;
          if (doFill) (*it_hist)->Fill(event_pfmet_mTZZ, event_wgt); it_hist++;
        }

        icutset++;
      }
      if (it_hist!=hlist.end()){
        MELAerr << "ERROR: Not all histograms are filled!" << endl;
      }
    }

    // Post-processing
    for (TH1F* htmp:hlist){
      HelperFunctions::wipeOverUnderFlows(htmp, false, true);
      double inthist = htmp->Integral(1, htmp->GetNbinsX());
      //MELAout << "Histogram " << htmp->GetName() << " has integral " << inthist << endl;
      if (TString(htmp->GetName()).Contains("mTZZ")){
        int x350 = htmp->GetXaxis()->FindBin(350.);
        float inthist_x350 = htmp->Integral(x350, htmp->GetNbinsX());
        MELAout << "Histogram " << htmp->GetName() << " has integral (mTZZ>350) " << inthist_x350 << endl;
      }
    }

    sample_hist_map[sampleList.at(isample).name]=hlist;
  }

  for (size_t iplot=0; iplot<hnames.size(); iplot++){
    curdir->cd();

    TString const& histname = hnames.at(iplot);
    int nbins = sample_hist_map[sampleList.front().name].at(iplot)->GetNbinsX();

    bool useLogY = true;

    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      if (sample.name.find("Sig")!=std::string::npos) continue;
      TH1F* hist = sample_hist_map[sample.name].at(iplot);

      for (size_t js=is+1; js<sampleList.size(); js++){
        auto& sample_j = sampleList.at(js);
        if (sample_j.name.find("Sig")!=std::string::npos) continue;
        TH1F* hist_j = sample_hist_map[sample_j.name].at(iplot);
        hist->Add(hist_j, 1.);
      }
      hist->SetFillStyle(3018);
    }

    size_t nplottables=0;
    double ymin = 0;
    double ymax = -1;
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample.name].at(iplot);

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
      TH1F* hist = sample_hist_map[sample.name].at(iplot);
      hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    TString canvasname = hnames.at(iplot);
    HelperFunctions::replaceString(canvasname, "h1D_", (useLogY ? "cLogY_" : "c_"));
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
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (size_t is=0; is<sampleList.size(); is++){
      auto& sample = sampleList.at(is);
      TH1F* hist = sample_hist_map[sample.name].at(iplot);

      std::vector<TString> tmplist;
      TString htitle = hist->GetTitle();
      HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
      TString hlabel = tmplist.front();

      hist->SetTitle("");
      if (sample.name=="DY_2l_M_50_HT") hlabel = "DY (HT-binned)";
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
  MELAout.close();
}
