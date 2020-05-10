#include "common_includes.h"
#include "offshell_cutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <DataFormats/MuonReco/interface/Muon.h>
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


bool testTagBaseSelection(ElectronObject const* part){
  return ParticleSelectionHelpers::isTightParticle(part);
}
bool testExtraTightTagSelection(ElectronObject const* part){
  //return (part->extras.id_cutBased_Fall17V2_Tight_Bits == 1023);
  return part->extras.id_MVA_Fall17V2_NoIso_pass_wp80;
}
bool testPreselectionId(ElectronObject const* part){
  return part->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_id);
}
bool testPreselectionIso(ElectronObject const* part){
  return part->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_iso);
}

bool testTagBaseSelection(MuonObject const* part){
  return ParticleSelectionHelpers::isTightParticle(part);
}
bool testExtraTightTagSelection(MuonObject const* part){
  return ((part->extras.POG_selector_bits & reco::Muon::CutBasedIdTight) == reco::Muon::CutBasedIdTight);
}
bool testPreselectionId(MuonObject const* part){
  return part->testSelectionBit(MuonSelectionHelpers::bit_preselectionTight_id);
}
bool testPreselectionIso(MuonObject const* part){
  return part->testSelectionBit(MuonSelectionHelpers::bit_preselectionTight_iso);
}

bool testTagBaseSelection(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testTagBaseSelection(muon);
  else if (electron) return testTagBaseSelection(electron);
  else return false;
}
bool testExtraTightTagSelection(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testExtraTightTagSelection(muon);
  else if (electron) return testExtraTightTagSelection(electron);
  else return false;
}
bool testPreselectionId(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testPreselectionId(muon);
  else if (electron) return testPreselectionId(electron);
  else return false;
}
bool testPreselectionIso(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testPreselectionIso(muon);
  else if (electron) return testPreselectionIso(electron);
  else return false;
}

float get_dxy(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.dxy_bestTrack_firstPV;
  else if (electron) return electron->extras.dxy_firstPV;
  else return 0;
}
float get_dz(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.dz_bestTrack_firstPV;
  else if (electron) return electron->extras.dz_firstPV;
  else return 0;
}
float get_etaSC(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (electron) return electron->etaSC();
  else return part->eta();
}
cms3_egamma_fid_type_mask_t get_fid_mask(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (electron) return electron->extras.fid_mask;
  else return 0;
}
bool get_tightCharge(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.pass_tightCharge;
  else if (electron) return HelperFunctions::test_bit(electron->extras.charge_consistency_bits, 2);
  else return 0;
}
bool testTiming(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  if (muon) return muon->extras.pass_muon_timing;
  else return true;
}

using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false
){
  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main =
    "output/LeptonEfficiencies/SkimTrees/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);
  ParticleSelectionHelpers::setUseProbeLeptonsInLooseSelection(true);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btag_loose_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{
    TriggerHelpers::kSingleMu,
    TriggerHelpers::kSingleEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_highpt{
    TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle_HighPt
  };
  // The difference between control and prescaled is that the HLT prescale for the 'prescaled' label could be 0 or 1 (hopefully).
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_prescaled{
    TriggerHelpers::kSingleMu_Prescaled
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_control{
    TriggerHelpers::kSingleMu_Control,
    TriggerHelpers::kSingleEle_Control
  };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);
  auto triggerPropsCheckList_highpt = TriggerHelpers::getHLTMenuProperties(requiredTriggers_highpt);
  auto triggerPropsCheckList_prescaled = TriggerHelpers::getHLTMenuProperties(requiredTriggers_prescaled);
  auto triggerPropsCheckList_control = TriggerHelpers::getHLTMenuProperties(requiredTriggers_control);

  TString strSystName = SystematicsHelpers::getSystName(theGlobalSyst).data();
  SystematicsHelpers::SystematicVariationTypes eleSyst = theGlobalSyst;
  SystematicsHelpers::SystematicVariationTypes muSyst = theGlobalSyst;
  switch (theGlobalSyst){
  case SystematicsHelpers::eEleScaleDn:
  case SystematicsHelpers::eMuScaleDn:
    eleSyst = SystematicsHelpers::eEleScaleDn;
    muSyst = SystematicsHelpers::eMuScaleDn;
    strSystName = "LepScaleDn";
    break;
  case SystematicsHelpers::eEleScaleUp:
  case SystematicsHelpers::eMuScaleUp:
    eleSyst = SystematicsHelpers::eEleScaleUp;
    muSyst = SystematicsHelpers::eMuScaleUp;
    strSystName = "LepScaleUp";
    break;
  case SystematicsHelpers::eEleResDn:
  case SystematicsHelpers::eMuResDn:
    eleSyst = SystematicsHelpers::eEleResDn;
    muSyst = SystematicsHelpers::eMuResDn;
    strSystName = "LepResDn";
    break;
  case SystematicsHelpers::eEleResUp:
  case SystematicsHelpers::eMuResUp:
    eleSyst = SystematicsHelpers::eEleResUp;
    muSyst = SystematicsHelpers::eMuResUp;
    strSystName = "LepResUp";
    break;
  default:
    break;
  }

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = false;
  bool isQCD = false;
  float pTG_true_range[2]={ -1, -1 }; bool haspTGRange = false;
  for (auto const& sname:sampledirs){
    isData = SampleHelpers::checkSampleIsData(sname);
    if (isData && theGlobalSyst!=sNominal) return;
    if (isData && nchunks>0) return;

    isQCD = sname.Contains("QCD") && sname.Contains("HT");
    if ((sname.Contains("ZGTo2NuG") || sname.Contains("ZGTo2LG")) && sname.Contains("amcatnloFXFX") && !sname.Contains("PtG-130")) pTG_true_range[1]=130;

    break;
  }
  haspTGRange = pTG_true_range[0]!=pTG_true_range[1];
  constexpr bool needGenParticleChecks = true; // Always turned on because we need to do gen. matching

  gSystem->mkdir(coutput_main, true);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;

  PhotonScaleFactorHandler photonSFHandler;
  BtagScaleFactorHandler btagSFHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(needGenParticleChecks);

  eventFilter.setCheckTriggerObjectsForHLTPaths(true);
  //eventFilter.setTrackTriggerObjects(true);
  //eventFilter.setVerbosity(TVar::DEBUG);

  bool isFirstInputFile=true;
  for (auto const& sname:sampledirs){
    TString coutput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(coutput, "_MINIAOD", "");

    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/Dilepton", "cms3ntuple/Dilepton_Control", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    const int nEntries = sample_tree.getNEvents();

    double sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      bool hasCounters = true;
      {
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
            hasCounters = false;
            sum_wgts = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          ftmp->Close();
        }
        if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters){
        MELAerr << "This script is designed to use skim ntuples. Aborting..." << endl;
        return;
      }
    }
    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;

    // Set data tracking options
    eventFilter.setTrackDataEvents(isData);
    eventFilter.setCheckUniqueDataEvent(isData && !isFirstInputFile);

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

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += Form("_%s", strSystName.Data());
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree* tout = new TTree("SkimTree", "");
    MELAout << "Created output file " << stroutput << "..." << endl;

#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
    // Event variables
    BRANCH_COMMAND(float, event_wgt);
    BRANCH_COMMAND(float, event_wgt_triggers);
    BRANCH_COMMAND(float, event_wgt_SFs);
    BRANCH_COMMAND(float, pfmet_pTmiss); BRANCH_COMMAND(float, pfmet_phimiss);
    BRANCH_COMMAND(float, puppimet_pTmiss); BRANCH_COMMAND(float, puppimet_phimiss);
    BRANCH_COMMAND(bool, isNominalTrigger);
    BRANCH_COMMAND(bool, isHighPtTrigger);
    BRANCH_COMMAND(bool, isPrescaledTrigger);
    BRANCH_COMMAND(bool, isControlTrigger);
    BRANCH_COMMAND(unsigned int, event_nvtxs_good);
    BRANCH_COMMAND(unsigned int, event_Njets);
    BRANCH_COMMAND(unsigned int, event_Njets20);
    BRANCH_COMMAND(unsigned int, event_Njets_btagged);
    BRANCH_COMMAND(unsigned int, event_Njets20_btagged);
    // LL
    BRANCH_COMMAND(float, pt_ll);
    BRANCH_COMMAND(float, eta_ll);
    BRANCH_COMMAND(float, phi_ll);
    BRANCH_COMMAND(float, mass_ll);
    // L1: Tag
    BRANCH_COMMAND(int, id_l1);
    BRANCH_COMMAND(float, pt_l1);
    BRANCH_COMMAND(float, eta_l1);
    BRANCH_COMMAND(float, phi_l1);
    BRANCH_COMMAND(bool, pass_extraTight_l1);
    BRANCH_COMMAND(float, dxy_l1);
    BRANCH_COMMAND(float, dz_l1);
    BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l1);
    BRANCH_COMMAND(float, etaSC_l1);
    BRANCH_COMMAND(bool, isGenMatched_l1);
    BRANCH_COMMAND(bool, hasTightCharge_l1);
    BRANCH_COMMAND(bool, passTiming_l1);
    // L2: Probe
    BRANCH_COMMAND(int, id_l2);
    BRANCH_COMMAND(float, pt_l2);
    BRANCH_COMMAND(float, eta_l2);
    BRANCH_COMMAND(float, phi_l2);
    BRANCH_COMMAND(bool, pass_preselectionId_l2);
    BRANCH_COMMAND(bool, pass_preselectionIso_l2);
    BRANCH_COMMAND(float, dxy_l2);
    BRANCH_COMMAND(float, dz_l2);
    BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l2);
    BRANCH_COMMAND(float, etaSC_l2);
    BRANCH_COMMAND(bool, isGenMatched_l2);
    BRANCH_COMMAND(bool, hasTightCharge_l2);
    BRANCH_COMMAND(bool, passTiming_l2);
#undef BRANCH_COMMAND

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
    size_t n_pass_genWeights=0;
    std::vector<size_t> n_pass_dilepton(2, 0);
    size_t n_pass_photonVeto=0;
    size_t n_pass_isotrackVeto=0;
    size_t n_pass_uniqueEvent=0;
    size_t n_pass_commonFilters=0;
    size_t n_pass_goodPVFilter=0;
    std::vector<size_t> n_pass_triggers(2, 0);
    std::vector<size_t> n_hasTO(2, 0);
    std::vector<size_t> n_pass_triggercheck(2, 0);
    std::vector<size_t> n_pass_TOMatch(2, 0);
    std::vector<size_t> n_pass_hasTag(2, 0);
    std::vector<size_t> n_pass_tagTOMatch(2, 0);
    size_t n_pass_HEMfilter=0;
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);

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
        auto const& genparticles = genInfoHandler.getGenParticles();

        if (needGenParticleChecks){
          if (isQCD){
            auto const& genparticles = genInfoHandler.getGenParticles();
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
                event_wgt = 0;
                break;
              }
            }
          }
          if (haspTGRange){
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
                if ((pTG_true_range[0]>=0.f && part->pt()<pTG_true_range[0]) || (pTG_true_range[1]>=0.f && part->pt()>=pTG_true_range[1])) event_wgt = 0;
                break;
              }
            }
          }
        }

        simEventHandler.constructSimEvent(theGlobalSyst);
        event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();

        if (event_wgt==0.f) continue;
      }
      n_pass_genWeights++;

      muonHandler.constructMuons(muSyst);
      electronHandler.constructElectrons(eleSyst);
      photonHandler.constructPhotons(theGlobalSyst);
      particleDisambiguator.disambiguateParticles(&muonHandler, &electronHandler, &photonHandler);

      ParticleObject::LorentzVector_t sump4_leptons(0, 0, 0, 0);
      auto const& muons = muonHandler.getProducts();
      std::vector<MuonObject*> muons_selected; muons_selected.reserve(muons.size());
      for (auto const& part:muons){
        if (part->testSelectionBit(MuonSelectionHelpers::kProbeId)){ muons_selected.push_back(part); sump4_leptons += part->p4(); }
        if (part->extras.pass_muon_timing ^ ((part->extras.POG_selector_bits & reco::Muon::InTimeMuon) == reco::Muon::InTimeMuon)) MELAerr << "Inconsistent timing flags..." << endl;
      }

      auto const& electrons = electronHandler.getProducts();
      std::vector<ElectronObject*> electrons_selected; electrons_selected.reserve(electrons.size());
      for (auto const& part:electrons){
        if (part->testSelectionBit(ElectronSelectionHelpers::kProbeId)){ electrons_selected.push_back(part); sump4_leptons += part->p4(); }
      }
      if (
        !(
          (muons_selected.size()==2 && electrons_selected.empty() && muons_selected.front()->pdgId() != muons_selected.back()->pdgId())
          ||
          (electrons_selected.size()==2 && muons_selected.empty() && electrons_selected.front()->pdgId() != electrons_selected.back()->pdgId())
          )
        ) continue;
      if (std::abs(sump4_leptons.M() - PDGHelpers::Zmass)>=30.f || sump4_leptons.Pt()>=10000.) continue;

      unsigned short const idx_emu = 0*(muons_selected.size()==2) + 1*(electrons_selected.size()==2);
      n_pass_dilepton[idx_emu]++;

      auto const& photons = photonHandler.getProducts();
      unsigned int n_photons_tight = 0;
      float SF_photons = 1;
      PhotonObject const* theChosenPhoton = nullptr;
      for (auto const& part:photons){
        float theSF = 1;
        if (!isData) photonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_photons *= theSF;

        if (ParticleSelectionHelpers::isTightParticle(part)){
          if (!theChosenPhoton) theChosenPhoton = part;
          n_photons_tight++;
        }
      }
      if (n_photons_tight!=0) continue;
      n_pass_photonVeto++;

      event_wgt_SFs = SF_photons;

      isotrackHandler.constructIsotracks(&muons_selected, &electrons_selected);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;
      n_pass_isotrackVeto++;

      jetHandler.constructJetMET(theGlobalSyst, &muons_selected, &electrons_selected, &photons);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      auto pfmet_p4 = pfmet->p4(true, true, true);
      pfmet_pTmiss = pfmet_p4.Pt();
      pfmet_phimiss = pfmet_p4.Phi();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto puppimet_p4 = puppimet->p4(true, true, true);
      puppimet_pTmiss = puppimet_p4.Pt();
      puppimet_phimiss = puppimet_p4.Phi();

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      n_pass_uniqueEvent++;

      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters() || !eventFilter.hasGoodVertex()) continue;
      n_pass_commonFilters++;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;
      n_pass_goodPVFilter++;

      isNominalTrigger = isHighPtTrigger = isPrescaledTrigger = isControlTrigger = false;
      HLTTriggerPathObject const* passingHLTPath = nullptr;
      event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList, &muons_selected, &electrons_selected, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
      if (event_wgt_triggers != 0.f) isNominalTrigger = true;
      else{
        event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_highpt, &muons_selected, &electrons_selected, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
        if (event_wgt_triggers != 0.f) isHighPtTrigger = true;
        else{
          event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_prescaled, &muons_selected, &electrons_selected, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
          if (event_wgt_triggers != 0.f) isPrescaledTrigger = true;
          else{
            // Pass the ak4 jets as well because the electron control triggers also require PFJet30.
            event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_control, &muons_selected, &electrons_selected, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
            if (event_wgt_triggers != 0.f) isControlTrigger = true;
          }
        }
      }
      if (!(isNominalTrigger || isPrescaledTrigger || isControlTrigger)) continue;
      n_pass_triggers[idx_emu]++;

      std::vector<MuonObject*> muons_tag;
      std::vector<ElectronObject*> electrons_tag;
      {
        auto const& passedTriggerObjects = passingHLTPath->getPassedTriggerObjects();
        if (!passedTriggerObjects.empty()) n_hasTO[idx_emu]++;

        std::vector<MuonObject*> muons_trigcheck_TOmatched;
        std::vector<MuonObject*> muons_trigcheck; muons_trigcheck.reserve(muons_selected.size());
        for (auto const& part:muons_selected){ if (ParticleSelectionHelpers::isParticleForTriggerChecking(part)) muons_trigcheck.push_back(part); }

        std::vector<ElectronObject*> electrons_trigcheck_TOmatched;
        std::vector<ElectronObject*> electrons_trigcheck; electrons_trigcheck.reserve(electrons_selected.size());
        for (auto const& part:electrons_selected){ if (ParticleSelectionHelpers::isParticleForTriggerChecking(part)) electrons_trigcheck.push_back(part); }

        if (muons_trigcheck.empty() && electrons_trigcheck.empty()) continue;
        n_pass_triggercheck[idx_emu]++;

        TriggerObject::getMatchedPhysicsObjects(
          passedTriggerObjects, { trigger::TriggerMuon }, 0.2,
          muons_trigcheck, muons_trigcheck_TOmatched
        );
        TriggerObject::getMatchedPhysicsObjects(
          passedTriggerObjects, { trigger::TriggerElectron, trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
          electrons_trigcheck, electrons_trigcheck_TOmatched
        );

        if (muons_trigcheck_TOmatched.empty() && electrons_trigcheck_TOmatched.empty()) continue;
        n_pass_TOMatch[idx_emu]++;

        for (auto const& part:muons_trigcheck_TOmatched){ if (testTagBaseSelection(part)) muons_tag.push_back(part); }
        for (auto const& part:electrons_trigcheck_TOmatched){ if (testTagBaseSelection(part)) electrons_tag.push_back(part); }
      }
      if (muons_tag.empty() && electrons_tag.empty()) continue;
      n_pass_hasTag[idx_emu]++;

      std::vector<MuonObject*> muons_probe; muons_probe.reserve(muons_selected.size() - muons_tag.size());
      for (auto const& part:muons_selected){ if (!HelperFunctions::checkListVariable(muons_tag, part)) muons_probe.push_back(part); }
      std::vector<ElectronObject*> electrons_probe; electrons_probe.reserve(electrons_selected.size() - electrons_tag.size());
      for (auto const& part:electrons_selected){ if (!HelperFunctions::checkListVariable(electrons_tag, part)) electrons_probe.push_back(part); }

      // Test dilepton OS with trigger matching
      ParticleObject const* lepton_tag = nullptr;
      ParticleObject* lepton_probe = nullptr;
      if (muons_tag.size()==2){
        lepton_tag = (muons_tag.front()->pt()>=muons_tag.back()->pt() ? muons_tag.front() : muons_tag.back());
        lepton_probe = (muons_tag.front()->pt()>=muons_tag.back()->pt() ? muons_tag.back() : muons_tag.front());
      }
      else if (muons_tag.size()==1){
        lepton_tag = muons_tag.front();
        lepton_probe = muons_probe.front();
      }
      if (electrons_tag.size()==2){
        lepton_tag = (electrons_tag.front()->pt()>=electrons_tag.back()->pt() ? electrons_tag.front() : electrons_tag.back());
        lepton_probe = (electrons_tag.front()->pt()>=electrons_tag.back()->pt() ? electrons_tag.back() : electrons_tag.front());
      }
      else if (electrons_tag.size()==1){
        lepton_tag = electrons_tag.front();
        lepton_probe = electrons_probe.front();
      }
      if (!lepton_tag || !lepton_probe) continue;
      n_pass_tagTOMatch[idx_emu]++;

      // LL
      auto const p4_dilepton = lepton_tag->p4() + lepton_probe->p4();
      pt_ll = p4_dilepton.Pt();
      eta_ll = p4_dilepton.Eta();
      phi_ll = p4_dilepton.Phi();
      mass_ll = p4_dilepton.M();
      // L1
      id_l1 = lepton_tag->pdgId();
      pt_l1 = lepton_tag->pt();
      eta_l1 = lepton_tag->eta();
      phi_l1 = lepton_tag->phi();
      pass_extraTight_l1 = testExtraTightTagSelection(lepton_tag);
      dxy_l1 = get_dxy(lepton_tag);
      dz_l1 = get_dz(lepton_tag);
      fid_mask_l1 = get_fid_mask(lepton_tag);
      etaSC_l1 = get_etaSC(lepton_tag);
      hasTightCharge_l1 = get_tightCharge(lepton_tag);
      passTiming_l1 = testTiming(lepton_tag);
      // L2
      id_l2 = lepton_probe->pdgId();
      pt_l2 = lepton_probe->pt();
      eta_l2 = lepton_probe->eta();
      phi_l2 = lepton_probe->phi();
      pass_preselectionId_l2 = testPreselectionId(lepton_probe);
      pass_preselectionIso_l2 = testPreselectionIso(lepton_probe);
      dxy_l2 = get_dxy(lepton_probe);
      dz_l2 = get_dz(lepton_probe);
      fid_mask_l2 = get_fid_mask(lepton_probe);
      etaSC_l2 = get_etaSC(lepton_probe);
      hasTightCharge_l2 = get_tightCharge(lepton_probe);
      passTiming_l2 = testTiming(lepton_probe);

      std::vector<ElectronObject*> electron_probe_container; electron_probe_container.reserve(1);
      if (std::abs(id_l2)==11) electron_probe_container.push_back(dynamic_cast<ElectronObject*>(lepton_probe));

      // Test HEM filter
      if (!eventFilter.test2018HEMFilter(&simEventHandler, &electron_probe_container, nullptr, &ak4jets, &ak8jets)) continue;
      n_pass_HEMfilter++;

      event_Njets = 0;
      event_Njets20 = 0;
      event_Njets_btagged = 0;
      event_Njets20_btagged = 0;
      float SF_btagging = 1;
      ParticleObject::LorentzVector_t ak4jets_sump4(0, 0, 0, 0);
      for (auto* jet:ak4jets){
        float theSF = 1;
        if (!isData) btagSFHandler.getSFAndEff(theGlobalSyst, jet, theSF, nullptr);
        if (theSF != 0.f) SF_btagging *= theSF;

        if (ParticleSelectionHelpers::isTightJet(jet)){
          event_Njets++;
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets_btagged++;
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
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets20_btagged++;
        }
      }
      event_wgt_SFs *= SF_btagging;

      sample_tree.getVal("vtxs_nvtxs_good", event_nvtxs_good);

      isGenMatched_l1 = isGenMatched_l2 = false;
      if (!isData){
        std::vector<ParticleObject const*> leptons; leptons.reserve(2);
        leptons.push_back(lepton_tag);
        leptons.push_back(lepton_probe);
        auto const& genparticles = genInfoHandler.getGenParticles();
        std::vector<GenParticleObject const*> genpromptleptons; genpromptleptons.reserve(genparticles.size());
        for (auto const& part:genparticles){
          if (part->extras.isPromptFinalState && PDGHelpers::isALepton(part->pdgId())) genpromptleptons.push_back(part);
        }
        std::unordered_map<ParticleObject const*, GenParticleObject const*> tmp_map;
        ParticleObjectHelpers::matchParticles(
          ParticleObjectHelpers::kMatchBy_DeltaR, 0.2,
          leptons.cbegin(), leptons.cend(),
          genpromptleptons.cbegin(), genpromptleptons.cend(),
          tmp_map
        );

        auto it_tag = tmp_map.find(leptons.front());
        if (it_tag!=tmp_map.cend() && it_tag->second) isGenMatched_l1 = (std::abs(it_tag->first->pdgId()) == std::abs(it_tag->second->pdgId()));

        auto it_probe = tmp_map.find(leptons.back());
        if (it_probe!=tmp_map.cend() && it_probe->second) isGenMatched_l2 = (std::abs(it_probe->first->pdgId()) == std::abs(it_probe->second->pdgId()));
      }

      tout->Fill();
      n_evts_acc++;
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc << " / " << (ev_end - ev_start) << endl;
    MELAout << "\t- Number of events passing each cut:\n"
      << "\t\t- Gen. weights!=0: " << n_pass_genWeights << '\n'
      << "\t\t- Dilepton base selection: " <<  n_pass_dilepton[0]+n_pass_dilepton[1] << " (" << n_pass_dilepton << ")" << '\n'
      << "\t\t- Photon veto: " <<  n_pass_photonVeto << '\n'
      << "\t\t- Isotrack veto: " <<  n_pass_isotrackVeto << '\n'
      << "\t\t- Unique event: " <<  n_pass_uniqueEvent << '\n'
      << "\t\t- Common filters: " <<  n_pass_commonFilters << '\n'
      << "\t\t- Good PV filter: " << n_pass_goodPVFilter << '\n'
      << "\t\t- Trigger: " <<  n_pass_triggers[0]+n_pass_triggers[1] << " (" << n_pass_triggers << ")" << '\n'
      << "\t\t- Has trigger objects: " <<  n_hasTO[0]+n_hasTO[1] << " (" << n_hasTO << ")" << '\n'
      << "\t\t- Trigger check: " <<  n_pass_triggercheck[0]+n_pass_triggercheck[1] << " (" << n_pass_triggercheck << ")" << '\n'
      << "\t\t- Pass TO matching: " <<  n_pass_TOMatch[0]+n_pass_TOMatch[1] << " (" << n_pass_TOMatch << ")" << '\n'
      << "\t\t- Has tag: " << n_pass_hasTag[0]+n_pass_hasTag[1] << " (" << n_pass_hasTag << ")" << '\n'
      << "\t\t- Tag trigger object matching: " <<  n_pass_tagTOMatch[0]+n_pass_tagTOMatch[1] << " (" << n_pass_tagTOMatch << ")" << '\n'
      << "\t\t- HEM15/16 veto: " << n_pass_HEMfilter
      << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;

    // Write output
    foutput->WriteTObject(tout);
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples list
}
