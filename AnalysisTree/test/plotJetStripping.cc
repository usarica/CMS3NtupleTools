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
  std::vector<TH1F*> const& hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts="hist",
  bool useLogY=true,
  bool adjustYLow=false,
  float factorYHigh=-1
);


using namespace SystematicsHelpers;
void plot(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool useDilepton = false, bool use_ee=false, bool forceInvertPTmissThr=false
){
  if (forceInvertPTmissThr && useDilepton) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  bool const periodIsLikeData = SampleHelpers::testDataPeriodIsLikeData();

  // Jet ID options
  bool applyPUIdToAK4Jets=true, applyTightLeptonVetoIdToAK4Jets=false;
  // MET options
  bool use_MET_Puppi=false;
  bool use_MET_XYCorr=true, use_MET_JERCorr=false, use_MET_ParticleMomCorr=true;
  bool use_MET_p4Preservation=false;
  bool use_MET_corrections=false;
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TDirectory* curdir = gDirectory;

  bool isData = false;
  {
    std::vector<TString> sampledirs;
    SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
    isData = SampleHelpers::checkSampleIsData(sampledirs.front());
  }
  if (isData && theGlobalSyst!=sNominal) return;
  bool useInvertedPTmissReq = (isData && useDilepton) || forceInvertPTmissThr;

  TString cinput = SampleHelpers::getSampleIdentifier(strSampleSet);
  HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
  HelperFunctions::replaceString(cinput, "_MINIAOD", "");
  bool containsPeriodString = false;
  for (auto const& pp:validDataPeriods){
    if (cinput.Contains(pp)){ containsPeriodString = true; break; }
  }
  bool isGGSig = cinput.Contains("GluGluH");
  bool isVBFSig = cinput.Contains("VBF");
  float xsec_scale = 1;
  if (isGGSig) xsec_scale = 0.00813704*0.541;
  if (isVBFSig) xsec_scale = 0.00846586*0.541;

  TString cinput_base_dir;
  /*if (!SampleHelpers::checkRunOnCondor()) cinput_base_dir = "output/";
  else*/ cinput_base_dir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/";

  // Set input directory
  std::vector<TChain*> tinlist;
  for (unsigned short it=0; it<3; it++){
    bool const useJetOverlapStripping = (it==2);

    TString cinput_main =
      cinput_base_dir + (useDilepton ? "DileptonEvents" : "SinglePhotonEvents") + "/SkimTrees/" + strdate
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

    tinlist.push_back(new TChain("SkimTree"));
    if (isData || containsPeriodString){
      for (auto const& pp:validDataPeriods){
        if (periodIsLikeData && pp!=SampleHelpers::getDataPeriod()) continue;

        TString cinput_tmp = cinput;
        if (containsPeriodString) HelperFunctions::replaceString(cinput_tmp, SampleHelpers::getDataPeriod(), pp);
        else if (isData && cinput_tmp==Form("Run%i", SampleHelpers::getDataYear())) cinput_tmp = Form("Run%s", pp.Data());

        TString strinput = Form("%s/%s/%s", cinput_main.Data(), pp.Data(), cinput_tmp.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
        strinput += ".root";

        MELAout << "Adding input " << strinput << " for tree " << it << endl;
        tinlist.back()->Add(strinput);
      }
    }
    else{
      TString strinput = Form("%s/%s/%s", cinput_main.Data(), period.Data(), cinput.Data());
      strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
      strinput += ".root";

      MELAout << "Adding input " << strinput << " for tree " << it << endl;
      tinlist.back()->Add(strinput);
    }

    MELAout << "Data tree has a total of " << tinlist.back()->GetEntries() << " entries..." << endl;
    curdir->cd();
  }

  TString coutput_main = 
    TString("output/JetStrippingChecks/") + (useDilepton ? "DileptonEvents" : "SinglePhotonEvents") + "/" + strdate
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  coutput_main = coutput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  if (!isData) coutput_main = coutput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  coutput_main = coutput_main + "/" + period;
  if (useDilepton) coutput_main = coutput_main + "/" + (use_ee ? "ee" : "mumu");
  else if (forceInvertPTmissThr) coutput_main = coutput_main + "/InvetedPTMissThr";
  coutput_main = coutput_main + "/" + cinput;

  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + "/histograms.root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

#define BRANCH_OPTIONAL_COMMANDS \
  BRANCH_COMMAND(float, KFactor_QCD_NNLO_ggZZ_Sig_Nominal) \
  BRANCH_COMMAND(float, p_Gen_CPStoBWPropRewgt) \
  BRANCH_COMMAND(float, p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM) \
  BRANCH_COMMAND(float, p_Gen_JJEW_SIG_ghv1_1_MCFM)
#define BRANCH_SCALAR_COMMON_COMMANDS \
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
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_SCALAR_SINGLEPHOTON_COMMANDS \
  BRANCH_COMMAND(float, photon_pt) \
  BRANCH_COMMAND(float, photon_eta) \
  BRANCH_COMMAND(float, photon_phi) \
  BRANCH_COMMAND(float, photon_mass) \
  BRANCH_COMMAND(bool, photon_is_conversionSafe) \
  BRANCH_COMMAND(bool, photon_is_inTime) \
  BRANCH_COMMAND(bool, photon_is_beamHaloSafe) \
  BRANCH_COMMAND(bool, photon_is_spikeSafe) \
  BRANCH_COMMAND(bool, photon_is_PFID) \
  BRANCH_COMMAND(bool, photon_is_METSafe) \
  BRANCH_COMMAND(bool, photon_pass_HGGSelection) \
  BRANCH_COMMAND(bool, photon_isGap) \
  BRANCH_COMMAND(bool, photon_isEB) \
  BRANCH_COMMAND(bool, photon_isEE) \
  BRANCH_COMMAND(bool, photon_isEBEEGap)
#define BRANCH_SCALAR_DILEPTON_COMMANDS \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass)
#define BRANCH_VECTOR_COMMON_COMMANDS \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched) \
  BRANCH_COMMAND(bool, ak4jets_is_genMatched_fullCone) \
  BRANCH_COMMAND(cms3_listSize_t, ak4jets_n_overlaps) \
  BRANCH_COMMAND(float, ak4jets_overlaps_pt) \
  BRANCH_COMMAND(float, ak4jets_original_pt) \
  BRANCH_COMMAND(unsigned char, ak4jets_btagWP_Bits) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(cms3_listSize_t, ak4jets_masked_n_overlaps) \
  BRANCH_COMMAND(float, ak4jets_masked_overlaps_pt) \
  BRANCH_COMMAND(unsigned char, ak4jets_masked_btagWP_Bits) \
  BRANCH_COMMAND(float, ak4jets_masked_pt) \
  BRANCH_COMMAND(float, ak4jets_masked_eta) \
  BRANCH_COMMAND(float, ak4jets_masked_phi) \
  BRANCH_COMMAND(float, ak4jets_masked_mass)
#define BRANCH_VECTOR_SINGLEPHOTON_COMMANDS 
#define BRANCH_VECTOR_DILEPTON_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt)
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_SCALAR_COMMON_COMMANDS \
  BRANCH_SCALAR_SINGLEPHOTON_COMMANDS \
  BRANCH_SCALAR_DILEPTON_COMMANDS
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_VECTOR_COMMON_COMMANDS \
  BRANCH_VECTOR_DILEPTON_COMMANDS
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS
#define BRANCH_DILEPTON_COMMANDS \
  BRANCH_SCALAR_COMMON_COMMANDS \
  BRANCH_SCALAR_DILEPTON_COMMANDS \
  BRANCH_VECTOR_COMMON_COMMANDS \
  BRANCH_VECTOR_DILEPTON_COMMANDS
#define BRANCH_SINGLEPHOTON_COMMANDS \
  BRANCH_SCALAR_COMMON_COMMANDS \
  BRANCH_SCALAR_SINGLEPHOTON_COMMANDS \
  BRANCH_VECTOR_COMMON_COMMANDS \
  BRANCH_VECTOR_SINGLEPHOTON_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 1;
  BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND


#define HISTOGRAM_COMMANDS \
  HISTOGRAM_COMMAND(event_n_ak4jets_pt30, "N_{jets}", "Event / bin", "", 7, 0, 7) \
  HISTOGRAM_COMMAND(foverlap_masked, "f^{l}_{j+l}", "Jets (removed) / bin", "", 50, 0, 1) \
  HISTOGRAM_COMMAND(foverlap_stripped, "f^{l}_{j+l}", "Jets (stripped) / bin", "", 50, 0, 1) \
  HISTOGRAM_COMMAND(mTZZ_inclusive, "m_{T}^{ZZ} (GeV)", "Event / bin", "", (1500-150)/50, 150, 1500) \
  HISTOGRAM_COMMAND(mTZZ_Njets_eq_0, "m_{T}^{ZZ} (GeV)", "Event / bin", "N_{jets}=0", (1500-150)/50, 150, 1500) \
  HISTOGRAM_COMMAND(boson_pt_Njets_eq_0, "p_{T}^{boson} (GeV)", "Event / bin", "N_{jets}=0", 100, 0, 1000) \
  HISTOGRAM_COMMAND(event_pTmiss_Njets_eq_0, "p_{T}^{miss} (GeV)", "Event / bin", "N_{jets}=0", 100, 0, 1000) \
  HISTOGRAM_COMMAND(mTZZ_Njets_eq_1, "m_{T}^{ZZ} (GeV)", "Event / bin", "N_{jets}=1", (1500-150)/50, 150, 1500) \
  HISTOGRAM_COMMAND(boson_pt_Njets_eq_1, "p_{T}^{boson} (GeV)", "Event / bin", "N_{jets}=1", 100, 0, 1000) \
  HISTOGRAM_COMMAND(event_pTmiss_Njets_eq_1, "p_{T}^{miss} (GeV)", "Event / bin", "N_{jets}=1", 100, 0, 1000) \
  HISTOGRAM_COMMAND(ak4jets_pt_Njets_eq_1, "p_{T}^{jet} (GeV)", "Jets / bin", "N_{jets}=1", (500-30)/10, 30, 500) \
  HISTOGRAM_COMMAND(ak4jets_pt_abseta_lt_2p5_Njets_eq_1, "p_{T}^{jet} (GeV)", "Jets (|#eta_{jet}|<2.5) / bin", "N_{jets}=1", (500-30)/10, 30, 500) \
  HISTOGRAM_COMMAND(mTZZ_Njets_eq_2, "m_{T}^{ZZ} (GeV)", "Event / bin", "N_{jets}=2", (1500-150)/50, 150, 1500) \
  HISTOGRAM_COMMAND(boson_pt_Njets_eq_2, "p_{T}^{boson} (GeV)", "Event / bin", "N_{jets}=2", 100, 0, 1000) \
  HISTOGRAM_COMMAND(event_pTmiss_Njets_eq_2, "p_{T}^{miss} (GeV)", "Event / bin", "N_{jets}=2", 100, 0, 1000) \
  HISTOGRAM_COMMAND(ak4jets_pt_Njets_eq_2, "p_{T}^{jet} (GeV)", "Jets / bin", "N_{jets}=2", (500-30)/10, 30, 500) \
  HISTOGRAM_COMMAND(ak4jets_pt_abseta_lt_2p5_Njets_eq_2, "p_{T}^{jet} (GeV)", "Jets (|#eta_{jet}|<2.5) / bin", "N_{jets}=2", (500-30)/10, 30, 500) \
  HISTOGRAM_COMMAND(mTZZ_Njets_gt_2, "m_{T}^{ZZ} (GeV)", "Event / bin", "N_{jets}>2", (1500-150)/50, 150, 1500) \
  HISTOGRAM_COMMAND(boson_pt_Njets_gt_2, "p_{T}^{boson} (GeV)", "Event / bin", "N_{jets}>2", 100, 0, 1000) \
  HISTOGRAM_COMMAND(event_pTmiss_Njets_gt_2, "p_{T}^{miss} (GeV)", "Event / bin", "N_{jets}>2", 100, 0, 1000) \
  HISTOGRAM_COMMAND(ak4jets_pt_Njets_gt_2, "p_{T}^{jet} (GeV)", "Jets / bin", "N_{jets}>2", (500-30)/10, 30, 500) \
  HISTOGRAM_COMMAND(ak4jets_pt_abseta_lt_2p5_Njets_gt_2, "p_{T}^{jet} (GeV)", "Jets (|#eta_{jet}|<2.5) / bin", "N_{jets}>2", (500-30)/10, 30, 500)

  foutput->cd();
#define HISTOGRAM_COMMAND(NAME, XLABEL, YLABEL, CUTLABEL, NBINS, XLOW, XHIGH) std::vector<TH1F*> hlist_##NAME; hlist_##NAME.reserve(tinlist.size());
  HISTOGRAM_COMMANDS;
#undef HISTOGRAM_COMMAND
  for (unsigned int it=0; it<tinlist.size(); it++){
    int icol = -1;
    switch (it){
    case 0:
      icol = (int) kBlue;
      break;
    case 1:
      icol = (int) kRed;
      break;
    case 2:
      icol = (int) kViolet;
      break;
    default:
      break;
    }
#define HISTOGRAM_COMMAND(NAME, XLABEL, YLABEL, CUTLABEL, NBINS, XLOW, XHIGH) \
    hlist_##NAME.emplace_back(new TH1F(Form("%s__%i", #NAME, it), "", NBINS, XLOW, XHIGH)); \
    hlist_##NAME.back()->Sumw2(); \
    hlist_##NAME.back()->SetMarkerColor(icol); \
    hlist_##NAME.back()->SetLineColor(icol); \
    hlist_##NAME.back()->SetLineWidth(2); \
    hlist_##NAME.back()->GetXaxis()->SetTitle(XLABEL); \
    hlist_##NAME.back()->GetYaxis()->SetTitle(YLABEL);
    HISTOGRAM_COMMANDS;
#undef HISTOGRAM_COMMAND
  }

  {
    using namespace OffshellCutflow;

    unsigned int it=0;
    for (auto& tin:tinlist){
#define BRANCH_COMMAND(TYPE, NAME) bool exists_##NAME = tin->GetBranchStatus(#NAME)==1;
      BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND

      tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) if (exists_##NAME){ MELAout << "Booking " << #NAME << "..." << endl; tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME); }
      BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND

#define BRANCH_COMMAND(TYPE, NAME) MELAout << "Booking " << #NAME << "..." << endl; tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
      if (useDilepton){
        BRANCH_DILEPTON_COMMANDS;
      }
      else{
        BRANCH_SINGLEPHOTON_COMMANDS;
      }
#undef BRANCH_COMMAND

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
      double sum_wgts_cuts_mTZZ_geq_350[ncuts]={ 0 };
      double sum_wgts_total = 0;
      const int nEntries = tin->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tin->GetEvent(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (!useDilepton) event_n_ak4jets_pt30_btagged_loose=0;

        float wgt = event_wgt*event_wgt_triggers*event_wgt_SFs*xsec_scale;
#define BRANCH_COMMAND(TYPE, NAME) if (exists_##NAME) wgt *= NAME;
        BRANCH_OPTIONAL_COMMANDS;
#undef BRANCH_COMMAND
        sum_wgts_total += wgt;

        if (
          (
            useDilepton
            &&
            (
            (use_ee && dilepton_id!=-11*11)
              ||
              (!use_ee && dilepton_id!=-13*13)
              )
            )
          ||
          (
            !useDilepton
            &&
            !(
              photon_is_spikeSafe
              &&
              photon_is_conversionSafe
              &&
              photon_is_inTime
              &&
              photon_is_beamHaloSafe
              &&
              photon_is_PFID
              &&
              photon_is_METSafe
              )
            )
          ) continue;
        sum_wgts_cuts[0] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[0] += wgt;

        if (useDilepton && !(check_pTl1(leptons_pt->front()) && check_pTl2(leptons_pt->back()))) continue;
        if (!useDilepton && !photon_isEB) continue;
        if (useDilepton && !check_mll(dilepton_mass, (dilepton_id!=-11*13))) continue;
        sum_wgts_cuts[1] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[1] += wgt;

        if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
        sum_wgts_cuts[2] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[2] += wgt;

        // If this is tree set 1, modify the counts to include some of the masked jets
        std::vector<ParticleObject::LorentzVector_t> extrajets;
        if (it==1){
          ParticleObject::LorentzVector_t p4_boson;
          p4_boson = ParticleObject::PolarLorentzVector_t(
            (useDilepton ? dilepton_pt : photon_pt),
            (useDilepton ? dilepton_eta : photon_eta),
            (useDilepton ? dilepton_phi : photon_phi),
            (useDilepton ? dilepton_mass : photon_mass)
          );
          ParticleObject::LorentzVector_t p4_jets;

          for (unsigned int ijet=0; ijet<ak4jets_pt->size(); ijet++){
            ParticleObject::LorentzVector_t p4_jet;
            p4_jet = ParticleObject::PolarLorentzVector_t(
              ak4jets_pt->at(ijet),
              ak4jets_eta->at(ijet),
              ak4jets_phi->at(ijet),
              ak4jets_mass->at(ijet)
            );
            p4_jets += p4_jet;
          }
          for (unsigned int ijet=0; ijet<ak4jets_masked_pt->size(); ijet++){
            if (ak4jets_masked_overlaps_pt->at(ijet)/ak4jets_masked_pt->at(ijet)>=0.8) continue;

            ParticleObject::LorentzVector_t p4_jet;
            p4_jet = ParticleObject::PolarLorentzVector_t(
              ak4jets_masked_pt->at(ijet),
              ak4jets_masked_eta->at(ijet),
              ak4jets_masked_phi->at(ijet),
              ak4jets_masked_mass->at(ijet)
            );
            p4_jets += p4_jet;
            extrajets.push_back(p4_jet);

            float dphi_tmp;
            HelperFunctions::deltaPhi(float(p4_jet.phi()), event_phimiss, dphi_tmp); dphi_tmp = std::abs(dphi_tmp);
            min_abs_dPhi_pTj_pTmiss = std::min(min_abs_dPhi_pTj_pTmiss, dphi_tmp);

            if (HelperFunctions::test_bit(ak4jets_masked_btagWP_Bits->at(ijet), 0)) event_n_ak4jets_pt30_btagged_loose++;
            event_n_ak4jets_pt30++;
          }

          HelperFunctions::deltaPhi(float((p4_boson+p4_jets).phi()), event_phimiss, dPhi_pTbosonjets_pTmiss);
        }

        if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
        sum_wgts_cuts[3] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[3] += wgt;

        if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss)) continue;
        sum_wgts_cuts[4] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[4] += wgt;

        if (event_n_ak4jets_pt30_btagged_loose>0) continue;
        sum_wgts_cuts[5] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[5] += wgt;

        bool pass_pTmiss = check_pTmiss(event_pTmiss);
        if (isData && useDilepton && pass_pTmiss) continue; // Blinded
        if (useInvertedPTmissReq) pass_pTmiss = !pass_pTmiss;

        bool pass_boson_pt = check_pTboson((useDilepton ? dilepton_pt : photon_pt));

        unsigned int index_Njets = std::min((unsigned int) 3, event_n_ak4jets_pt30);

        TH1F* h_event_n_ak4jets_pt30 = hlist_event_n_ak4jets_pt30.at(it);
        TH1F* h_foverlap_masked = hlist_foverlap_masked.at(it);
        TH1F* h_foverlap_stripped = hlist_foverlap_stripped.at(it);
        TH1F* h_mTZZ_inclusive = hlist_mTZZ_inclusive.at(it);
        TH1F* h_mTZZ = nullptr;
        TH1F* h_boson_pt = nullptr;
        TH1F* h_pTmiss = nullptr;
        TH1F* h_ak4jets_pt = nullptr;
        TH1F* h_ak4jets_pt_abseta_lt_2p5 = nullptr;
        switch (index_Njets){
        case 0:
          h_mTZZ = hlist_mTZZ_Njets_eq_0.at(it);
          h_boson_pt = hlist_boson_pt_Njets_eq_0.at(it);
          h_pTmiss = hlist_event_pTmiss_Njets_eq_0.at(it);
          break;
        case 1:
          h_mTZZ = hlist_mTZZ_Njets_eq_1.at(it);
          h_boson_pt = hlist_boson_pt_Njets_eq_1.at(it);
          h_pTmiss = hlist_event_pTmiss_Njets_eq_1.at(it);
          h_ak4jets_pt = hlist_ak4jets_pt_Njets_eq_1.at(it);
          h_ak4jets_pt_abseta_lt_2p5 = hlist_ak4jets_pt_abseta_lt_2p5_Njets_eq_1.at(it);
          break;
        case 2:
          h_mTZZ = hlist_mTZZ_Njets_eq_2.at(it);
          h_boson_pt = hlist_boson_pt_Njets_eq_2.at(it);
          h_pTmiss = hlist_event_pTmiss_Njets_eq_2.at(it);
          h_ak4jets_pt = hlist_ak4jets_pt_Njets_eq_2.at(it);
          h_ak4jets_pt_abseta_lt_2p5 = hlist_ak4jets_pt_abseta_lt_2p5_Njets_eq_2.at(it);
          break;
        case 3:
          h_mTZZ = hlist_mTZZ_Njets_gt_2.at(it);
          h_boson_pt = hlist_boson_pt_Njets_gt_2.at(it);
          h_pTmiss = hlist_event_pTmiss_Njets_gt_2.at(it);
          h_ak4jets_pt = hlist_ak4jets_pt_Njets_gt_2.at(it);
          h_ak4jets_pt_abseta_lt_2p5 = hlist_ak4jets_pt_abseta_lt_2p5_Njets_gt_2.at(it);
          break;
        default:
          MELAerr << "index_Njets = " << index_Njets << " is undefined!" << endl;
          assert(0);
          break;
        }

        if (pass_pTmiss) h_boson_pt->Fill((useDilepton ? dilepton_pt : photon_pt), wgt);
        if (pass_boson_pt) h_pTmiss->Fill(event_pTmiss, wgt);

        if (!pass_pTmiss || !pass_boson_pt) continue;

        h_event_n_ak4jets_pt30->Fill(event_n_ak4jets_pt30, wgt);
        h_mTZZ->Fill(event_mTZZ, wgt);
        h_mTZZ_inclusive->Fill(event_mTZZ, wgt);
        if (it!=1){
          for (unsigned int ijet=0; ijet<ak4jets_masked_pt->size(); ijet++){
            h_foverlap_masked->Fill(ak4jets_masked_overlaps_pt->at(ijet)/ak4jets_masked_pt->at(ijet), wgt);
          }
        }
        if (it==2){
          for (unsigned int ijet=0; ijet<ak4jets_pt->size(); ijet++){
            h_foverlap_stripped->Fill(ak4jets_overlaps_pt->at(ijet)/ak4jets_original_pt->at(ijet), wgt);
          }
        }

        if (index_Njets>=1){
          for (unsigned int ijet=0; ijet<ak4jets_pt->size(); ijet++){
            h_ak4jets_pt->Fill(ak4jets_pt->at(ijet), wgt);
            if (std::abs(ak4jets_eta->at(ijet))<2.5) h_ak4jets_pt_abseta_lt_2p5->Fill(ak4jets_pt->at(ijet), wgt);
          }
          for (auto const& jp4:extrajets) h_ak4jets_pt->Fill(jp4.Pt(), wgt);
        }

        sum_wgts_cuts[6] += wgt;
        if (event_mTZZ>=350.f) sum_wgts_cuts_mTZZ_geq_350[6] += wgt;
      }
      MELAout << "Accumulated sum of weights:" << sum_wgts_total << endl;
      for (unsigned int icut=0; icut<ncuts; icut++){
        MELAout << "\t- " << cutlabels.at(icut) << ": " << sum_wgts_cuts[icut]
          << ", efficiency: " << sum_wgts_cuts[icut]/sum_wgts_cuts[0]
          << ",  recursive eff.: " << sum_wgts_cuts[icut]/sum_wgts_cuts[(icut==0 ? 0 : icut-1)]
          << ",  mTZZ>=350 GeV: " << sum_wgts_cuts_mTZZ_geq_350[icut]
          << endl;
      }

      delete tin;
      foutput->cd();
      it++;
    }
  }

  std::vector<TString> hlabels{
    "#DeltaR<0.4 cleaning",
    "#DeltaR<0.4, f^{l}_{j+l}#geq0.8 cleaning",
    "#DeltaR<0.4, PF cand. stripping"
  };
  assert(hlabels.size() == tinlist.size());
  TString strChannelLabel = (useDilepton ? (use_ee ? "ee channel" : "#mu#mu channel") : "Single #gamma");
  TString strChannelTitle = (useDilepton ? (use_ee ? "ee" : "mumu") : "SinglePhoton");

  TString selectionLabels;
  TString label_Nj;

#define HISTOGRAM_COMMAND(NAME, XLABEL, YLABEL, CUTLABEL, NBINS, XLOW, XHIGH) \
  std::vector<TH1F*> hplot_##NAME; hplot_##NAME.reserve(hlist_##NAME.size()); \
  std::vector<TString> hlabels_##NAME; hlabels_##NAME.reserve(hlist_##NAME.size()); \
  for (unsigned short it=0; it<hlist_##NAME.size(); it++){ \
    if (it==1 && TString(#NAME).Contains("foverlap_masked")) continue; \
    if (it!=2 && TString(#NAME).Contains("foverlap_stripped")) continue; \
    foutput->WriteTObject(hlist_##NAME.at(it)); \
    hplot_##NAME.push_back(hlist_##NAME.at(it)); \
    hlabels_##NAME.push_back(hlabels.at(it)); \
  } \
  if (!useInvertedPTmissReq){ \
    if (TString(#NAME).Contains("boson_pt")) selectionLabels = Form("%s|All selections except p_{T}^{boson}|%s", strChannelLabel.Data(), CUTLABEL); \
    else if (TString(#NAME).Contains("pTmiss")) selectionLabels = Form("%s|All selections except p_{T}^{miss}|%s", strChannelLabel.Data(), CUTLABEL); \
    else selectionLabels = Form("%s|All selections|%s", strChannelLabel.Data(), CUTLABEL); \
  } \
  else{ \
    if (TString(#NAME).Contains("boson_pt")) selectionLabels = Form("%s|All selections except p_{T}^{boson}, and p_{T}^{miss}<125 GeV|%s", strChannelLabel.Data(), CUTLABEL); \
    else if (TString(#NAME).Contains("pTmiss")) selectionLabels = Form("%s|All selections but p_{T}^{miss}<125 GeV|%s", strChannelLabel.Data(), CUTLABEL); \
    else selectionLabels = Form("%s|All selections but p_{T}^{miss}<125 GeV|%s", strChannelLabel.Data(), CUTLABEL); \
  } \
  makePlot( \
  coutput_main, lumi, \
    Form("c_%s_%s", #NAME, strChannelTitle.Data()), \
    hplot_##NAME, hlabels_##NAME, \
    selectionLabels, \
    "e1p", true, true, 15 \
    );

  HISTOGRAM_COMMANDS;
#undef HISTOGRAM_COMMAND

  foutput->Close();

#undef HISTOGRAM_COMMANDS
}

void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> const& hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool useLogY,
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
  if (adjustYLow || useLogY) ymin = 9e9;
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
      if ((adjustYLow || useLogY) && (!useLogY || bclow>0.)) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //MELAout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : (!useLogY ? 1.5 : 15.));
  else ymax /= (factorYHigh>0.f ? factorYHigh : (!useLogY ? 1.5 : 15.));
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (auto& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

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
  if (useLogY) canvas->SetLogy(true);

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
