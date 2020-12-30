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
  srcDir(srcDir_),
  hist(name, title, nbins, xlow, xhigh)
{
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
  for (unsigned int imet=0; imet<3; imet++){
    bool do_met_low, do_met_high;
    float met_low, met_high;
    if (imet==0){
      do_met_low=false; do_met_high=true;
      met_low=-1; met_high=20;
    }
    else if (imet==1){
      do_met_low=true; do_met_high=true;
      met_low=20; met_high=125;
    }
    else{
      do_met_low=true; do_met_high=false;
      met_low=125; met_high=-1;
    }
    // Nj==0
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(4);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 0, 0
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    cutsets.back().emplace_back(
      "pTl1", "p_{T}^{l1}",
      true, false, 25, 0
    );
    cutsets.back().emplace_back(
      "pTl2", "p_{T}^{l2}",
      true, false, 20, 0
    );
    // Nj==1
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(4);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 1, 1
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    cutsets.back().emplace_back(
      "pTl1", "p_{T}^{l1}",
      true, false, 25, 0
    );
    cutsets.back().emplace_back(
      "pTl2", "p_{T}^{l2}",
      true, false, 20, 0
    );
    // Nj==2
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(4);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 2, 2
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    cutsets.back().emplace_back(
      "pTl1", "p_{T}^{l1}",
      true, false, 25, 0
    );
    cutsets.back().emplace_back(
      "pTl2", "p_{T}^{l2}",
      true, false, 20, 0
    );
    // Nj>=3
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(4);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, false, 3, -1
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    cutsets.back().emplace_back(
      "pTl1", "p_{T}^{l1}",
      true, false, 25, 0
    );
    cutsets.back().emplace_back(
      "pTl2", "p_{T}^{l2}",
      true, false, 20, 0
    );
  }
}



void getTrees(int procsel, int ichunk, int nchunks, TString strdate){
  constexpr bool doOldSelection = true;

  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  if (procsel<0 || procsel>1) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  constexpr int nchannels = 2; // ichannel=0, 1, 2, 3 for ee, mumu, emu, ee+mumu

  TString const coutput_main = "output/LepEffFromData/SkimTrees/" + strdate;
  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", "store:200101");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);

  std::vector<std::string> triggerCheckList = procsel==0 ?
    OffshellTriggerHelpers::getHLTMenus(
      {
        OffshellTriggerHelpers::kDoubleEle, OffshellTriggerHelpers::kSingleEle, OffshellTriggerHelpers::kMuEle
      }
    )
    :
    OffshellTriggerHelpers::getHLTMenus(
      {
        OffshellTriggerHelpers::kDoubleMu, OffshellTriggerHelpers::kSingleMu, OffshellTriggerHelpers::kMuEle
      }
    )
    ;

  std::vector<SampleSpecs> sampleList;
  sampleList.emplace_back("EGamma_2018B", "EGamma 2018B", "EGamma_2018B", -1, HistogramProperties((int) kYellow-3, 1, 2));
  sampleList.emplace_back("DoubleMuon_2018B", "DoubleMuon 2018B", "DoubleMuon_2018B", -1, HistogramProperties((int) kTeal-1, 1, 2));

  // Get handlers
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  ParticleDisambiguator particleDisambiguator;
  DileptonHandler dileptonHandler;

  eventFilter.setTrackDataEvents(false);
  eventFilter.setCheckUniqueDataEvent(false);

  std::vector< std::vector<CutSpecs> > cutsets; getCutSets(cutsets);

  for (size_t isample=0; isample<sampleList.size(); isample++){
    if (procsel>=0 && isample!=static_cast<size_t>(procsel)) continue;

    auto& sample = sampleList.at(isample);

    bool const isData = SampleHelpers::checkSampleIsData(sample.path);
    if (!isData) continue;

    std::vector<TString> sampledirs;
    SampleHelpers::constructSamplesList(sample.path, theGlobalSyst, sampledirs);

    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sampledirs.front()), EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sampledirs.front());

    TString stroutput = Form("%s/%s", coutput_main.Data(), sample.name.data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    TTree* tout = new TTree("SkimTree", "");
#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
    // Event variables
    BRANCH_COMMAND(float, event_pTmiss);
    BRANCH_COMMAND(float, event_phimiss);
    BRANCH_COMMAND(unsigned int, event_Njets);
    // LL
    BRANCH_COMMAND(bool, isOS_ll);
    BRANCH_COMMAND(float, pt_ll);
    BRANCH_COMMAND(float, eta_ll);
    BRANCH_COMMAND(float, phi_ll);
    BRANCH_COMMAND(float, mass_ll);
    BRANCH_COMMAND(bool, pass_cutbasedloose);
    BRANCH_COMMAND(bool, pass_cutbasedmedium);
    BRANCH_COMMAND(bool, pass_cutbasedtight);
    BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp90);
    BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp80);
    BRANCH_COMMAND(bool, pass_fall17v2mvaisowp90);
    BRANCH_COMMAND(bool, pass_fall17v2mvaisowp80);
    BRANCH_COMMAND(bool, pass_hzzmvaisowphzz);
    // L1
    BRANCH_COMMAND(float, pt_l1);
    BRANCH_COMMAND(float, eta_l1);
    BRANCH_COMMAND(float, phi_l1);
    BRANCH_COMMAND(float, mass_l1);
    BRANCH_COMMAND(float, relpfiso03_l1);
    BRANCH_COMMAND(float, relpfiso04_l1);
    // L2
    BRANCH_COMMAND(float, pt_l2);
    BRANCH_COMMAND(float, eta_l2);
    BRANCH_COMMAND(float, phi_l2);
    BRANCH_COMMAND(float, mass_l2);
    BRANCH_COMMAND(float, relpfiso03_l2);
    BRANCH_COMMAND(float, relpfiso04_l2);
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

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    foutput->cd();

    const int nevents = sample_tree.getSelectedNEvents();
    size_t nacc=0;
    int ev_start = 0;
    int ev_end = nevents;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nevents, (ichunk == nchunks-1 ? nevents : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nevents << " events, starting from " << ev_start << " and ending at " << ev_end << "..." << endl;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nevents);
      if (ev % 10000 == 0) MELAout << "Number of accumulated events: " << nacc << '/' << ev << '/' << nevents << endl;
      //if (nacc>1000) break;
      sample_tree.getSelectedEvent(ev);

      float wgt = 1;

      eventFilter.constructFilters();
      if (!eventFilter.passMETFilters() || !eventFilter.passCommonSkim() || !eventFilter.hasGoodVertex()) continue;

      float triggerWeight = eventFilter.getTriggerWeight(triggerCheckList);
      wgt *= triggerWeight;
      if (wgt==0.f) continue;

      muonHandler.constructMuons(theGlobalSyst);
      auto const& muons = muonHandler.getProducts();
      // FIXME: In the future, loose should be even looser if other ids are implemented
      std::vector<MuonObject*> looseMuons;
      for (auto* muon:muons){
        if (ParticleSelectionHelpers::isLooseParticle(muon) || muon->testSelectionBit(MuonSelectionHelpers::kLooseId)) looseMuons.push_back(muon);
      }
      //MELAout << "N loose muons: " << looseMuons.size() << endl;

      electronHandler.constructElectrons(theGlobalSyst);
      auto const& electrons = electronHandler.getProducts();
      std::vector<ElectronObject*> looseElectrons;
      for (auto* electron:electrons){
        auto const& extras = electron->extras;
        bool isCutBasedLoose = ((extras.id_cutBased_Fall17V2_Loose_Bits & 895) == 895);
        bool isFall17V2MVANoIsoLoose = extras.id_MVA_Fall17V2_NoIso_pass_wpLoose;
        bool isFall17V2MVAIsoLoose = extras.id_MVA_Fall17V2_Iso_pass_wpLoose;
        bool isHZZMVAIsoLoose = extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
        if (isCutBasedLoose || isFall17V2MVANoIsoLoose || isFall17V2MVAIsoLoose || isHZZMVAIsoLoose) electron->setSelectionBit(ElectronSelectionHelpers::kLooseId, true);
        else continue;
        bool passLooseIso = (isFall17V2MVAIsoLoose || isHZZMVAIsoLoose);
        if (!passLooseIso){
          if (isCutBasedLoose) passLooseIso |= HelperFunctions::test_bit(extras.id_cutBased_Fall17V2_Loose_Bits, 7);
          else if (isFall17V2MVANoIsoLoose) passLooseIso |= ElectronSelectionHelpers::relPFIso_DR0p3(*electron)<0.15;
        }
        electron->setSelectionBit(ElectronSelectionHelpers::kLooseIso, passLooseIso);
        looseElectrons.push_back(electron);
      }
      //MELAout << "N loose electrons: " << looseElectrons.size() << endl;
      // Disambiguate electrons and muons
      size_t nLooseElectrons_precleaning = looseElectrons.size();
      {
        std::vector<ElectronObject*> looseElectrons_new; looseElectrons_new.reserve(looseElectrons.size());
        for (auto*& product:looseElectrons){
          bool doRemove = false;
          for (auto const* part:looseMuons){
            if (reco::deltaR(product->p4(), part->p4())<0.05){ doRemove=true; break; }
          }
          if (!doRemove) looseElectrons_new.push_back(product);
        }
        looseElectrons = looseElectrons_new;
      }
      size_t nLooseElectrons_postcleaning = looseElectrons.size();
      //MELAout << "N loose electrons after disambiguation: " << looseElectrons.size() << endl;

      photonHandler.constructPhotons(theGlobalSyst);
      auto const& photons = photonHandler.getProducts();
      bool hasLoosePhotons = false;
      for (auto const& photon:photons){
        if (ParticleSelectionHelpers::isLooseParticle(photon)){
          bool doSkip = false;
          for (auto const& part:looseElectrons){
            if (reco::deltaR(photon->p4(), part->p4())<0.1){
              doSkip = true;
              break;
            }
          }
          for (auto const& part:looseMuons){
            if (reco::deltaR(photon->p4(), part->p4())<0.1){
              doSkip = true;
              break;
            }
          }
          if (doSkip) continue;
          hasLoosePhotons = true;
          break;
        }
      }
      if (hasLoosePhotons) continue;
      //MELAout << "Pass photon veto" << endl;

      if (
          !(
            (looseMuons.size()==2 && looseElectrons.size()==0)
            ||
            (looseMuons.size()==0 && looseElectrons.size()==2)
          )
        ) continue;
      //MELAout << "Pass lepton count veto" << endl;

      bool is_ee = looseElectrons.size()==2;
      bool is_mumu = looseMuons.size()==2;
      constexpr bool is_emu=false;

      if (procsel==0 && !is_ee) continue;
      else if (procsel==1 && !is_mumu) continue;

      ParticleObject* leadingLepton = (is_ee ? dynamic_cast<ParticleObject*>(looseElectrons.front()) : dynamic_cast<ParticleObject*>(looseMuons.front()));
      ParticleObject* subleadingLepton = (is_ee ? dynamic_cast<ParticleObject*>(looseElectrons.back()) : dynamic_cast<ParticleObject*>(looseMuons.back()));
      isOS_ll = (leadingLepton->pdgId()*subleadingLepton->pdgId()<0);
      if (leadingLepton->pdgId()*subleadingLepton->pdgId()==0){
        MELAerr << "Leading lepton id (" << leadingLepton->pdgId() << ") or subleading lepton id (" << subleadingLepton->pdgId() << ") are invalid." << endl;
        assert(0);
      }
      //MELAout << "Pass OSSF veto" << endl;
      //if (leadingLepton->pt()<25. || subleadingLepton->pt()<20.) continue;
      //MELAout << "Pass lepton pT veto" << endl;

      ParticleObject::LorentzVector_t p4_ll = leadingLepton->p4() + subleadingLepton->p4();
      float mll = p4_ll.M();
      if (mll<60. || mll>=125.) continue;

      // Now jets
      jetHandler.constructJetMET(theGlobalSyst, &looseMuons, &looseElectrons, nullptr); // Since events with loose photons are skipped, no longer need to pass photons here
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();

      std::vector<AK4JetObject*> ak4jets_tight; ak4jets_tight.reserve(ak4jets.size());
      std::vector<AK4JetObject*> ak4jets_tight_btagged; ak4jets_tight_btagged.reserve(ak4jets.size());
      for (auto* jet:ak4jets){
        if (ParticleSelectionHelpers::isTightJet(jet)){
          ak4jets_tight.push_back(jet);
          if (jet->getBtagValue()>=btagvalue_thr) ak4jets_tight_btagged.push_back(jet);
        }
      }
      if (!ak4jets_tight_btagged.empty()) continue;
      size_t n_ak4jets_tight = ak4jets_tight.size();
      //MELAout << "Pass b veto" << endl;

      if (!eventFilter.test2018HEMFilter(nullptr, &looseElectrons, nullptr, &ak4jets, &ak8jets)) continue; // Since events with loose photons are skipped, no longer need to pass photons here
      //MELAout << "Pass HEM filter" << endl;

      bool pass_ee_cutbasedLoose = false;
      bool pass_ee_cutbasedMedium = false;
      bool pass_ee_cutbasedTight = false;
      bool pass_ee_fall17v2mvanoisowp90 = false;
      bool pass_ee_fall17v2mvanoisowp80 = false;
      bool pass_ee_fall17v2mvaisowp90 = false;
      bool pass_ee_fall17v2mvaisowp80 = false;
      bool pass_ee_hzzmvaisowpHZZ = false;
      if (is_ee){
        pass_ee_cutbasedLoose = true;
        pass_ee_cutbasedMedium = true;
        pass_ee_cutbasedTight = true;
        pass_ee_fall17v2mvanoisowp90 = true;
        pass_ee_fall17v2mvanoisowp80 = true;
        pass_ee_fall17v2mvaisowp90 = true;
        pass_ee_fall17v2mvaisowp80 = true;
        pass_ee_hzzmvaisowpHZZ = true;
        for (auto const& electron:looseElectrons){
          auto const& extras = electron->extras;
          //float relpfiso03 = ElectronSelectionHelpers::relPFIso_DR0p3(*electron);
          pass_ee_cutbasedLoose &= ((extras.id_cutBased_Fall17V2_Loose_Bits & 1023) == 1023);
          pass_ee_cutbasedMedium &= ((extras.id_cutBased_Fall17V2_Medium_Bits & 1023) == 1023);
          pass_ee_cutbasedTight &= ((extras.id_cutBased_Fall17V2_Tight_Bits & 1023) == 1023);
          pass_ee_fall17v2mvanoisowp90 &= extras.id_MVA_Fall17V2_NoIso_pass_wp90;
          pass_ee_fall17v2mvanoisowp80 &= extras.id_MVA_Fall17V2_NoIso_pass_wp80;
          //pass_ee_fall17v2mvanoisowp90 &= extras.id_MVA_Fall17V2_NoIso_pass_wp90 && relpfiso03<ElectronSelectionHelpers::isoThr_medium;
          //pass_ee_fall17v2mvanoisowp80 &= extras.id_MVA_Fall17V2_NoIso_pass_wp80 && relpfiso03<ElectronSelectionHelpers::isoThr_tight;
          pass_ee_fall17v2mvaisowp90 &= extras.id_MVA_Fall17V2_Iso_pass_wp90;
          pass_ee_fall17v2mvaisowp80 &= extras.id_MVA_Fall17V2_Iso_pass_wp80;
          pass_ee_hzzmvaisowpHZZ &= extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
        }
        if (
          //pass_ee_cutbasedLoose ||
          pass_ee_cutbasedMedium ||
          pass_ee_cutbasedTight ||
          pass_ee_fall17v2mvanoisowp90 ||
          pass_ee_fall17v2mvanoisowp80 ||
          pass_ee_fall17v2mvaisowp90 ||
          pass_ee_fall17v2mvaisowp80 ||
          pass_ee_hzzmvaisowpHZZ
        ) nacc++;
      }

      bool pass_mumu_cutbasedLoose = false;
      bool pass_mumu_cutbasedMedium = false;
      bool pass_mumu_cutbasedTight = false;
      if (is_mumu){
        pass_mumu_cutbasedLoose = true;
        pass_mumu_cutbasedMedium = true;
        pass_mumu_cutbasedTight = true;
        for (auto const& muon:looseMuons){
          auto const& extras = muon->extras;
          //float relpfiso03 = MuonSelectionHelpers::relPFIso_DR0p3(*muon);
          pass_mumu_cutbasedLoose &= ((extras.POG_selector_bits & Muon::CutBasedIdLoose) == Muon::CutBasedIdLoose);
          pass_mumu_cutbasedMedium &= ((extras.POG_selector_bits & Muon::CutBasedIdMediumPrompt) == Muon::CutBasedIdMediumPrompt);
          pass_mumu_cutbasedTight &= ((extras.POG_selector_bits & Muon::CutBasedIdTight) == Muon::CutBasedIdTight);
          //pass_mumu_cutbasedMedium &= ((extras.POG_selector_bits & Muon::CutBasedIdMediumPrompt) == Muon::CutBasedIdMediumPrompt) && relpfiso03<MuonSelectionHelpers::isoThr_medium;
          //pass_mumu_cutbasedTight &= ((extras.POG_selector_bits & Muon::CutBasedIdTight) == Muon::CutBasedIdTight) && relpfiso03<MuonSelectionHelpers::isoThr_tight;
        }
        if (
          //pass_mumu_cutbasedLoose ||
          pass_mumu_cutbasedMedium ||
          pass_mumu_cutbasedTight
          ) nacc++;
      }

      // MET
      event_pTmiss = pfmet->pt();
      event_phimiss = pfmet->phi();
      event_Njets = n_ak4jets_tight;
      // LL
      pt_ll = p4_ll.Pt();
      eta_ll = p4_ll.Eta();
      phi_ll = p4_ll.Phi();
      mass_ll = mll;
      pass_cutbasedloose = (is_ee ? pass_ee_cutbasedLoose : pass_mumu_cutbasedLoose);
      pass_cutbasedmedium = (is_ee ? pass_ee_cutbasedMedium : pass_mumu_cutbasedMedium);
      pass_cutbasedtight = (is_ee ? pass_ee_cutbasedTight : pass_mumu_cutbasedTight);
      pass_fall17v2mvanoisowp90 = pass_ee_fall17v2mvanoisowp90;
      pass_fall17v2mvanoisowp80 = pass_ee_fall17v2mvanoisowp80;
      pass_fall17v2mvaisowp90 = pass_ee_fall17v2mvaisowp90;
      pass_fall17v2mvaisowp80 = pass_ee_fall17v2mvaisowp80;
      pass_hzzmvaisowphzz = pass_ee_hzzmvaisowpHZZ;
      // L1
      pt_l1 = leadingLepton->pt();
      eta_l1 = leadingLepton->eta();
      phi_l1 = leadingLepton->phi();
      mass_l1 = leadingLepton->m();
      relpfiso03_l1 = (is_ee ? ElectronSelectionHelpers::relPFIso_DR0p3(*(dynamic_cast<ElectronObject*>(leadingLepton))) : MuonSelectionHelpers::relPFIso_DR0p3(*(dynamic_cast<MuonObject*>(leadingLepton))));
      relpfiso04_l1 = (is_ee ? ElectronSelectionHelpers::relPFIso_DR0p4(*(dynamic_cast<ElectronObject*>(leadingLepton))) : MuonSelectionHelpers::relPFIso_DR0p4(*(dynamic_cast<MuonObject*>(leadingLepton))));
      // L2
      pt_l2 = subleadingLepton->pt();
      eta_l2 = subleadingLepton->eta();
      phi_l2 = subleadingLepton->phi();
      mass_l2 = subleadingLepton->m();
      relpfiso03_l2 = (is_ee ? ElectronSelectionHelpers::relPFIso_DR0p3(*(dynamic_cast<ElectronObject*>(subleadingLepton))) : MuonSelectionHelpers::relPFIso_DR0p3(*(dynamic_cast<MuonObject*>(subleadingLepton))));
      relpfiso04_l2 = (is_ee ? ElectronSelectionHelpers::relPFIso_DR0p4(*(dynamic_cast<ElectronObject*>(subleadingLepton))) : MuonSelectionHelpers::relPFIso_DR0p4(*(dynamic_cast<MuonObject*>(subleadingLepton))));

      tout->Fill();
    } // End loop over events

    foutput->WriteTObject(tout);
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples
}

void makePlots(int do_ee_mumu, bool doOS=true, int doDR0p4=false, bool doRatio=false, TString strdate=""){
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString const cinput_main = "output/LepEffFromData/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Plots" + (doOS? "/OS" : "/SS") + (do_ee_mumu == 0 ? "ee" : "mumu") + (!doDR0p4 ? "/DR0p3" : "/DR0p4") + (doRatio ? "/Ratios" : "/Nevents");
  gSystem->mkdir(coutput_main, true);

  TString const strRelPFIso = (doDR0p4 ? "relPFIso04" : "relPFIso03");
  float const relPFIsoCut = (doDR0p4 ? 0.4 : 0.3);

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

  TString stroutput_txt = Form("%s/Integrals.txt", coutput_main.Data());
  MELAout.open(stroutput_txt.Data());

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(bool, isOS_ll) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(bool, pass_cutbasedloose) \
  BRANCH_COMMAND(bool, pass_cutbasedmedium) \
  BRANCH_COMMAND(bool, pass_cutbasedtight) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp80) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp80) \
  BRANCH_COMMAND(bool, pass_hzzmvaisowphzz) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(float, phi_l1) \
  BRANCH_COMMAND(float, mass_l1) \
  BRANCH_COMMAND(float, relpfiso03_l1) \
  BRANCH_COMMAND(float, relpfiso04_l1) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, phi_l2) \
  BRANCH_COMMAND(float, mass_l2) \
  BRANCH_COMMAND(float, relpfiso03_l2) \
  BRANCH_COMMAND(float, relpfiso04_l2)

#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<TString> sampleList;
  float lumi = SampleHelpers::getIntegratedLuminosity("2018B");
  if (do_ee_mumu==0) sampleList.push_back("EGamma_2018B");
  else sampleList.push_back("DoubleMuon_2018B");
  /*
  sampleList.emplace_back("DY_2l_M_10to50");
  sampleList.emplace_back("DY_2l_M_50_HT");
  sampleList.emplace_back("TT2L2Nu");
  sampleList.emplace_back("ZZ2L2Nu");
  sampleList.emplace_back("WW2L2Nu");
  */
  const size_t nsamples = sampleList.size();
  std::vector<TFile*> finputlist;
  std::vector<TTree*> treelist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample + ".root";
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

  TString strChannel, strChannelLabel;
  getChannelTitleLabel(do_ee_mumu, strChannel, strChannelLabel);
  if (!doOS){
    HelperFunctions::replaceString(strChannel, "OS", "SS");
    HelperFunctions::replaceString(strChannelLabel, "OS", "SS");
  }
  std::vector< std::vector<CutSpecs> > cutsets; getCutSets(cutsets);
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

    const int nbins = ((!doOS && (doRatio || do_ee_mumu==1)) ? 20 : 60);
    constexpr float xinf = 61;
    constexpr float xsup = 121;
    TH1F defaultSelection("h1D_default", "", nbins, xinf, xsup); defaultSelection.Sumw2();
    std::vector<HistogramObject> hspecs; std::vector<int> hcolors;
    if (do_ee_mumu==0){
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_fall17v2mvanoisowp90_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kYellow-3);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_fall17v2mvanoisowp90_%s_0p1", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.1|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kViolet);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "mll_hzzmvawphzz", cuttitle.Data()),
        Form("M_{18} (no iso., HZZ)|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kAzure-2);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedloose", cuttitle.Data()),
        Form("C. B. L.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      /*
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedmedium", cuttitle.Data()),
        Form("C. B. M.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      */
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedtight", cuttitle.Data()),
        Form("C. B. T.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    else if (do_ee_mumu==1){
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedloose_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedloose_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedmedium_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedmedium_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      if (!doDR0p4){
        hspecs.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedtight_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
          Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
          "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
        ); hcolors.push_back((int) kBlack);
      }
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("mll_cutbasedtight_%s_0p15", "relPFIso04"), cuttitle.Data()),
        Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", 0.4, strChannelLabel.Data(), cutlabel.Data()),
        "m_{ll} (GeV)", (!doRatio ? "Events / 1 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    std::vector<TString> hnames; std::vector<TString> htitles;
    for (auto const& hspec:hspecs){
      hnames.push_back(hspec.name);
      htitles.push_back(hspec.title);
    }
    MELAout << "Histogram names: " << hnames << endl;

    bool firstHistogram = true;
    int ibin_sb1=1;
    int jbin_sb1=0, ibin_sr=0, jbin_sr=0, ibin_sb2=0, jbin_sb2=0;
    std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;
    for (unsigned int isample=0; isample<nsamples; isample++){
      std::vector<TH1F*> hlist; hlist.reserve(hspecs.size());
      TFile* const& finput = finputlist.at(isample);
      TTree* const& tin = treelist.at(isample);
      finput->cd();

      for (auto& hspec:hspecs){
        TH1F* htmp = &(hspec.hist); hlist.push_back(htmp);
        htmp->SetLineWidth(2);
        htmp->Sumw2();
      }

      for (int ev=0; ev<tin->GetEntries(); ev++){
        tin->GetEntry(ev);

        if (isOS_ll ^ doOS) continue;

        bool doFill=true;
        for (auto const& cut:cutset){
          TString const& cutvar = cut.cutvar;
          float cutval=0;
          if (cutvar == "Nj") cutval = event_Njets;
          else if (cutvar == "pfmet_pTmiss") cutval = event_pTmiss;
          else if (cutvar == "pTl1") cutval = pt_l1;
          else if (cutvar == "pTl2") cutval = pt_l2;
          doFill &= cut.testCut(cutval);
        }
        if (!doFill) continue;

        const float& relpfiso_l1 = (doDR0p4 ? relpfiso04_l1 : relpfiso03_l1);
        const float& relpfiso_l2 = (doDR0p4 ? relpfiso04_l2 : relpfiso03_l2);

        if (do_ee_mumu==0){
          if (pass_cutbasedtight) defaultSelection.Fill(mass_ll);
        }
        else{
          if (pass_cutbasedtight && std::max(relpfiso04_l1, relpfiso04_l2)<0.15) defaultSelection.Fill(mass_ll);
        }
        auto it_hist = hlist.begin();
        if (do_ee_mumu==0){
          if (pass_fall17v2mvanoisowp90 && relpfiso_l1<0.2 && relpfiso_l2<0.2) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_fall17v2mvanoisowp90 && relpfiso_l1<0.1 && relpfiso_l2<0.1) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_hzzmvaisowphzz) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_cutbasedloose) (*it_hist)->Fill(mass_ll); it_hist++;
          /*
          if (pass_cutbasedmedium) (*it_hist)->Fill(mass_ll); it_hist++;
          */
          if (pass_cutbasedtight) (*it_hist)->Fill(mass_ll); it_hist++;
        }
        else{
          if (pass_cutbasedloose && relpfiso_l1<0.2 && relpfiso_l2<0.2) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_cutbasedloose && relpfiso_l1<0.15 && relpfiso_l2<0.15) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_cutbasedmedium && relpfiso_l1<0.2 && relpfiso_l2<0.2) (*it_hist)->Fill(mass_ll); it_hist++;
          if (pass_cutbasedmedium && relpfiso_l1<0.15 && relpfiso_l2<0.15) (*it_hist)->Fill(mass_ll); it_hist++;
          if (!doDR0p4){
            if (pass_cutbasedtight && relpfiso_l1<0.15 && relpfiso_l2<0.15) (*it_hist)->Fill(mass_ll); it_hist++;
          }
          if (pass_cutbasedtight && relpfiso04_l1<0.15 && relpfiso04_l2<0.15) (*it_hist)->Fill(mass_ll); it_hist++;
        }
      }

      unsigned int ihist=0;
      for (auto const& hspec:hspecs){
        TH1F* htmp = hlist.at(ihist);
        if (firstHistogram){
          firstHistogram=false;
          TAxis const* xaxis = htmp->GetXaxis();
          jbin_sb1 = xaxis->FindBin(75);
          ibin_sr = xaxis->FindBin(76.2);
          jbin_sr = xaxis->FindBin(105.8);
          ibin_sb2 = xaxis->FindBin(106.2);
          jbin_sb2 = nbins;
        }
        htmp->SetFillStyle(0);
        htmp->SetLineColor(hcolors.at(ihist));
        htmp->SetMarkerColor(hcolors.at(ihist));
        if (do_ee_mumu==0){
        }
        else{
          if (ihist%2==0) htmp->SetLineStyle(7);
        }
        ihist++;
      }
      sample_hist_map[sampleList.at(isample)]=hlist;
    }


    double ymin = 0;
    double ymax = -1;
    for (auto& it:sample_hist_map){
      TH1F* hdivide=nullptr;
      if (doRatio){
        hdivide = &defaultSelection;
      }
      for (auto& hist:it.second){
        if (hdivide) HelperFunctions::divideHistograms(hist, hdivide, hist, true);
        for (int ix=1; ix<=nbins; ix++){
          double bc = hist->GetBinContent(ix);
          double be = hist->GetBinError(ix);
          if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
          ymax = std::max(ymax, bc+be);
          if (bc>0.){
            if (ymin<=0.) ymin = bc;
            else ymin = std::min(ymin, bc);
          }
        }
        // Now print integrals
        if (!doRatio){
          double int_sb1 = hist->Integral(ibin_sb1, jbin_sb1);
          double int_sb2 = hist->Integral(ibin_sb2, jbin_sb2);
          double int_sr = hist->Integral(ibin_sr, jbin_sr);
          MELAout << "Integrals[" << it.first << "::" << hist->GetName() << "](sb1, sb2, sb1+sb2, sr, sr/sb1+sb2) = ( " << int_sb1 << ", " << int_sb2 << ", " << int_sb1+int_sb2 << ", " << int_sr << ", " << int_sr/(int_sb1+int_sb2) << " )" << endl;
        }
      }
    }
    double yinf=ymin;
    double ysup=ymax;
    //ymax *= (!doRatio ? 15. : 1.2);
    ymax *= (!doRatio ? 1.7 : 1.5);
    if (doRatio){
      ymin *= 0.95;
      if (!doOS) ymax=std::min(ymax, 5.);
    }
    else ymin=0;
    for (auto const& it:sample_hist_map){
      for (auto const& hist:it.second) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }


    TString canvasname = Form("%s_%s_%s", strChannel.Data(), "mll", cuttitle.Data());
    if (doRatio) canvasname = Form("ratios_%s", canvasname.Data());
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
    //if (!doRatio) canvas->SetLogy();

    TLegend* legend = new TLegend(
      0.55,
      0.90-0.10/4.*2.*float(hnames.size()),
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
    text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (auto& it:sample_hist_map){
      for (auto& hist:it.second){
        //if (doRatio && hist==it.second.back()) continue;
        std::vector<TString> tmplist;
        TString htitle = hist->GetTitle();
        HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
        TString hlabel = tmplist.front();

        hist->SetTitle("");
        legend->AddEntry(hist, hlabel, "l");

        if (firstHist){
          for (size_t is=1; is<tmplist.size(); is++) selectionList.push_back(tmplist.at(is));
          hist->Draw("hist");
          firstHist = false;
        }
        else{
          hist->Draw("histsame");
        }
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
    if (!doRatio){
      for (auto& it:sample_hist_map){
        for (auto const& hist:it.second) hist->GetYaxis()->SetRangeUser(yinf*0.5, ysup*15.);
      }
      canvas->Modified();
      canvas->Update();
      canvas->SetLogy();
      canvas->Modified();
      canvas->Update();
      canvas->SaveAs(coutput_main + "/logY_" + canvasname + ".pdf");
      canvas->SaveAs(coutput_main + "/logY_" + canvasname + ".png");
    }

    for (auto*& ptSel:ptSelectionList) delete ptSel;
    delete pt;
    delete legend;
    canvas->Close();
  }

  for (TFile*& finput:finputlist) finput->Close();
  MELAout.close();
}

void makePlots_pTll(int do_ee_mumu, bool doOS=true, bool doSB=false, bool doRatio=false, TString strdate=""){
  if (do_ee_mumu==1 && !doOS && doSB) return;
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString const cinput_main = "output/LepEffFromData/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Plots" + (doOS ? "/OS" : "/SS") + (do_ee_mumu==0 ? "ee" : "mumu") + "/pTllComparisons" + (!doSB ? (do_ee_mumu==1 && !doOS ? "/SRandSideband" : "/SignalRegion") : "/Sidebands") + (doRatio ? "/Ratios" : "/Nevents");
  gSystem->mkdir(coutput_main, true);

  TString strRelPFIso; float relPFIsoCut;

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

  TString stroutput_txt = Form("%s/Integrals.txt", coutput_main.Data());
  MELAout.open(stroutput_txt.Data());

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(bool, isOS_ll) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(bool, pass_cutbasedloose) \
  BRANCH_COMMAND(bool, pass_cutbasedmedium) \
  BRANCH_COMMAND(bool, pass_cutbasedtight) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp80) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp80) \
  BRANCH_COMMAND(bool, pass_hzzmvaisowphzz) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(float, phi_l1) \
  BRANCH_COMMAND(float, mass_l1) \
  BRANCH_COMMAND(float, relpfiso03_l1) \
  BRANCH_COMMAND(float, relpfiso04_l1) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, phi_l2) \
  BRANCH_COMMAND(float, mass_l2) \
  BRANCH_COMMAND(float, relpfiso03_l2) \
  BRANCH_COMMAND(float, relpfiso04_l2)

#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<TString> sampleList;
  float lumi = SampleHelpers::getIntegratedLuminosity("2018B");
  if (do_ee_mumu==0) sampleList.push_back("EGamma_2018B");
  else sampleList.push_back("DoubleMuon_2018B");

  const size_t nsamples = sampleList.size();
  std::vector<TFile*> finputlist;
  std::vector<TTree*> treelist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample + ".root";
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

  TString strChannel, strChannelLabel;
  getChannelTitleLabel(do_ee_mumu, strChannel, strChannelLabel);
  if (!doOS){
    HelperFunctions::replaceString(strChannel, "OS", "SS");
    HelperFunctions::replaceString(strChannelLabel, "OS", "SS");
  }
  std::vector< std::vector<CutSpecs> > cutsets; getCutSets(cutsets);
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
    if (do_ee_mumu==1 && !doOS && !doSB){}
    else if (doSB){
      cutlabel = cutlabel + '|' + "Sidebands";
      cuttitle = cuttitle + '_' + "sidebands";
    }
    else{
      cutlabel = cutlabel + '|' + "Z peak";
      cuttitle = cuttitle + '_' + "Zpeak";
    }

    const int nbins = (doRatio ? 20 : 50);
    constexpr float xinf = 0;
    constexpr float xsup = 150;
    TH1F defaultSelection("h1D_default", "", nbins, xinf, xsup); defaultSelection.Sumw2();
    std::vector<HistogramObject> hspecs; std::vector<int> hcolors;
    if (do_ee_mumu==0){
      strRelPFIso = "relPFIso03";
      relPFIsoCut = 0.3;

      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_fall17v2mvanoisowp90_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kYellow-3);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_fall17v2mvanoisowp90_%s_0p1", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.1|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kViolet);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTll_hzzmvawphzz", cuttitle.Data()),
        Form("M_{18} (no iso., HZZ)|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kAzure-2);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTll_cutbasedloose", cuttitle.Data()),
        Form("C. B. L.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTll_cutbasedtight", cuttitle.Data()),
        Form("C. B. T.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    else if (do_ee_mumu==1){
      strRelPFIso = "relPFIso03";
      relPFIsoCut = 0.3;

      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedloose_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedloose_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedmedium_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedmedium_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedtight_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);

      strRelPFIso = "relPFIso04";
      relPFIsoCut = 0.4;
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTll_cutbasedtight_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{ll} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    std::vector<TString> hnames; std::vector<TString> htitles;
    for (auto const& hspec:hspecs){
      hnames.push_back(hspec.name);
      htitles.push_back(hspec.title);
    }
    MELAout << "Histogram names: " << hnames << endl;

    bool firstHistogram = true;
    int ibin_sb1=1;
    int jbin_sb1=0, ibin_sr=0, jbin_sr=0, ibin_sb2=0, jbin_sb2=0;
    std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;
    for (unsigned int isample=0; isample<nsamples; isample++){
      std::vector<TH1F*> hlist; hlist.reserve(hspecs.size());
      TFile* const& finput = finputlist.at(isample);
      TTree* const& tin = treelist.at(isample);
      finput->cd();

      for (auto& hspec:hspecs){
        TH1F* htmp = &(hspec.hist); hlist.push_back(htmp);
        htmp->SetLineWidth(2);
        htmp->Sumw2();
      }

      for (int ev=0; ev<tin->GetEntries(); ev++){
        tin->GetEntry(ev);

        if (isOS_ll ^ doOS) continue;
        if (do_ee_mumu==1 && !doOS && !doSB){}
        else if (!doSB && std::abs(mass_ll-91.2)>=15.) continue;
        else if (doSB && std::abs(mass_ll-91.2)<15.) continue;

        bool doFill=true;
        for (auto const& cut:cutset){
          TString const& cutvar = cut.cutvar;
          float cutval=0;
          if (cutvar == "Nj") cutval = event_Njets;
          else if (cutvar == "pfmet_pTmiss") cutval = event_pTmiss;
          else if (cutvar == "pTl1") cutval = pt_l1;
          else if (cutvar == "pTl2") cutval = pt_l2;
          doFill &= cut.testCut(cutval);
        }
        if (!doFill) continue;

        float max_relpfiso04 = std::max(relpfiso04_l1, relpfiso04_l2);
        float max_relpfiso03 = std::max(relpfiso03_l1, relpfiso03_l2);

        if (do_ee_mumu==0){
          if (pass_cutbasedtight) defaultSelection.Fill(pt_ll);
        }
        else{
          if (pass_cutbasedtight && max_relpfiso04<0.15) defaultSelection.Fill(pt_ll);
        }

        auto it_hist = hlist.begin();
        if (do_ee_mumu==0){
          if (pass_fall17v2mvanoisowp90 && max_relpfiso03<0.2) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_fall17v2mvanoisowp90 && max_relpfiso03<0.1) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_hzzmvaisowphzz) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedloose) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedtight) (*it_hist)->Fill(pt_ll); it_hist++;
        }
        else{
          if (pass_cutbasedloose && max_relpfiso03<0.2) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedloose && max_relpfiso03<0.15) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedmedium && max_relpfiso03<0.2) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedmedium && max_relpfiso03<0.15) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedtight && max_relpfiso03<0.15) (*it_hist)->Fill(pt_ll); it_hist++;
          if (pass_cutbasedtight && max_relpfiso04<0.15) (*it_hist)->Fill(pt_ll); it_hist++;
        }
      }

      unsigned int ihist=0;
      for (auto const& hspec:hspecs){
        TH1F* htmp = hlist.at(ihist);
        if (firstHistogram){
          firstHistogram=false;
          TAxis const* xaxis = htmp->GetXaxis();
          jbin_sb1 = xaxis->FindBin(75);
          ibin_sr = xaxis->FindBin(76.2);
          jbin_sr = xaxis->FindBin(105.8);
          ibin_sb2 = xaxis->FindBin(106.2);
          jbin_sb2 = nbins;
        }
        htmp->SetFillStyle(0);
        htmp->SetLineColor(hcolors.at(ihist));
        htmp->SetMarkerColor(hcolors.at(ihist));
        if (do_ee_mumu==0){
          //if (ihist%3==1) htmp->SetLineStyle(7);
          //else if (ihist%3==2) htmp->SetLineStyle(2);
        }
        else{
          if (ihist%2==0) htmp->SetLineStyle(7);
        }
        ihist++;
      }
      sample_hist_map[sampleList.at(isample)]=hlist;
    }


    double ymin = 0;
    double ymax = -1;
    for (auto& it:sample_hist_map){
      TH1F* hdivide=nullptr;
      if (doRatio){
        hdivide = &defaultSelection;
      }
      for (auto& hist:it.second){
        if (hdivide) HelperFunctions::divideHistograms(hist, hdivide, hist, true);
        for (int ix=1; ix<=nbins; ix++){
          double bc = hist->GetBinContent(ix);
          double be = hist->GetBinError(ix);
          if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
          ymax = std::max(ymax, bc+be);
          if (bc>0.){
            if (ymin<=0.) ymin = bc;
            else ymin = std::min(ymin, bc);
          }
        }
        // Now print integrals
        if (!doRatio){
          double int_sb1 = hist->Integral(ibin_sb1, jbin_sb1);
          double int_sb2 = hist->Integral(ibin_sb2, jbin_sb2);
          double int_sr = hist->Integral(ibin_sr, jbin_sr);
          MELAout << "Integrals[" << it.first << "::" << hist->GetName() << "](sb1, sb2, sb1+sb2, sr, sr/sb1+sb2) = ( " << int_sb1 << ", " << int_sb2 << ", " << int_sb1+int_sb2 << ", " << int_sr << ", " << int_sr/(int_sb1+int_sb2) << " )" << endl;
        }
      }
    }
    //ymax *= (!doRatio ? 15. : 1.5);
    ymax *= (!doRatio ? 1.7 : 1.5);
    if (doRatio){
      ymin *= 0.95;
      if (!doOS) ymax=std::min(ymax, 5.);
    }
    else ymin=0;
    for (auto const& it:sample_hist_map){
      for (auto const& hist:it.second) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }


    TString canvasname = Form("%s_%s_%s", strChannel.Data(), "pTll", cuttitle.Data());
    if (doRatio) canvasname = Form("ratios_%s", canvasname.Data());
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
    //if (!doRatio) canvas->SetLogy();

    TLegend* legend = new TLegend(
      0.55,
      0.90-0.10/4.*2.*float(hnames.size()),
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
    text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (auto& it:sample_hist_map){
      for (auto& hist:it.second){
        std::vector<TString> tmplist;
        TString htitle = hist->GetTitle();
        HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
        TString hlabel = tmplist.front();

        hist->SetTitle("");
        legend->AddEntry(hist, hlabel, "l");

        if (firstHist){
          for (size_t is=1; is<tmplist.size(); is++) selectionList.push_back(tmplist.at(is));
          hist->Draw("hist");
          firstHist = false;
        }
        else{
          hist->Draw("histsame");
        }
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

void makePlots_pTl(int do_ee_mumu, bool doOS=true, bool doSB=false, bool doRatio=false, TString strdate=""){
  if (do_ee_mumu==1 && !doOS && doSB) return;
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TString const cinput_main = "output/LepEffFromData/SkimTrees/" + strdate;
  TString coutput_main = cinput_main + "/Plots" + (doOS ? "/OS" : "/SS") + (do_ee_mumu==0 ? "ee" : "mumu") + "/pTlComparisons" + (!doSB ? (do_ee_mumu==1 && !doOS ? "/SRandSideband" : "/SignalRegion") : "/Sidebands") + (doRatio ? "/Ratios" : "/Nevents");
  gSystem->mkdir(coutput_main, true);

  TString strRelPFIso; float relPFIsoCut;

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

  TString stroutput_txt = Form("%s/Integrals.txt", coutput_main.Data());
  MELAout.open(stroutput_txt.Data());

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(bool, isOS_ll) \
  BRANCH_COMMAND(float, pt_ll) \
  BRANCH_COMMAND(float, eta_ll) \
  BRANCH_COMMAND(float, phi_ll) \
  BRANCH_COMMAND(float, mass_ll) \
  BRANCH_COMMAND(bool, pass_cutbasedloose) \
  BRANCH_COMMAND(bool, pass_cutbasedmedium) \
  BRANCH_COMMAND(bool, pass_cutbasedtight) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvanoisowp80) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp90) \
  BRANCH_COMMAND(bool, pass_fall17v2mvaisowp80) \
  BRANCH_COMMAND(bool, pass_hzzmvaisowphzz) \
  BRANCH_COMMAND(float, pt_l1) \
  BRANCH_COMMAND(float, eta_l1) \
  BRANCH_COMMAND(float, phi_l1) \
  BRANCH_COMMAND(float, mass_l1) \
  BRANCH_COMMAND(float, relpfiso03_l1) \
  BRANCH_COMMAND(float, relpfiso04_l1) \
  BRANCH_COMMAND(float, pt_l2) \
  BRANCH_COMMAND(float, eta_l2) \
  BRANCH_COMMAND(float, phi_l2) \
  BRANCH_COMMAND(float, mass_l2) \
  BRANCH_COMMAND(float, relpfiso03_l2) \
  BRANCH_COMMAND(float, relpfiso04_l2)

#define BRANCH_COMMAND(type, name) type name=0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND

  std::vector<TString> sampleList;
  float lumi = SampleHelpers::getIntegratedLuminosity("2018B");
  if (do_ee_mumu==0) sampleList.push_back("EGamma_2018B");
  else sampleList.push_back("DoubleMuon_2018B");

  const size_t nsamples = sampleList.size();
  std::vector<TFile*> finputlist;
  std::vector<TTree*> treelist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample + ".root";
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

  TString strChannel, strChannelLabel;
  getChannelTitleLabel(do_ee_mumu, strChannel, strChannelLabel);
  if (!doOS){
    HelperFunctions::replaceString(strChannel, "OS", "SS");
    HelperFunctions::replaceString(strChannelLabel, "OS", "SS");
  }
  std::vector< std::vector<CutSpecs> > cutsets; getCutSets(cutsets);
  for (auto const& cutset:cutsets){
    TString cutlabel, cuttitle;
    for (auto it_cut = cutset.cbegin(); it_cut != cutset.cend(); it_cut++){
      TString const& cutvar = it_cut->cutvar;
      if (cutvar == "pTl1" || cutvar == "pTl2") continue;
      if (it_cut == cutset.cbegin()){
        cutlabel = it_cut->getLabel();
        cuttitle = it_cut->getTitle();
      }
      else{
        cutlabel = cutlabel + '|' + it_cut->getLabel();
        cuttitle = cuttitle + '_' + it_cut->getTitle();
      }
    }
    if (do_ee_mumu==1 && !doOS && !doSB){}
    else if (doSB){
      cutlabel = cutlabel + '|' + "Sidebands";
      cuttitle = cuttitle + '_' + "sidebands";
    }
    else{
      cutlabel = cutlabel + '|' + "Z peak";
      cuttitle = cuttitle + '_' + "Zpeak";
    }

    const int nbins = (doRatio ? 20 : 50);
    constexpr float xinf = 20;
    constexpr float xsup = 170;
    TH1F defaultSelection("h1D_default", "", nbins, xinf, xsup); defaultSelection.Sumw2();
    std::vector<HistogramObject> hspecs; std::vector<int> hcolors;
    if (do_ee_mumu==0){
      strRelPFIso = "relPFIso03";
      relPFIsoCut = 0.3;

      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_fall17v2mvanoisowp90_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kYellow-3);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_fall17v2mvanoisowp90_%s_0p1", strRelPFIso.Data()), cuttitle.Data()),
        Form("M_{17}^{V2} (no iso., 90), I_{rel}^{%.1f}<0.1|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kViolet);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTl_hzzmvawphzz", cuttitle.Data()),
        Form("M_{18} (no iso., HZZ)|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Events / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kAzure-2);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTl_cutbasedloose", cuttitle.Data()),
        Form("C. B. L.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), "pTl_cutbasedtight", cuttitle.Data()),
        Form("C. B. T.|%s|%s", strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    else if (do_ee_mumu==1){
      strRelPFIso = "relPFIso03";
      relPFIsoCut = 0.3;

      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedloose_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedloose_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. L., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kRed);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedmedium_%s_0p2", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.2|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedmedium_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. M., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlue);
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedtight_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);

      strRelPFIso = "relPFIso04";
      relPFIsoCut = 0.4;
      hspecs.emplace_back(
        Form("h1D_%s_%s_%s", strChannel.Data(), Form("pTl_cutbasedtight_%s_0p15", strRelPFIso.Data()), cuttitle.Data()),
        Form("C. B. T., I_{rel}^{%.1f}<0.15|%s|%s", relPFIsoCut, strChannelLabel.Data(), cutlabel.Data()),
        "p_{T}^{l1,l2} (GeV)", (!doRatio ? "Leptons / 3 GeV" : "Ratio to default id+iso"), nbins, xinf, xsup, nullptr
      ); hcolors.push_back((int) kBlack);
    }
    std::vector<TString> hnames; std::vector<TString> htitles;
    for (auto const& hspec:hspecs){
      hnames.push_back(hspec.name);
      htitles.push_back(hspec.title);
    }
    MELAout << "Histogram names: " << hnames << endl;

    bool firstHistogram = true;
    int ibin_sb1=1;
    int jbin_sb1=0, ibin_sr=0, jbin_sr=0, ibin_sb2=0, jbin_sb2=0;
    std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;
    for (unsigned int isample=0; isample<nsamples; isample++){
      std::vector<TH1F*> hlist; hlist.reserve(hspecs.size());
      TFile* const& finput = finputlist.at(isample);
      TTree* const& tin = treelist.at(isample);
      finput->cd();

      for (auto& hspec:hspecs){
        TH1F* htmp = &(hspec.hist); hlist.push_back(htmp);
        htmp->SetLineWidth(2);
        htmp->Sumw2();
      }

      for (int ev=0; ev<tin->GetEntries(); ev++){
        tin->GetEntry(ev);

        if (isOS_ll ^ doOS) continue;
        if (do_ee_mumu==1 && !doOS && !doSB){}
        else if (!doSB && std::abs(mass_ll-91.2)>=15.) continue;
        else if (doSB && std::abs(mass_ll-91.2)<15.) continue;

        bool doFill=true;
        for (auto const& cut:cutset){
          TString const& cutvar = cut.cutvar;
          float cutval=0;
          if (cutvar == "Nj") cutval = event_Njets;
          else if (cutvar == "pfmet_pTmiss") cutval = event_pTmiss;
          else if (cutvar == "pTl1" || cutvar == "pTl2") continue; // Skip these cuts
          doFill &= cut.testCut(cutval);
        }
        if (!doFill) continue;

        float max_relpfiso04 = std::max(relpfiso04_l1, relpfiso04_l2);
        float max_relpfiso03 = std::max(relpfiso03_l1, relpfiso03_l2);

        if (do_ee_mumu==0){
          if (pass_cutbasedtight){ defaultSelection.Fill(pt_l1); defaultSelection.Fill(pt_l2); }
        }
        else{
          if (pass_cutbasedtight && max_relpfiso04<0.15){ defaultSelection.Fill(pt_l1); defaultSelection.Fill(pt_l2); }
        }

        auto it_hist = hlist.begin();
        if (do_ee_mumu==0){
          if (pass_fall17v2mvanoisowp90 && max_relpfiso03<0.2){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_fall17v2mvanoisowp90 && max_relpfiso03<0.1){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_hzzmvaisowphzz){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedloose){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedtight){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
        }
        else{
          if (pass_cutbasedloose && max_relpfiso03<0.2){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedloose && max_relpfiso03<0.15){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedmedium && max_relpfiso03<0.2){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedmedium && max_relpfiso03<0.15){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedtight && max_relpfiso03<0.15){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
          if (pass_cutbasedtight && max_relpfiso04<0.15){ (*it_hist)->Fill(pt_l1); (*it_hist)->Fill(pt_l2); }; it_hist++;
        }
      }

      unsigned int ihist=0;
      for (auto const& hspec:hspecs){
        TH1F* htmp = hlist.at(ihist);
        if (firstHistogram){
          firstHistogram=false;
          TAxis const* xaxis = htmp->GetXaxis();
          jbin_sb1 = xaxis->FindBin(75);
          ibin_sr = xaxis->FindBin(76.2);
          jbin_sr = xaxis->FindBin(105.8);
          ibin_sb2 = xaxis->FindBin(106.2);
          jbin_sb2 = nbins;
        }
        htmp->SetFillStyle(0);
        htmp->SetLineColor(hcolors.at(ihist));
        htmp->SetMarkerColor(hcolors.at(ihist));
        if (do_ee_mumu==0){
          //if (ihist%3==1) htmp->SetLineStyle(7);
          //else if (ihist%3==2) htmp->SetLineStyle(2);
        }
        else{
          if (ihist%2==0) htmp->SetLineStyle(7);
        }
        ihist++;
      }
      sample_hist_map[sampleList.at(isample)]=hlist;
    }


    double ymin = 0;
    double ymax = -1;
    for (auto& it:sample_hist_map){
      TH1F* hdivide=nullptr;
      if (doRatio){
        hdivide = &defaultSelection;
      }
      for (auto& hist:it.second){
        if (hdivide) HelperFunctions::divideHistograms(hist, hdivide, hist, true);
        for (int ix=1; ix<=nbins; ix++){
          double bc = hist->GetBinContent(ix);
          double be = hist->GetBinError(ix);
          if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
          ymax = std::max(ymax, bc+be);
          if (bc>0.){
            if (ymin<=0.) ymin = bc;
            else ymin = std::min(ymin, bc);
          }
        }
        // Now print integrals
        if (!doRatio){
          double int_sb1 = hist->Integral(ibin_sb1, jbin_sb1);
          double int_sb2 = hist->Integral(ibin_sb2, jbin_sb2);
          double int_sr = hist->Integral(ibin_sr, jbin_sr);
          MELAout << "Integrals[" << it.first << "::" << hist->GetName() << "](sb1, sb2, sb1+sb2, sr, sr/sb1+sb2) = ( " << int_sb1 << ", " << int_sb2 << ", " << int_sb1+int_sb2 << ", " << int_sr << ", " << int_sr/(int_sb1+int_sb2) << " )" << endl;
        }
      }
    }
    //ymax *= (!doRatio ? 15. : 1.5);
    ymax *= (!doRatio ? 1.7 : 1.5);
    if (doRatio){
      ymin *= 0.95;
      if (!doOS) ymax=std::min(ymax, 5.);
    }
    else ymin=0;
    for (auto const& it:sample_hist_map){
      for (auto const& hist:it.second) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }


    TString canvasname = Form("%s_%s_%s", strChannel.Data(), "pTl", cuttitle.Data());
    if (doRatio) canvasname = Form("ratios_%s", canvasname.Data());
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
    //if (!doRatio) canvas->SetLogy();

    TLegend* legend = new TLegend(
      0.55,
      0.90-0.10/4.*2.*float(hnames.size()),
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
    text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (auto& it:sample_hist_map){
      for (auto& hist:it.second){
        std::vector<TString> tmplist;
        TString htitle = hist->GetTitle();
        HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
        TString hlabel = tmplist.front();

        hist->SetTitle("");
        legend->AddEntry(hist, hlabel, "l");

        if (firstHist){
          for (size_t is=1; is<tmplist.size(); is++) selectionList.push_back(tmplist.at(is));
          hist->Draw("hist");
          firstHist = false;
        }
        else{
          hist->Draw("histsame");
        }
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
