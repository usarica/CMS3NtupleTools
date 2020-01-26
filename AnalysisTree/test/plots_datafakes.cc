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

  MELAout << "Created histogram " << hist.GetName() << " [" << hist.GetTitle() << "]" << endl;
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

void getHistograms(int procsel, int ichunk, int nchunks, TString strdate){
  constexpr bool doOldSelection = true;

  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  if (procsel<0 || procsel>1) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  constexpr int nchannels = 2; // ichannel=0, 1, 2, 3 for ee, mumu, emu, ee+mumu

  TString const coutput_main = "output/LepEffFromData/" + strdate;

  gSystem->mkdir(coutput_main, true);

  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  SampleHelpers::configure("2018", "hadoop:200101");

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btagvalue_thr = BtagHelpers::getBtagWP(false);
  float lumi = SampleHelpers::getIntegratedLuminosity("2018");

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

  std::vector< std::vector<CutSpecs> > cutsets;
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
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 0, 0
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    // Nj==1
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 1, 1
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    // Nj==2
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, true, 2, 2
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
    // Nj>=3
    cutsets.push_back(std::vector<CutSpecs>());
    cutsets.back().reserve(2);
    cutsets.back().emplace_back(
      "Nj", "N_{j}",
      true, false, 3, -1
    );
    cutsets.back().emplace_back(
      "pfmet_pTmiss", "p_{T}^{miss,PF}",
      do_met_low, do_met_high, met_low, met_high
    );
  }

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

    for (int ichannel=0; ichannel<nchannels; ichannel++){
      if (procsel!=ichannel) continue;

      TString strChannel, strChannelLabel;
      getChannelTitleLabel(ichannel, strChannel, strChannelLabel);

      TDirectory* channeldir = foutput->mkdir(strChannel);
      channeldir->cd();

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

        sample.hlist_1D.emplace_back(
          Form("h1D_%s_%s_%s", strChannel.Data(), "mll_anyloose", cuttitle.Data()), Form("Any loose|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
          "m_{ll} (GeV)", "",
          32, 31., 191.,
          channeldir
        );

        if (ichannel==0){
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedmedium", cuttitle.Data()), Form("Cut-based medium|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedtight", cuttitle.Data()), Form("Cut-based tight|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_fall17v2mvanoisowp90", cuttitle.Data()), Form("Fall17V2 MVA (no iso.) WP90|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_fall17v2mvanoisowp80", cuttitle.Data()), Form("Fall17V2 MVA (no iso.) WP80|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_fall17v2mvaisowp90", cuttitle.Data()), Form("Fall17V2 MVA (w/ iso.) WP90|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_fall17v2mvaisowp80", cuttitle.Data()), Form("Fall17V2 MVA (w/ iso.) WP80|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_hzzmvawphzz", cuttitle.Data()), Form("HZZ MVA WPHZZ|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
        }
        else{
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedmedium", cuttitle.Data()), Form("Cut-based medium|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
          sample.hlist_1D.emplace_back(
            Form("h1D_%s_%s_%s", strChannel.Data(), "mll_cutbasedtight", cuttitle.Data()), Form("Cut-based tight|%s (%s)|%s", strChannelLabel.Data(), sample.label.data(), cutlabel.Data()),
            "m_{ll} (GeV)", "",
            32, 31., 191.,
            channeldir
          );
        }

      }

      foutput->cd();
    }

    // Configure histograms
    sample.setup();

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
        if (ParticleSelectionHelpers::isLooseParticle(muon)) looseMuons.push_back(muon);
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
        if (isCutBasedLoose || isFall17V2MVANoIsoLoose || isFall17V2MVAIsoLoose || isHZZMVAIsoLoose) electron->setSelectionBit(ElectronSelectionHelpers::kLooseId);
        else continue;
        bool passLooseIso = (isFall17V2MVAIsoLoose || isHZZMVAIsoLoose);
        if (!passLooseIso){
          if (isCutBasedLoose) passLooseIso |= HelperFunctions::test_bit(extras.id_cutBased_Fall17V2_Loose_Bits, 7);
          else if (isFall17V2MVANoIsoLoose) passLooseIso |= ElectronSelectionHelpers::relPFIso_DR0p3(*electron)<ElectronSelectionHelpers::isoThr_loose;
        }
        if (passLooseIso){
          electron->setSelectionBit(ElectronSelectionHelpers::kLooseIso);
          looseElectrons.push_back(electron);
        }
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
      if (leadingLepton->pdgId()*subleadingLepton->pdgId()>0) continue;
      else if (leadingLepton->pdgId()*subleadingLepton->pdgId()==0){
        MELAerr << "Leading lepton id (" << leadingLepton->pdgId() << ") or subleading lepton id (" << subleadingLepton->pdgId() << ") are invalid." << endl;
        assert(0);
      }
      //MELAout << "Pass OSSF veto" << endl;
      if (leadingLepton->pt()<25. || subleadingLepton->pt()<20.) continue;
      //MELAout << "Pass lepton pT veto" << endl;

      ParticleObject::LorentzVector_t p4_ll = leadingLepton->p4() + subleadingLepton->p4();
      float mll = p4_ll.M();

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

      bool pass_ee_cutbasedMedium = false;
      bool pass_ee_cutbasedTight = false;
      bool pass_ee_fall17v2mvanoisowp90 = false;
      bool pass_ee_fall17v2mvanoisowp80 = false;
      bool pass_ee_fall17v2mvaisowp90 = false;
      bool pass_ee_fall17v2mvaisowp80 = false;
      bool pass_ee_hzzmvaisowpHZZ = false;
      if (is_ee){
        pass_ee_cutbasedMedium = true;
        pass_ee_cutbasedTight = true;
        pass_ee_fall17v2mvanoisowp90 = true;
        pass_ee_fall17v2mvanoisowp80 = true;
        pass_ee_fall17v2mvaisowp90 = true;
        pass_ee_fall17v2mvaisowp80 = true;
        pass_ee_hzzmvaisowpHZZ = true;
        for (auto const& electron:looseElectrons){
          auto const& extras = electron->extras;
          float relpfiso03 = ElectronSelectionHelpers::relPFIso_DR0p3(*electron);
          pass_ee_cutbasedMedium &= ((extras.id_cutBased_Fall17V2_Medium_Bits & 1023) == 1023);
          pass_ee_cutbasedTight &= ((extras.id_cutBased_Fall17V2_Tight_Bits & 1023) == 1023);
          pass_ee_fall17v2mvanoisowp90 &= extras.id_MVA_Fall17V2_NoIso_pass_wp90 && relpfiso03<ElectronSelectionHelpers::isoThr_medium;
          pass_ee_fall17v2mvanoisowp80 &= extras.id_MVA_Fall17V2_NoIso_pass_wp80 && relpfiso03<ElectronSelectionHelpers::isoThr_tight;
          pass_ee_fall17v2mvaisowp90 &= extras.id_MVA_Fall17V2_Iso_pass_wp90;
          pass_ee_fall17v2mvaisowp80 &= extras.id_MVA_Fall17V2_Iso_pass_wp80;
          pass_ee_hzzmvaisowpHZZ &= extras.id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ;
        }
        if (
          pass_ee_cutbasedMedium ||
          pass_ee_cutbasedTight ||
          pass_ee_fall17v2mvanoisowp90 ||
          pass_ee_fall17v2mvanoisowp80 ||
          pass_ee_fall17v2mvaisowp90 ||
          pass_ee_fall17v2mvaisowp80 ||
          pass_ee_hzzmvaisowpHZZ
        ) nacc++;
      }

      bool pass_mumu_cutbasedMedium = false;
      bool pass_mumu_cutbasedTight = false;
      if (is_mumu){
        pass_mumu_cutbasedMedium = true;
        pass_mumu_cutbasedTight = true;
        for (auto const& muon:looseMuons){
          auto const& extras = muon->extras;
          float relpfiso03 = MuonSelectionHelpers::relPFIso_DR0p3(*muon);
          pass_mumu_cutbasedMedium &= ((extras.POG_selector_bits & Muon::CutBasedIdMediumPrompt) == Muon::CutBasedIdMediumPrompt) && relpfiso03<MuonSelectionHelpers::isoThr_medium;
          pass_mumu_cutbasedTight &= ((extras.POG_selector_bits & Muon::CutBasedIdTight) == Muon::CutBasedIdTight) && relpfiso03<MuonSelectionHelpers::isoThr_tight;
        }
        if (
          pass_mumu_cutbasedMedium ||
          pass_mumu_cutbasedTight
          ) nacc++;
      }

      /*
      if (nacc%100 == 0 || n_ak4jets_tight>0){
        MELAout << "Accumulated event: " << nacc << '/' << ev << '/' << nevents << endl;
        MELAout << "\t- Number of electrons (raw, precleaning, postcleaning) = ( " << electrons.size() << ", " << nLooseElectrons_precleaning << ", " << nLooseElectrons_postcleaning << " )" << endl;
        MELAout << "\t- Number of muons (raw, loose) = ( " << muons.size() << ", " << looseMuons.size() << " )" << endl;
        if (is_ee){
          if (pass_ee_cutbasedMedium) MELAout << "\t- pass_ee_cutbasedMedium" << endl;
          if (pass_ee_cutbasedTight) MELAout << "\t- pass_ee_cutbasedTight" << endl;
          if (pass_ee_fall17v2mvanoisowp90) MELAout << "\t- pass_ee_fall17v2mvanoisowp90" << endl;
          if (pass_ee_fall17v2mvanoisowp80) MELAout << "\t- pass_ee_fall17v2mvanoisowp80" << endl;
          if (pass_ee_fall17v2mvaisowp90) MELAout << "\t- pass_ee_fall17v2mvaisowp90" << endl;
          if (pass_ee_fall17v2mvaisowp80) MELAout << "\t- pass_ee_fall17v2mvaisowp80" << endl;
          if (pass_ee_hzzmvaisowpHZZ) MELAout << "\t- pass_ee_hzzmvaisowpHZZ" << endl;
        }
        else{
          if (pass_mumu_cutbasedMedium) MELAout << "\t- pass_mumu_cutbasedMedium" << endl;
          if (pass_mumu_cutbasedTight) MELAout << "\t- pass_mumu_cutbasedTight" << endl;
        }
        MELAout << "\t- Leading lepton: [" << leadingLepton->pdgId() << "] ( " << leadingLepton->pt() << ", " << leadingLepton->eta() << ", " << leadingLepton->phi() << ", " << leadingLepton->m() << ")" << endl;
        MELAout << "\t- Subleading lepton: [" << subleadingLepton->pdgId() << "] ( " << subleadingLepton->pt() << ", " << subleadingLepton->eta() << ", " << subleadingLepton->phi() << ", " << subleadingLepton->m() << ")" << endl;
        MELAout << "\t- Njets = " << n_ak4jets_tight << ", mll = " << mll << ", MET = " << pfmet->pt() << endl;
      }
      */
      // Fill histograms
      // Enclosed around braces to localize it_hist
      {
        auto it_hist = sample.hlist_1D.begin();
        for (int ichannel=0; ichannel<nchannels; ichannel++){
          if (ichannel!=procsel) continue;
          bool isCorrectChannel = (
            ichannel==-1
            || (ichannel==0 && is_ee)
            || (ichannel==1 && is_mumu)
            || (ichannel==2 && is_emu)
            || (ichannel==3 && (is_ee || is_mumu))
            );

          for (auto const& cutset:cutsets){
            bool doFill=true;
            for (auto it_cut = cutset.cbegin(); it_cut != cutset.cend(); it_cut++){
              TString const& cutvar = it_cut->cutvar;
              float cutval=0;
              if (cutvar == "Nj") cutval = n_ak4jets_tight;
              else if (cutvar == "pfmet_pTmiss") cutval = pfmet->pt();
              doFill &= it_cut->testCut(cutval);
            }

            if (isCorrectChannel && doFill) it_hist->hist.Fill(mll, wgt); it_hist++;
            if (procsel == 0){ // Z->ee
              if (isCorrectChannel && doFill && pass_ee_cutbasedMedium) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_cutbasedTight) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_fall17v2mvanoisowp90) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_fall17v2mvanoisowp80) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_fall17v2mvaisowp90) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_fall17v2mvaisowp80) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_ee_hzzmvaisowpHZZ) it_hist->hist.Fill(mll, wgt); it_hist++;
            }
            else if (procsel == 1){ // Z->mumu
              if (isCorrectChannel && doFill && pass_mumu_cutbasedMedium) it_hist->hist.Fill(mll, wgt); it_hist++;
              if (isCorrectChannel && doFill && pass_mumu_cutbasedTight) it_hist->hist.Fill(mll, wgt); it_hist++;
            }
          }
        } // End loop over channels
      } // End fill
    } // End loop over events

    sample.writeHistograms();

    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples
}

void makePlots(
  int doZZWW, bool doOldSelection, bool usePuppiMETForSelection, TString strdate="",
  bool useLogY=true, bool isStacked=true, int nfoci=0, int ifocus=0
){
  if (nfoci>0 && (ifocus>=nfoci || ifocus<0)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  //SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  //constexpr int nchannels = 4;

  TString const cinput_main = "output/" + strdate + (doZZWW==0 ? "/ZZCuts" : "/WWCuts") + (doOldSelection ? "/OldCuts" : "/NewCuts") + (usePuppiMETForSelection ? "/PUPPIMETCuts" : "/PFMETCuts");;
  TString coutput_main = cinput_main + "/Plots";
  if (isStacked) coutput_main += "/Stacked";
  else coutput_main += "/Unstacked";
  if (useLogY) coutput_main += "/LogY";
  else coutput_main += "/Linear";
  if (nfoci>0) coutput_main += Form("/Divide%i/Focus%i", nfoci, ifocus);
  gSystem->mkdir(coutput_main, true);

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


  std::vector<std::string> sampleList;
  //sampleList.emplace_back("DY_2l_M_10to50");
  //sampleList.emplace_back("DY_2l_M_50");
  sampleList.emplace_back("DY_2l_M_50_HT");
  sampleList.emplace_back("TT2L2Nu");
  sampleList.emplace_back("ZZ2L2Nu");
  sampleList.emplace_back("WW2L2Nu");
  sampleList.emplace_back("ggZZ_BSI");
  sampleList.emplace_back("VBF_BSI");
  sampleList.emplace_back("ggZZ_Sig");
  sampleList.emplace_back("VBF_Sig");

  std::vector<std::string> samplePlottedList;

  std::unordered_map<TString, std::vector<TH1F*>> sample_hist_map;

  bool firstFile = true;
  std::vector<TString> hnames;
  std::vector<TFile*> finputlist;
  for (auto const& sample:sampleList){
    TString cinput = cinput_main + '/' + sample.data() + ".root";
    if (!HostHelpers::FileReadable(cinput)) continue;

    TFile* finput = TFile::Open(cinput, "read");
    finputlist.push_back(finput);

    std::vector<TH1F*> hlist;
    HelperFunctions::extractHistogramsFromDirectory(finput, hlist);
    if (firstFile){
      for (TH1F* hh:hlist) hnames.push_back(hh->GetName());
    }
    if (hnames.size() == hlist.size()){
      sample_hist_map[sample] = hlist;
      samplePlottedList.push_back(sample);
    }

    if (firstFile) firstFile=false;
  }

  for (size_t iplot=0; iplot<hnames.size(); iplot++){
    int nbins=sample_hist_map[samplePlottedList.front()].at(0)->GetNbinsX();
    int binstart=1;
    int binend=nbins;
    if (nfoci>0){
      int bininc = nbins / nfoci;
      binstart = (ifocus * bininc) + 1;
      binend = (ifocus==nfoci-1 ? nbins : binstart + bininc);
    }
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);
      hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetBinLowEdge(binstart), hist->GetXaxis()->GetBinUpEdge(binend));
    }
    if (isStacked){
      for (size_t is=0; is<samplePlottedList.size(); is++){
        auto& sample = samplePlottedList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        bool isSignal = sample.find("Sig")!=std::string::npos;
        for (size_t js=is+1; js<samplePlottedList.size(); js++){
          auto& sample_j = samplePlottedList.at(js);
          TH1F* hist_j = sample_hist_map[sample_j].at(iplot);
          bool isSignal_j = sample_j.find("Sig")!=std::string::npos;
          if (isSignal_j!=isSignal) continue;
          hist->Add(hist_j);
        }
      }
    }
    else{
      for (size_t is=0; is<samplePlottedList.size(); is++){
        auto& sample = samplePlottedList.at(is);
        TH1F* hist = sample_hist_map[sample].at(iplot);
        double inthist = hist->Integral(binstart, binend);
        hist->Scale(1. / inthist);
      }
    }

    size_t nplottables=0;
    double ymin = 0;
    double ymax = -1;
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (doZZWW==0 && TString(hist->GetName()).Contains("emu")) continue;
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
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);
      hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    TString canvasname = hnames.at(iplot);
    TString canvasname_core = (useLogY ? (isStacked ? "cLogY_Stacked_" : "cLogY_") : (isStacked ? "c_Stacked_" : "c_"));
    if (nfoci>0) canvasname_core += Form("Focus_%i_of_%i_", ifocus, nfoci);
    HelperFunctions::replaceString(canvasname, "h1D_", "");
    canvasname = canvasname_core + canvasname;
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
    TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", 59.7, 13);
    text = pt->AddText(0.82, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    bool firstHist = true;
    std::vector<TString> selectionList;
    for (size_t is=0; is<samplePlottedList.size(); is++){
      auto& sample = samplePlottedList.at(is);
      TH1F* hist = sample_hist_map[sample].at(iplot);

      if (doZZWW==0 && TString(hist->GetName()).Contains("emu")) continue;

      std::vector<TString> tmplist;
      TString htitle = hist->GetTitle();
      HelperFunctions::splitOptionRecursive(htitle, tmplist, '|');
      TString hlabel = tmplist.front();

      hist->SetTitle("");
      if (sample.find("Sig")!=string::npos){
        hist->SetFillStyle(3354);
        hist->SetFillColor(hist->GetLineColor());
      }
      else{
        hist->SetLineColor(kBlack);
        hist->SetFillStyle(1001);
      }
      if (sample=="DY_2l_M_50_HT") hlabel = "DY (HT-binned)";
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
}
