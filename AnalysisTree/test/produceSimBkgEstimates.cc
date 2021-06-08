#include <thread>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


constexpr bool useJetOverlapStripping=false;



void produceSystematicsReweighting_Pythia(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  using namespace SystematicsHelpers;

  if (!(theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp)) return;
  if (period=="2016") return; // We are already running to patch 2016!

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  MELAout << "PS weights will be used for Pythia reweighting..." << endl;

  // Acquire the nominal / syst tree pairs
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;

  // Create the output file
  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString const coutput_main = "output/SystematicsReweighting/" + strdate + "/" + period;

  // Set output directory
  gSystem->mkdir(coutput_main, true);

  TString stroutput = Form("%s/%s_%s.root", coutput_main.Data(), strSampleSet.Data(), strSyst.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");

  foutput->cd();

  ExtendedBinning binning_mass({ 0, 100, 200, 400, 600, 800, 1000 }, "KFactor_EW_NLO_qqVV_Bkg_arg_mass");
  ExtendedBinning binning_pt({ 0, 20, 50, 100, 200, 300 }, "KFactor_EW_NLO_qqVV_Bkg_arg_pthat");
  ExtendedBinning binning_that({ -1000, -300, -150, -100, -50, -20, 0 }, "KFactor_EW_NLO_qqVV_Bkg_arg_that");

  TH3F h_nominal(
    "h_nominal", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
    binning_pt.getNbins(), binning_pt.getBinning(),
    binning_that.getNbins(), binning_that.getBinning()
  );
  TH3F h_syst(
    "h_syst", "",
    binning_mass.getNbins(), binning_mass.getBinning(),
    binning_pt.getNbins(), binning_pt.getBinning(),
    binning_that.getNbins(), binning_that.getBinning()
  );
  h_nominal.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_syst.GetXaxis()->SetTitle(binning_mass.getLabel());
  h_nominal.GetYaxis()->SetTitle(binning_pt.getLabel());
  h_syst.GetYaxis()->SetTitle(binning_pt.getLabel());
  h_nominal.GetZaxis()->SetTitle(binning_that.getLabel()); h_nominal.Sumw2();
  h_syst.GetZaxis()->SetTitle(binning_that.getLabel()); h_syst.Sumw2();

  curdir->cd();

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(false);
  genInfoHandler.setAcquireGenAK4Jets(false);

  for (auto const& sdir:sampledirs){
    TString sid = SampleHelpers::getSampleIdentifier(sdir);
    
    MELAout << "Looping over " << sid << "..." << endl;

    BaseTree* sample_tree = new BaseTree(SampleHelpers::getDatasetFileName(sdir), "cms3ntuple/Events", "", "");
    sample_tree->sampleIdentifier = sid;

    genInfoHandler.bookBranches(sample_tree);

    sample_tree->silenceUnused();

    int nEntries = 0;
    TH3F* htmp_nominal = (TH3F*) h_nominal.Clone("htmp_nom"); htmp_nominal->Reset("ICESM");
    TH3F* htmp_syst = (TH3F*) h_nominal.Clone("htmp_syst"); htmp_syst->Reset("ICESM");

    double sum_wgts = 0, sum_wgts_syst = 0;
    nEntries = sample_tree->getNEvents();
    genInfoHandler.wrapTree(sample_tree);
    MELAout << "\t- Looping over " << nEntries << " entries:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      HelperFunctions::progressbar(ev, nEntries);
      sample_tree->getEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal);
      auto const& genInfo = genInfoHandler.getGenInfo();

      double genwgt = genInfo->getGenWeight(true);
      double genwgt_syst = genwgt;

      double genwgt_adjustment_NNPDF30 = (genwgt==0. ? 1. : genInfo->getGenWeight(false)/genwgt);
      double genwgt_adjustment_NNPDF30_syst = genwgt_adjustment_NNPDF30;
      if (theGlobalSyst==tPythiaScaleDn) genwgt_adjustment_NNPDF30_syst *= genInfo->extras.PythiaWeight_isr_muR0p25 * genInfo->extras.PythiaWeight_fsr_muR0p25;
      else if (theGlobalSyst==tPythiaScaleUp) genwgt_adjustment_NNPDF30_syst *= genInfo->extras.PythiaWeight_isr_muR4 * genInfo->extras.PythiaWeight_fsr_muR4;

      if (std::abs(genwgt_adjustment_NNPDF30)>10.) genwgt_adjustment_NNPDF30 = 10. / std::abs(genwgt_adjustment_NNPDF30);
      if (std::abs(genwgt_adjustment_NNPDF30_syst)>10.) genwgt_adjustment_NNPDF30_syst = 10. / std::abs(genwgt_adjustment_NNPDF30_syst);

      genwgt *= genwgt_adjustment_NNPDF30;
      genwgt_syst *= genwgt_adjustment_NNPDF30_syst;

      sum_wgts += genwgt;
      sum_wgts_syst += genwgt_syst;

      auto const& var_mass = genInfo->extras.Kfactors.find(binning_mass.getName())->second;
      auto const& var_pt = genInfo->extras.Kfactors.find(binning_pt.getName())->second;
      auto const& var_that = genInfo->extras.Kfactors.find(binning_that.getName())->second;

      htmp_nominal->Fill(
        var_mass, var_pt, var_that,
        genwgt
      );
      htmp_syst->Fill(
        var_mass, var_pt, var_that,
        genwgt_syst
      );
    }
    htmp_nominal->Scale(double(nEntries)/sum_wgts);
    htmp_syst->Scale(double(nEntries)/sum_wgts); // Use sum_wgts here.
    HelperFunctions::wipeOverUnderFlows(htmp_nominal, false, true);
    HelperFunctions::wipeOverUnderFlows(htmp_syst, false, true);
    h_nominal.Add(htmp_nominal, 1.);
    h_syst.Add(htmp_syst, 1.);
    delete htmp_nominal;
    delete htmp_syst;

    delete sample_tree;
  }

  foutput->WriteTObject(&h_nominal);
  foutput->WriteTObject(&h_syst);

  HelperFunctions::conditionalizeHistogram<TH3F>(&h_nominal, 0, nullptr, false, false);
  HelperFunctions::conditionalizeHistogram<TH3F>(&h_syst, 0, nullptr, false, false);

  TH3F* h_ratio = (TH3F*) h_nominal.Clone("h_ratio");
  HelperFunctions::divideHistograms(&h_syst, &h_nominal, h_ratio, false);
  foutput->WriteTObject(h_ratio);

  for (unsigned int isl=0; isl<binning_mass.getNbins(); isl++){
    TH2F* h_nominal_slice = HelperFunctions::getHistogramSlice(&h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal.GetName(), isl));
    TH2F* h_syst_slice = HelperFunctions::getHistogramSlice(&h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst.GetName(), isl));
    TH2F* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));

    foutput->WriteTObject(h_nominal_slice); delete h_nominal_slice;
    foutput->WriteTObject(h_syst_slice); delete h_syst_slice;
    foutput->WriteTObject(h_ratio_slice); delete h_ratio_slice;
  }

  delete h_ratio;

  foutput->Close();

  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}

// Alias function to do various stuff.
void produceSystematicsReweighting(
  TString strSampleSet,
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  produceSystematicsReweighting_Pythia(strSampleSet, period, prodVersion, strdate, theGlobalSyst);
}

void combineSystematicsReweightings(
  TString strSampleSet, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  using namespace SystematicsHelpers;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  if (strdate=="") strdate = HelperFunctions::todaysdate();
  if (!(theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp)) return;
  if (
    !(
      strSampleSet == "qqZZ"
      ||
      strSampleSet == "qqWZ"
      ||
      strSampleSet == "qqWW"
      )
    ) return;

  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();

  std::vector< std::pair<TString, TString> > input_period_sname_pairs;
  if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
    std::vector<TString> periods{ "2017", "2018" };
    for (auto const& period:periods){
      std::vector<TString> slist;
      if (period == "2017"){
        if (strSampleSet == "qqZZ") slist = std::vector<TString>{ "qqZZ_2l2nu_mZ_18-inf", "qqZZ_2l2nu", "qqZZ_4l", "qqZZ_4l_ext" };
        else if (strSampleSet == "qqWZ") slist = std::vector<TString>{ "qqWZ_3lnu_POWHEG_mll_0p1-inf", "qqWZ_3lnu_POWHEG" };
        else if (strSampleSet == "qqWW") slist = std::vector<TString>{ "qqWW_lnu2q", "qqWW_lnu2q_ext", "qqWW_2l2nu", "qqWW_2l2nu_ext" };
      }
      else{
        if (strSampleSet == "qqZZ") slist = std::vector<TString>{ "qqZZ_2l2nu_mZ_18-inf", "qqZZ_2l2nu", "qqZZ_2l2nu_ext", "qqZZ_4l" };
        else if (strSampleSet == "qqWZ") slist = std::vector<TString>{ "qqWZ_3lnu_POWHEG_mll_0p1-inf", "qqWZ_3lnu_POWHEG" };
        else if (strSampleSet == "qqWW") slist = std::vector<TString>{ "qqWW_lnu2q", "qqWW_2l2nu" };
      }

      for (auto const& ss:slist) input_period_sname_pairs.emplace_back(period, ss);
    }
  }

  std::vector< std::pair<TString, TString> > output_period_sname_pairs;
  if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
    std::vector<TString> slist;
    if (strSampleSet == "qqZZ") slist = std::vector<TString>{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext", "qqZZ_4l", "qqZZ_4l_ext" };
    else if (strSampleSet == "qqWZ") slist = std::vector<TString>{ "qqWZ_3lnu_POWHEG_mll_0p1-inf", "qqWZ_3lnu_POWHEG" };
    else if (strSampleSet == "qqWW") slist = std::vector<TString>{ "qqWW_lnu2q", "qqWW_2l2nu" };

    for (auto const& ss:slist) output_period_sname_pairs.emplace_back("2016", ss);
  }

  //foutput->cd();

  TH3F* h_nominal=nullptr;
  TH3F* h_syst=nullptr;

  TH2F* h_nominal_reduced=nullptr;
  TH2F* h_syst_reduced=nullptr;

  TH1F* h_nominal_1D=nullptr;
  TH1F* h_syst_1D=nullptr;

  for (auto const& period_sname_pair:input_period_sname_pairs){
    TString const cinput_main = "output/SystematicsReweighting/" + strdate + "/" + period_sname_pair.first;

    TString strinput = Form("%s/%s_%s.root", cinput_main.Data(), period_sname_pair.second.Data(), strSyst.Data());
    if (!SampleHelpers::checkFileOnWorker(strinput)) continue;

    TFile* finput = TFile::Open(strinput, "read");
    finput->cd();

    if (theGlobalSyst==tPythiaScaleDn || theGlobalSyst==tPythiaScaleUp){
      TH3F* htmp_nominal = (TH3F*) finput->Get("h_nominal");
      TH3F* htmp_syst = (TH3F*) finput->Get("h_syst");
      TH3F* hratio = (TH3F*) finput->Get("h_ratio");

      bool hasValidRatio = false;
      for (int ix=1; ix<=hratio->GetNbinsX(); ix++){
        for (int iy=1; iy<=hratio->GetNbinsY(); iy++){
          for (int iz=1; iz<=hratio->GetNbinsZ(); iz++){
            double bc = hratio->GetBinContent(ix, iy, iz);
            if (std::abs(bc-1.)>1e-4 && std::abs(bc)>1e-4){
              hasValidRatio = true;
              break;
            }
          }
        }
      }
      if (!hasValidRatio){
        output_period_sname_pairs.push_back(period_sname_pair);
        MELAout << period_sname_pair << " has invalid ratios." << endl;
      }
      else{
        curdir->cd();
        if (!h_nominal) h_nominal = (TH3F*) htmp_nominal->Clone(htmp_nominal->GetName());
        else h_nominal->Add(htmp_nominal, 1.);
        if (!h_syst) h_syst = (TH3F*) htmp_syst->Clone(htmp_syst->GetName());
        else h_syst->Add(htmp_syst, 1.);
      }
    }

    finput->Close();
    curdir->cd();
  }

  TH3F* h_ratio = nullptr;
  std::vector<TH2F*> h_nominal_slices;
  std::vector<TH2F*> h_syst_slices;
  std::vector<TH2F*> h_ratio_slices;
  if (h_nominal){
    h_ratio = (TH3F*) h_nominal->Clone("h_ratio");
    HelperFunctions::divideHistograms(h_syst, h_nominal, h_ratio, false);

    for (int isl=0; isl<h_nominal->GetNbinsX(); isl++){
      TH2F* h_nominal_slice = HelperFunctions::getHistogramSlice(h_nominal, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_nominal->GetName(), isl));
      TH2F* h_syst_slice = HelperFunctions::getHistogramSlice(h_syst, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_syst->GetName(), isl));
      TH2F* h_ratio_slice = HelperFunctions::getHistogramSlice(h_ratio, 1, 2, isl+1, isl+1, Form("%s_slice%u", h_ratio->GetName(), isl));

      h_nominal_slices.push_back(h_nominal_slice);
      h_syst_slices.push_back(h_syst_slice);
      h_ratio_slices.push_back(h_ratio_slice);
    }
  }

  for (auto const& period_sname_pair:output_period_sname_pairs){
    TString const& strSample = period_sname_pair.second;
    TString strOutPeriod = Form("Combined/%s", period_sname_pair.first.Data());
    TString const coutput_main = "output/SystematicsReweighting/" + strdate + "/" + strOutPeriod;


    // Set output directory
    gSystem->mkdir(coutput_main, true);

    TString stroutput = Form("%s/%s_%s.root", coutput_main.Data(), strSample.Data(), strSyst.Data());
    SampleHelpers::addToCondorTransferList(stroutput);

    TFile* foutput = TFile::Open(stroutput, "recreate");

    foutput->WriteTObject(h_nominal);
    foutput->WriteTObject(h_syst);
    foutput->WriteTObject(h_ratio);

    for (auto const& hh:h_ratio_slices) foutput->WriteTObject(hh);
    for (auto const& hh:h_syst_slices) foutput->WriteTObject(hh);
    for (auto const& hh:h_nominal_slices) foutput->WriteTObject(hh);

    foutput->Close();
  }

  curdir->cd();

  for (auto& hh:h_ratio_slices) delete hh;
  for (auto& hh:h_syst_slices) delete hh;
  for (auto& hh:h_nominal_slices) delete hh;

  delete h_ratio;
  delete h_syst;
  delete h_nominal;

  curdir->cd();
}



// Dummy for now, but can be extended to veto certain samples for certain systematics
bool checkSystematicForSampleGroup(TString const& sgroup, SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;
  switch (syst){
  case tEWDn:
  case tEWUp:
    return (sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW"));
  default:
    return true;
  }
}


void getMCSampleSet_ZZTo2L2Nu(std::vector< std::pair< TString, std::vector<TString> > >& sampleSpecs);
void getMCSampleSet_ZWTo3L1Nu(std::vector< std::pair< TString, std::vector<TString> > >& sampleSpecs);

void getMCSampleSet(std::vector< std::pair< TString, std::vector<TString> > >& sampleSpecs){
  OffshellCutflow::FinalStateType const& fstype = OffshellCutflow::activeFinalState;
  if (fstype==OffshellCutflow::fs_ZZ_2l2nu) getMCSampleSet_ZZTo2L2Nu(sampleSpecs);
  else if (fstype==OffshellCutflow::fs_ZW_3l1nu) getMCSampleSet_ZWTo3L1Nu(sampleSpecs);
  else{
    MELAerr << "getMCSampleSet: Final state " << fstype << " is not defined." << endl;
    exit(1);
  }
}

void getMCSampleDirs(
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> >& strsamples,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  std::unordered_map<TString, TString>& externalCorrections,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  using namespace SystematicsHelpers;

  TString const strSyst_real = SystematicsHelpers::getSystName(theGlobalSyst).data();
  std::vector<SystematicsHelpers::SystematicVariationTypes> const disallowedSysts{
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,
    tEWDn, tEWUp,

    eEleEffDn, eEleEffUp,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,

    eMuEffDn, eMuEffUp,
    eMuEffStatDn, eMuEffStatUp,
    eMuEffSystDn, eMuEffSystUp,
    eMuEffAltMCDn, eMuEffAltMCUp,

    ePhoEffDn, ePhoEffUp,

    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    ePUDn, ePUUp,
    eL1PrefiringDn, eL1PrefiringUp,

    eTriggerEffDn, eTriggerEffUp
  };
  if (HelperFunctions::checkListVariable(disallowedSysts, theGlobalSyst)) theGlobalSyst = sNominal;

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString period = SampleHelpers::getDataPeriod();

  std::vector< std::pair< TString, std::vector<TString> > > sampleSpecs;
  getMCSampleSet(sampleSpecs);

  TString cinput_main =
    TString("AK4Jets")
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  cinput_main = cinput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  cinput_main = cinput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  cinput_main = cinput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  cinput_main = cinput_main + "/" + period;

  for (auto const& s:sampleSpecs){
    if (!checkSystematicForSampleGroup(s.first, theGlobalSyst)) continue;

    std::vector<TString> sdirs;
    std::vector<std::pair<TString, TString>> sname_dir_pairs;
    for (auto const& strSampleSet:s.second){
      SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sdirs);
      TString strExtCorr;
      
      TString strinput_customSyst = ANALYSISTREEPKGDATAPATH + "ScaleFactors/SystematicsCustomReweighting/";
      HostHelpers::ExpandEnvironmentVariables(strinput_customSyst);
      strinput_customSyst = strinput_customSyst + Form("%i/%s_%s%s", SampleHelpers::getDataYear(), strSampleSet.Data(), strSyst_real.Data(), ".root");
      MELAout << "getMCSampleDirs: Checking for a corrections file named " << strinput_customSyst << "...";
      if (HostHelpers::FileExists(strinput_customSyst)){
        strExtCorr = strinput_customSyst;
        MELAout << " FOUND...";
      }
      else MELAout << " DNE";
      MELAout << endl;

      if (strExtCorr!="") externalCorrections[sdirs.back()] = strExtCorr;
    }
    if (sdirs.size()!=s.second.size()){
      MELAerr << "getMCSampleDirs: Size of sdirs " << sdirs.size() << " does not match the input sample set size = " << s.second.size() << " for group " << s.first << "." << endl;
      exit(1);
    }
    sname_dir_pairs.reserve(sdirs.size());
    for (auto const& sname:sdirs){
      TString cinput = SampleHelpers::getSampleIdentifier(sname);
      HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
      HelperFunctions::replaceString(cinput, "_MINIAOD", "");
      cinput = cinput + "_" + strSyst + "*.root";
      cinput = cinput_main + "/" + cinput;
      sname_dir_pairs.emplace_back(sname, cinput);
    }
    strsamples.emplace_back(s.first, sname_dir_pairs);
  }
}


using namespace SystematicsHelpers;

#include "produceSimBkgEstimates_ZZTo2L2Nu.h"

#include "produceSimBkgEstimates_ZWTo3L1Nu.h"
