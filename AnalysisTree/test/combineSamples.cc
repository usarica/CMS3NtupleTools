#include "common_includes.h"
#include "BulkReweightingBuilder.h"
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"
#include "TChain.h"
#include "TString.h"


using namespace std;
using namespace MELAStreamHelpers;
using namespace PDGHelpers;

struct SampleSpecs{
  float mass;
  TString path;
  SampleSpecs(float mass_, TString path_) : mass(mass_), path(path_){}
};

void combineSamples(TString inpath, TString outpath){
  std::vector<SampleSpecs> strsamples_POWHEG_ZZ{
    SampleSpecs(200, inpath + "/GluGluHToZZTo2L2Nu_M200_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1"),
    SampleSpecs(300, inpath + "/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1"),
    SampleSpecs(400, inpath + "/GluGluHToZZTo2L2Nu_M400_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(500, inpath + "/GluGluHToZZTo2L2Nu_M500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(700, inpath + "/GluGluHToZZTo2L2Nu_M700_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(800, inpath + "/GluGluHToZZTo2L2Nu_M800_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(600, inpath + "/GluGluHToZZTo2L2Nu_M600_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(900, inpath + "/GluGluHToZZTo2L2Nu_M900_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(1000, inpath + "/GluGluHToZZTo2L2Nu_M1000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(1500, inpath + "/GluGluHToZZTo2L2Nu_M1500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(2500, inpath + "/GluGluHToZZTo2L2Nu_M2500_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(2000, inpath + "/GluGluHToZZTo2L2Nu_M2000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2"),
    SampleSpecs(3000, inpath + "/GluGluHToZZTo2L2Nu_M3000_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2")
  };

  TString strxsec = "xsec";
  TString strWgtNominal = "genHEPMCweight_default";
  TString strbinningvar = "LHECandMass";
  ExtendedBinning LHECandMassBinning(strbinningvar);
  for (unsigned int is=0; is<strsamples_POWHEG_ZZ.size()-1; is++){
    if (strsamples_POWHEG_ZZ.at(is).mass==strsamples_POWHEG_ZZ.at(is+1).mass) continue;
    float boundary = (strsamples_POWHEG_ZZ.at(is).mass + strsamples_POWHEG_ZZ.at(is+1).mass)/2.;
    LHECandMassBinning.addBinBoundary(boundary);
  }
  LHECandMassBinning.addBinBoundary(0);
  LHECandMassBinning.addBinBoundary(13000.);

  TString cpsname;
  std::vector<TString> melist;
  std::vector<BaseTree> samplelist; samplelist.reserve(strsamples_POWHEG_ZZ.size());
  bool firstSample = true;
  for (auto const& ss:strsamples_POWHEG_ZZ){
    samplelist.emplace_back(ss.path+"/allevents_modified.root", "cms3ntuple/Events", "", "");
    BaseTree& sample_tree = samplelist.back();

    std::vector<TString> allbranchnames; sample_tree.getValidBranchNamesWithoutAlias(allbranchnames, false);
    for (TString const& bname : allbranchnames){
      if (bname.Contains("p_Gen")){
        if (firstSample && bname.Contains("CPS")) cpsname = bname;
        else if (firstSample) melist.push_back(bname);
        sample_tree.bookBranch<float>(bname, 0.f);
      }
      else if (bname == strbinningvar || bname == strxsec || bname == strWgtNominal) sample_tree.bookBranch<float>(bname, 0.f);
    }

    firstSample = false;
  }

  std::vector<TString> strNominalWeights{ strWgtNominal };
  std::vector<TString> strCrossSectionWeights{ strxsec };
  BulkReweightingBuilder melarewgtBuilder(
    strNominalWeights,
    strCrossSectionWeights,
    ReweightingFunctions::getSimpleWeight,
    ReweightingFunctions::getSimpleWeight
  );
  melarewgtBuilder.rejectNegativeWeights(false);
  melarewgtBuilder.setWeightBinning(LHECandMassBinning);
  for (TString const& strme:melist){
    std::vector<TString> strReweightingWeights{ strme };
    if (cpsname != "") strReweightingWeights.push_back(cpsname);
    melarewgtBuilder.addReweightingWeights(
      strme,
      strReweightingWeights,
      ReweightingFunctions::getSimpleWeight
    );
  }
  for (auto& ss:samplelist) melarewgtBuilder.setupWeightVariables(&ss, 0.999, 250);
  melarewgtBuilder.setupCaches();

  std::vector<TString> newsamplelist;
  for (size_t is = 0; is<samplelist.size(); is++){
    BaseTree& ss = samplelist.at(is);
    // Unmute all branches that were silenced before
    ss.unmuteAllBranches();
    TString sample_path = strsamples_POWHEG_ZZ.at(is).path + "/allevents_modified_modifiedMEs.root";
    newsamplelist.push_back(sample_path);

    TFile* foutput = TFile::Open(sample_path, "recreate");
    foutput->cd();
    TDirectory* curdir = gDirectory;
    TDirectory* subdir = foutput->mkdir("cms3ntuple");
    subdir->cd();

    TTree* outputTree = ss.getSelectedTree()->CloneTree(0);
    outputTree->SetAutoSave(0);

    int nevents = ss.getSelectedNEvents();
    for (int ev=0; ev<nevents; ev++){
      HelperFunctions::progressbar(ev, nevents);
      ss.getSelectedEvent(ev);

      for (unsigned int ime=0; ime<melist.size(); ime++){
        float* MEval;
        ss.getValRef(melist.at(ime), MEval);
        float finalWeight = melarewgtBuilder.getFinalEventWeight(&ss, melist.at(ime));
        float nominalWeights = melarewgtBuilder.eval_nominalweights(&ss) * melarewgtBuilder.eval_xsecweights(&ss);
        if (nominalWeights!=0.) *MEval = finalWeight / nominalWeights;
      }

      outputTree->Fill();
    }

    subdir->WriteTObject(outputTree);
    subdir->Close();
    foutput->Close();
  }

  {
    TChain* intree = new TChain("cms3ntuple/Events");
    for (auto const& s:newsamplelist) intree->Add(s);

    gSystem->mkdir(outpath, true);

    TFile* foutput = TFile::Open(outpath+"/allevents.root", "recreate");
    foutput->cd();
    TDirectory* curdir = gDirectory;
    TDirectory* subdir = foutput->mkdir("cms3ntuple");
    subdir->cd();

    TTree* outputTree = intree->CloneTree(-1, "fast");

    subdir->WriteTObject(outputTree);
    subdir->Close();
    foutput->Close();

    delete intree;
  }
}
