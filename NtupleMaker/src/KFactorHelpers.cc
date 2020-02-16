#include <cassert>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/KFactorHelpers.h>
#include <CMSDataTools/AnalysisTree/interface/HostHelpersCore.h>
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"


using namespace std;
using namespace MELAStreamHelpers;


namespace KFactorHelpers{

  const std::string KFactorHandler_QCD_ggZZ_Sig::KFactorArgName = "KFactor_QCD_ggZZ_Sig_arg";
  KFactorHandler_QCD_ggZZ_Sig::KFactorHandler_QCD_ggZZ_Sig(int year) :
    KFactorHandlerBase(),
    kfFile_NNLO(nullptr),
    kfFile_NLO(nullptr)
  {
    TString strKFactorDir = "${CMSSW_BASE}/src/CMS3/NtupleMaker/data/Kfactors/";
    HostHelpers::ExpandEnvironmentVariables(strKFactorDir);

    TString strSqrts="";
    switch (year){
    case 2015:
    case 2016:
    case 2017:
    case 2018:
      strSqrts="13TeV";
      break;
    default:
      MELAerr << "KFactorHandler_QCD_ggZZ_Sig::KFactorHandler_QCD_ggZZ_Sig: Year " << year << " is not implemented. Aborting..." << endl;
      assert(0);
    }
    strKFactorDir += strSqrts;

    TString strKFactorNNLOFile = strKFactorDir + "/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root";
    TString strKFactorNLOFile = strKFactorDir + "/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root";

    TDirectory* curdir = gDirectory;
    kfFile_NNLO = std::make_shared<TFile>(strKFactorNNLOFile, "read"); curdir->cd();
    kfFile_NLO = std::make_shared<TFile>(strKFactorNLOFile, "read"); curdir->cd();

    this->setup();
  }
  KFactorHandler_QCD_ggZZ_Sig::KFactorHandler_QCD_ggZZ_Sig(KFactorHandler_QCD_ggZZ_Sig const& other) :
    KFactorHandlerBase(other),
    kfFile_NNLO(other.kfFile_NNLO),
    kfFile_NLO(other.kfFile_NLO)
  {
    // Do not copy vectors of splines, get them from scratch again
    this->setup();
  }
  void KFactorHandler_QCD_ggZZ_Sig::setup(){
    TDirectory* curdir = gDirectory;

    kfFile_NNLO->cd();
    sp_NNLO.reserve(9);
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_Nominal"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_QCDScaleDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_QCDScaleUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFScaleDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFScaleUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFReplicaDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_PDFReplicaUp"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_AsDn"));
    sp_NNLO.push_back((TSpline3*) kfFile_NNLO->Get("sp_kfactor_AsUp"));
    for (size_t ikf=0; ikf<sp_NNLO.size(); ikf++){
      if (!sp_NNLO.at(ikf)) throw cms::Exception(Form("KFactorHandler_QCD_ggZZ_Sig::setup: NNLO K factor at location %lu cannot be found.", ikf));
    }

    kfFile_NLO->cd();
    sp_NLO.reserve(9);
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_Nominal"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_QCDScaleDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_QCDScaleUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFScaleDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFScaleUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFReplicaDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_PDFReplicaUp"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_AsDn"));
    sp_NLO.push_back((TSpline3*) kfFile_NLO->Get("sp_kfactor_AsUp"));
    for (size_t ikf=0; ikf<sp_NLO.size(); ikf++){
      if (!sp_NLO.at(ikf)) throw cms::Exception(Form("KFactorHandler_QCD_ggZZ_Sig::setup: NLO K factor at location %lu cannot be found.", ikf));
    }

    curdir->cd();
  }
  void KFactorHandler_QCD_ggZZ_Sig::eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map){
    static const std::vector<std::string> kfactornames{
      "KFactor_QCD_NNLO_ggZZ_Sig_Nominal",
      "KFactor_QCD_NNLO_ggZZ_Sig_QCDScaleDn",
      "KFactor_QCD_NNLO_ggZZ_Sig_QCDScaleUp",
      "KFactor_QCD_NNLO_ggZZ_Sig_PDFScaleDn",
      "KFactor_QCD_NNLO_ggZZ_Sig_PDFScaleUp",
      "KFactor_QCD_NNLO_ggZZ_Sig_PDFReplicaDn",
      "KFactor_QCD_NNLO_ggZZ_Sig_PDFReplicaUp",
      "KFactor_QCD_NNLO_ggZZ_Sig_AsDn",
      "KFactor_QCD_NNLO_ggZZ_Sig_AsUp",
      "KFactor_QCD_NLO_ggZZ_Sig_Nominal",
      "KFactor_QCD_NLO_ggZZ_Sig_QCDScaleDn",
      "KFactor_QCD_NLO_ggZZ_Sig_QCDScaleUp",
      "KFactor_QCD_NLO_ggZZ_Sig_PDFScaleDn",
      "KFactor_QCD_NLO_ggZZ_Sig_PDFScaleUp",
      "KFactor_QCD_NLO_ggZZ_Sig_PDFReplicaDn",
      "KFactor_QCD_NLO_ggZZ_Sig_PDFReplicaUp",
      "KFactor_QCD_NLO_ggZZ_Sig_AsDn",
      "KFactor_QCD_NLO_ggZZ_Sig_AsUp"
    };
    assert(kfactornames.size()==(this->sp_NNLO.size() + this->sp_NLO.size()));

    auto it_arg = kfactors_map.find(KFactorHandler_QCD_ggZZ_Sig::KFactorArgName);
    if (it_arg == kfactors_map.end()) throw cms::Exception(Form("KFactorHandler_QCD_ggZZ_Sig::eval: K factor evaluation argument, candidate mass with name %s, cannot be found.", KFactorHandler_QCD_ggZZ_Sig::KFactorArgName.data()));
    float const& kfactor_arg = it_arg->second;

    float kfactor_denominator = 1;
    switch (denominator){
    case kf_QCD_NNLO_GGZZ_SIG:
      kfactor_denominator = sp_NNLO.front()->Eval(kfactor_arg);
      break;
    case kf_QCD_NLO_GGZZ_SIG:
      kfactor_denominator = sp_NLO.front()->Eval(kfactor_arg);
      break;
    default:
      break;
    }

    size_t ikf=0;
    if (type == kf_QCD_NNLO_GGZZ_SIG){
      for (auto const& spkf:sp_NNLO){
        kfactors_map[kfactornames.at(ikf)] = spkf->Eval(kfactor_arg) / kfactor_denominator;
        ikf++;
      }
    }
    ikf = sp_NNLO.size();
    if (type == kf_QCD_NLO_GGZZ_SIG){
      for (auto const& spkf:sp_NLO){
        kfactors_map[kfactornames.at(ikf)] = spkf->Eval(kfactor_arg) / kfactor_denominator;
        ikf++;
      }
    }
  }

}
