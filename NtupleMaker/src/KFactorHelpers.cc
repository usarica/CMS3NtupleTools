#include <cassert>
#include <sstream>
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
  void KFactorHandler_QCD_ggZZ_Sig::eval(KFactorHelpers::KFactorType type, KFactorHelpers::KFactorType denominator, std::unordered_map<std::string, float>& kfactors_map) const{
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


  //const std::string KFactorHandler_EW_qqVV_Bkg::KFactorArgName = "KFactor_EW_qqVV_Sig_arg";
  KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg() :
    KFactorHandlerBase()
  {
    TString strKFactorDir = "${CMSSW_BASE}/src/CMS3/NtupleMaker/data/Kfactors/";
    HostHelpers::ExpandEnvironmentVariables(strKFactorDir);

    TString strKFactorFile_ZZ = strKFactorDir + "/ZZ_EWCorrections.dat"; readTableFromFile(strKFactorFile_ZZ, table_ZZ);
    TString strKFactorFile_WZ = strKFactorDir + "/WZ_EWCorrections.dat"; readTableFromFile(strKFactorFile_WZ, table_WZ);
    TString strKFactorFile_WW = strKFactorDir + "/WW_EWCorrections.dat"; readTableFromFile(strKFactorFile_WW, table_WW);

    this->setup();
  }
  KFactorHandler_EW_qqVV_Bkg::KFactorHandler_EW_qqVV_Bkg(KFactorHandler_EW_qqVV_Bkg const& other) :
    KFactorHandlerBase(other),
    table_ZZ(other.table_ZZ),
    table_WZ(other.table_WZ),
    table_WW(other.table_WW)
  {
    this->setup();
  }
  void KFactorHandler_EW_qqVV_Bkg::readTableFromFile(TString const& fname, std::vector< std::vector<double> >& table){
    ifstream fin(fname.Data(), ios_base::in);
    while (!fin.eof()){
      std::string strline;
      std::getline(fin, strline);
      if (strline.empty()) continue;

      table.push_back(std::vector<double>());
      std::vector<double>& tbl_line = table.back(); tbl_line.reserve(5);

      std::stringstream ss(strline);
      std::string strtmp;
      while (ss >> strtmp){
        if (strtmp.empty()) continue;
        tbl_line.push_back(std::stod(strtmp));
      }
    }
  }
  void KFactorHandler_EW_qqVV_Bkg::findTableEntry(double const& shat, double const& that, std::vector< std::vector<double> > const& table, std::vector< std::vector<std::vector<double>>::const_iterator > const& table_sqrtShatBegin, double& val_uc, double& val_ds, double& val_b){
    std::vector<std::vector<double>>::const_iterator const table_end = table.cend();

    double bestShatDiff = -1;
    std::vector<std::vector<double>>::const_iterator it_bestShat = table_end;
    std::vector<std::vector<double>>::const_iterator itNext_bestShat = table_end;
    for (size_t is=0; is<table_sqrtShatBegin.size()-1;is++){
      auto const& itFirst = table_sqrtShatBegin.at(is);
      auto const& itSecond = table_sqrtShatBegin.at(is+1);
      double tmpShatDiff = std::abs(itFirst->front() - std::sqrt(std::abs(shat)));
      if (bestShatDiff<0. || tmpShatDiff < bestShatDiff){
        bestShatDiff = tmpShatDiff;
        it_bestShat = itFirst;
        itNext_bestShat = itSecond;
      }
    }
    assert(it_bestShat != table_end);

    std::vector<std::vector<double>>::const_iterator it_That = it_bestShat;
    std::vector<std::vector<double>>::const_iterator it_bestThat = table_end;
    double bestThatDiff = -1;
    while (it_That != itNext_bestShat){
      double tmpThatDiff = std::abs(it_That->at(1) - that);
      if (bestThatDiff<0. || tmpThatDiff<bestThatDiff){
        bestThatDiff = tmpThatDiff;
        it_bestThat = it_That;
      }
      it_That++;
    }
    assert(it_bestThat != table_end);

    val_uc = it_bestThat->at(2);
    val_ds = it_bestThat->at(3);
    val_b = it_bestThat->at(4);
  }

  void KFactorHandler_EW_qqVV_Bkg::setup(){
    double sqrtShat;

    sqrtShat = -1;
    for (std::vector<std::vector<double>>::const_iterator it = table_ZZ.cbegin(); it!=table_ZZ.cend(); it++){
      double const& tmpSqrtShat = it->front();
      if (tmpSqrtShat!=sqrtShat){
        sqrtShat = tmpSqrtShat;
        table_sqrtShatBegin_ZZ.push_back(it);
      }
    }
    table_sqrtShatBegin_ZZ.push_back(table_ZZ.cend());

    sqrtShat = -1;
    for (std::vector<std::vector<double>>::const_iterator it = table_WZ.cbegin(); it!=table_WZ.cend(); it++){
      double const& tmpSqrtShat = it->front();
      if (tmpSqrtShat!=sqrtShat){
        sqrtShat = tmpSqrtShat;
        table_sqrtShatBegin_WZ.push_back(it);
      }
    }
    table_sqrtShatBegin_WZ.push_back(table_WZ.cend());

    sqrtShat = -1;
    for (std::vector<std::vector<double>>::const_iterator it = table_WW.cbegin(); it!=table_WW.cend(); it++){
      double const& tmpSqrtShat = it->front();
      if (tmpSqrtShat!=sqrtShat){
        sqrtShat = tmpSqrtShat;
        table_sqrtShatBegin_WW.push_back(it);
      }
    }
    table_sqrtShatBegin_WW.push_back(table_WW.cend());
  }
  void KFactorHandler_EW_qqVV_Bkg::eval(KFactorHelpers::KFactorType type, pat::PackedGenParticleCollection const& genparticles, std::unordered_map<std::string, float>& kfactors_map) const{
    static const std::vector<std::string> kfactornames{
      "KFactor_EW_NLO_qqVV_Bkg_Nominal",
      "KFactor_EW_NLO_qqVV_Bkg_EWScaleDn",
      "KFactor_EW_NLO_qqVV_Bkg_EWScaleUp"
    };

    std::vector< std::vector<double> > const* kfactor_table = nullptr;
    double mPole_V1=-1;
    double mPole_V2=-1;
    switch (type){
    case kf_EW_NLO_QQZZ_BKG:
      kfactor_table = &(this->table_ZZ);
      mPole_V1 = mPole_V2 = PDGHelpers::Zmass;
      break;
    case kf_EW_NLO_QQWZ_BKG:
      kfactor_table = &(this->table_WZ);
      mPole_V1 = PDGHelpers::Wmass;
      mPole_V2 = PDGHelpers::Zmass;
      break;
    case kf_EW_NLO_QQWW_BKG:
      kfactor_table = &(this->table_WW);
      mPole_V1 = mPole_V2 = PDGHelpers::Wmass;
      break;
    default:
      // Do nothing
      break;
    }
    if (!kfactor_table) cms::Exception(Form("KFactorHandler_EW_qqVV_Bkg::eval: No K factor table for type %i.", static_cast<int>(type)));

  }

}
