#include "MuonScaleFactorHandler.h"
#include "ParticleSelectionHelpers.h"
#include "SamplesCore.h"
#include "TDirectory.h"


using namespace std;
using namespace SampleHelpers;
using namespace ParticleSelectionHelpers;


MuonScaleFactorHandler::MuonScaleFactorHandler() :
  ScaleFactorHandlerBase(),
  finput_SF_id(nullptr),
  finput_SF_iso(nullptr),
  h_SF_id(nullptr),
  h_SF_iso(nullptr)
{
  this->setup();
}

MuonScaleFactorHandler::~MuonScaleFactorHandler(){ this->reset(); }


bool MuonScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  // Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF
  if (theDataYear == 2016){
    // ID/Iso./IP and tracking SF files
    finput_SF_id = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2016/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root", "read");
    finput_SF_iso = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2016/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "read");
    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF_id, "SF")
      && getHistogram(h_SF_iso, finput_SF_iso, "SF")
      );
  }
  else if (theDataYear == 2017){
    // ID/Iso./IP and tracking SF files
    finput_SF_id = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2017/RunBCDEF_SF_ID.root", "read");
    finput_SF_iso = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2017/sf_mu_iso_susy_2017.root", "read");
    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF_id, "NUM_MediumPromptID_DEN_genTracks_pt_abseta")
      && getHistogram(h_SF_iso, finput_SF_iso, "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta")
      );
  }
  else if (theDataYear == 2018){
    // ID/Iso./IP and tracking SF files
    finput_SF_id = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2018/RunABCD_SF_ID.root", "read");
    finput_SF_iso = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Muons/2018/sf_mu_iso_susy_2017.root", "read");
    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF_id, "NUM_MediumPromptID_DEN_genTracks_pt_abseta")
      && getHistogram(h_SF_iso, finput_SF_iso, "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta")
      );
  }

  return res;
}
void MuonScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(finput_SF_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_iso);

  h_SF_id = nullptr;
  h_SF_iso = nullptr;
}

void MuonScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, MuonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const{
  if (!hist) return;
  if (!obj) return;

  float const pt = obj->pt();
  float const eta = obj->eta();

  int ix, iy;
  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();
  if (!etaOnY){
    ix = hist->GetXaxis()->FindBin((!useAbsEta ? eta : fabs(eta)));
    iy = hist->GetYaxis()->FindBin(pt);
  }
  else{
    ix = hist->GetXaxis()->FindBin(pt);
    iy = hist->GetYaxis()->FindBin((!useAbsEta ? eta : fabs(eta)));
  }
  if (ix==0) ix=1;
  else if (ix==nbinsx+1) ix=nbinsx;
  if (iy==0) iy=1;
  else if (iy==nbinsy+1) iy=nbinsy;

  float bc = hist->GetBinContent(ix, iy);
  float be = hist->GetBinError(ix, iy);
  if (bc!=0.f) be /= bc;
  if (be<0.f) be=0;

  theSF *= bc; theSFRelErr = sqrt(pow(theSFRelErr, 2)+pow(be, 2));
}

void MuonScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, bool useFastSim) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = ParticleSelectionHelpers::isTightParticle(obj);

  if (!useFastSim){
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_id, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_iso, true, true);
    }
  }
  else{
  }
}
