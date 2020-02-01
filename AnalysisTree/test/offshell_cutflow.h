#include "common_includes.h"
#include "ParticleSelectionHelpers.h"


using namespace std;


bool check_dPhi_pTll_pTmiss(float const& val, size_t const& n_ak4jets_tight, bool use_old=true){
  float thr;
  if (use_old) thr=0.5;
  else{
    if (n_ak4jets_tight==0) thr=2.6;
    else if (n_ak4jets_tight==1) thr=1.6;
    else /*if (n_ak4jets_tight==2)*/ thr=1.0;
  }
  return std::abs(val)>=thr;
}

bool check_dPhi_pTlljets_pTmiss(float const& val, bool use_old=true){
  if (use_old) return true;
  float thr=2.6;
  return std::abs(val)>=thr;
}

bool check_min_abs_dPhi_pTj_pTmiss(float const& val, bool use_old=true){
  float thr;
  if (use_old) thr=0.5;
  else thr=0.5/*0.6*/;
  return std::abs(val)>=thr;
}

#define _old_ZZ_met_thr_ 125.f
#define _new_ZZ_met_thr_ 125.f
#define _old_WW_met_thr_ 20.f
#define _new_WW_met_thr_ 20.f
bool check_pTmiss_over_pTlljets(float const& pTmiss, float const& pTlljets, bool use_old=true, bool useZZ=true){
  return true;
  if (use_old || pTlljets==0.f) return true;
  float thr=std::pow((useZZ ? _new_ZZ_met_thr_ : _new_WW_met_thr_)/pTlljets, 1.5);
  return pTmiss/pTlljets>=thr;
}

bool check_pTmiss(float const& val, bool use_old=true, bool useZZ=true){
  return val>=(use_old ? (useZZ ? _old_ZZ_met_thr_ : _old_WW_met_thr_) : (useZZ ? _new_ZZ_met_thr_ : _new_WW_met_thr_));
}

bool check_Nb_veto(size_t const& nbs){
  return (nbs==0);
}

bool check_pTl1(float const& val, bool use_old=true){
  return val>=(use_old ? 25.f : 25.f);
}

bool check_pTl2(float const& val, bool use_old=true){
  return val>=(use_old ? 25.f : 25.f);
}

bool check_pTll(float const& val, bool use_old=true){
  return val>=(use_old ? 55.f : 55.f);
}

bool check_ml1(float const& val, bool use_old=true, bool useZZ=true){
  if (useZZ) return val>=(use_old ? 76.f : 81.f) && val<(use_old ? 106.f : 101.f);
  else return val>=(use_old ? 106.f : 101.f);
}

bool check_VBF_category(float const& kd, std::vector<AK4JetObject*> const& ak4jets_tight, ParticleObject const* theChosenCand, bool use_old=true){
  if (!use_old) return kd>=0.8;
  else{
    if (ak4jets_tight.size()<2) return false;
    auto itFirstJet = ak4jets_tight.cbegin();
    AK4JetObject* leading_jet = *itFirstJet; itFirstJet++;
    AK4JetObject* subleading_jet = *itFirstJet; itFirstJet++;
    float leading_eta = leading_jet->eta();
    float subleading_eta = subleading_jet->eta();
    if (leading_eta<subleading_eta) std::swap(leading_eta, subleading_eta);
    if (std::abs(leading_eta - subleading_eta)<4.f) return false;
    if ((leading_jet->p4() + subleading_jet->p4()).M()<500.f) return false;
    float eta_cand = theChosenCand->eta();
    if (eta_cand<=subleading_eta || eta_cand>=leading_eta) return false;
    for (auto it=itFirstJet; it!=ak4jets_tight.cend(); it++){
      float eta_jet = (*it)->eta();
      if (eta_jet>subleading_eta && eta_jet<leading_eta) return false;
    }
    return true;
  }
}

bool check_NoExtraLeptons(std::vector<MuonObject*> const& muons, std::vector<ElectronObject*> const& electrons, ParticleObject const* theChosenCand, bool use_old=true){
  typedef bool (*MuonVetoFcn)(MuonObject const*);
  typedef bool (*ElectronVetoFcn)(ElectronObject const*);

  MuonVetoFcn muonVetoFcn=nullptr;
  ElectronVetoFcn electronVetoFcn=nullptr;
  //if (use_old){
    muonVetoFcn = &ParticleSelectionHelpers::isLooseParticle<MuonObject>;
    electronVetoFcn = &ParticleSelectionHelpers::isLooseParticle<ElectronObject>;
  //}
  //else{
  //  muonVetoFcn = &ParticleSelectionHelpers::isVetoParticle<MuonObject>;
  //  electronVetoFcn = &ParticleSelectionHelpers::isVetoParticle<ElectronObject>;
  //}

  std::vector<ParticleObject*> const& daughters = theChosenCand->getDaughters();
  size_t nleptons = muons.size()+electrons.size();
  std::vector<ParticleObject*> extraVetoLeptons; if (nleptons>0) extraVetoLeptons.reserve(nleptons);

  for (auto const& muon:muons){
    if (ParticleObject::checkParticleExists(dynamic_cast<ParticleObject*>(muon), daughters)) continue;
    if (muonVetoFcn(muon)) extraVetoLeptons.push_back(muon);
  }
  for (auto const& electron:electrons){
    if (ParticleObject::checkParticleExists(dynamic_cast<ParticleObject*>(electron), daughters)) continue;
    if (electronVetoFcn(electron)) extraVetoLeptons.push_back(electron);
  }

  return extraVetoLeptons.empty();
}
