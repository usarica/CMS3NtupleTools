#include <cassert>
#include "DileptonObject.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


namespace OffshellCutflow{
  constexpr float MZ_VAL_CUTS = 91.2f;

  enum FinalStateType{
    fs_ZZ_2l2nu,
    fs_WW_2l2nu,
    fs_ZW_3l1nu,
    nFinalStateTypes
  };

  // Final state configuration
  FinalStateType activeFinalState = nFinalStateTypes;

  void setActiveFinalState(FinalStateType const& type){ activeFinalState = type; }

  // Selection functions
  bool check_dPhi_pTll_pTmiss(float const& val){
    constexpr float thr = 1.f;
    return std::abs(val)>=thr;
  }

  bool check_dPhi_pTlljets_pTmiss(float const& val){
    constexpr float thr=2.5f;
    return std::abs(val)>=thr;
  }

  bool check_min_abs_dPhi_pTj_pTmiss(float const& val, int njets=-1){
    const float thr = (njets<2 || activeFinalState==fs_ZW_3l1nu ? 0.25f : 0.5f);
    return std::abs(val)>=thr;
  }

  bool check_pTmiss(float const& val, int njets=-1){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
      thr=(njets<2 ? 125.f : 140.f);
      break;
    case fs_WW_2l2nu:
    case fs_ZW_3l1nu:
      thr=20.f;
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTmiss: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }

  bool check_Nb_veto(size_t const& nbs){
    return (nbs==0);
  }

  bool check_pTl1(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
    case fs_ZW_3l1nu:
      thr=25.f;
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTl1: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }

  bool check_pTl2(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
      thr=25.f;
      break;
    case fs_ZW_3l1nu:
      thr=20.f;
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTl2: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }

  bool check_pTl3(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
      IVYerr << "OffshellCutflow::check_pTl3: This function should not be called for the active final state " << activeFinalState << "." << endl;
      return false;
      break;
    case fs_ZW_3l1nu:
      thr=10.f;
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTl3: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }

  bool check_pTZ1(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZW_3l1nu:
      thr=30.f;
      break;
    case fs_ZZ_2l2nu:
      return check_pTl1(val);
    case fs_WW_2l2nu:
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTZ1: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }
  bool check_pTZ2(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZW_3l1nu:
      thr=20.f;
      break;
    case fs_ZZ_2l2nu:
      return check_pTl2(val);
    case fs_WW_2l2nu:
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTZ2: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }
  bool check_pTlW(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZW_3l1nu:
      thr=20.f;
      break;
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTlW: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }
  bool check_mTlW_pTmiss(cms3_id_t const& id_l3, float const& mTlW, float const& pTmiss){
    cms3_absid_t abs_id = std::abs(id_l3);
    const double A = (abs_id==11 ? 90. : 75.);
    const double B = 120.;
    const double thr_mTlW = (abs_id==11 ? 10. : 20.);
    switch (activeFinalState){
    case fs_ZW_3l1nu:
      return (mTlW/A + pTmiss/B)>=1. && mTlW>=thr_mTlW;
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
      return true;
    default:
      IVYerr << "OffshellCutflow::check_mTlW_pTmiss: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      return true;
    }
  }

  // Please don't rename this function to check_pTll... It becomes painful to distinguish from the pTl1 version (and multiple bugs appeared in the past because of this).
  bool check_pTboson(float const& val){
    float thr=-1;
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
      thr=55.f;
      break;
    case fs_WW_2l2nu:
    case fs_ZW_3l1nu:
      break;
    default:
      IVYerr << "OffshellCutflow::check_pTboson: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return val>=thr;
  }

  bool check_mll_QCDsuppression(float const& val){
    bool res = true;
    switch (activeFinalState){
    case fs_ZW_3l1nu:
      res = (val>4.f);
      break;
    case fs_ZZ_2l2nu:
    case fs_WW_2l2nu:
      break;
    default:
      IVYerr << "OffshellCutflow::check_mll_QCDsuppression: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return res;
  }

  bool check_mll(float const& val, bool const& isSF){
    float thr[2]={ -1, -1 };
    switch (activeFinalState){
    case fs_ZZ_2l2nu:
    case fs_ZW_3l1nu:
      thr[0]=MZ_VAL_CUTS-15.f;
      thr[1]=MZ_VAL_CUTS+15.f;
      break;
    case fs_WW_2l2nu:
      thr[0]=MZ_VAL_CUTS+15.f;
      break;
    default:
      IVYerr << "OffshellCutflow::check_mll: The active final state " << activeFinalState << " is not specified." << endl;
      assert(0);
      break;
    }
    return (thr[0]<0.f || val>=thr[0]) && (thr[1]<0.f || val<thr[1]);
  }

  bool check_mll(DileptonObject const& dilepton){
    return check_mll(dilepton.m(), dilepton.isSF());
  }

  bool check_cutbased_VBF_category(std::vector<AK4JetObject*> const& ak4jets_tight, ParticleObject const* theChosenCand){
    if (!theChosenCand) return false;

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
