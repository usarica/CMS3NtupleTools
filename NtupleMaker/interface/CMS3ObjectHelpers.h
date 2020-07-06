#ifndef CMS3OBJECTHELPERS_H
#define CMS3OBJECTHELPERS_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>
#include <DataFormats/JetReco/interface/GenJet.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>


namespace CMS3ObjectHelpers{
  enum ObjectType{
    kMuon = 0,
    kElectron,
    kPhoton,
    kAK4Jet,
    kAK8Jet,

    kGenPrunedParticle,
    kGenPackedParticle,

    kGenAK4Jet,
    kGenAK8Jet,

    nObjectTypes
  };

  enum ObjectMatchType{
    kMatchBy_DeltaR,
    kMatchBy_DeltaR_PDGid,
    kMatchBy_FourMomentum,
    kMatchBy_FourMomentum_PDGid
  };

  // edm::View<obj_type> holds const obj_type references, whereas std::vector<const obj_type*> holds pointers...
  template<typename T> void getObjectPointer(typename edm::View<T>::const_iterator const& it, T const*& ptr);
  template<typename T> void getObjectPointer(typename std::vector<T>::const_iterator const& it, T const*& ptr);
  template<typename T> void getObjectPointer(typename std::vector<T const*>::const_iterator const& it, T const*& ptr);

  template<typename T> bool objHasGreaterPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterPt(T const* earlier, T const* later);

  template<typename T> bool objHasGreaterPz(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterPz(T const* earlier, T const* later);

  template<typename T> void sortByGreaterPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterPt(std::vector<T*>& vec);

  template<typename T> void sortByGreaterPz(std::vector<T>& vec);
  template<typename T> void sortByGreaterPz(std::vector<T*>& vec);

  // T and U are pointer types
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles(
    CMS3ObjectHelpers::ObjectMatchType type,
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map<T, U>& key_val_map
  );
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles(
    CMS3ObjectHelpers::ObjectMatchType type, double match_thr,
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map<T, U>& key_val_map
  );
  // The values in the map below are ordered in the smallest matching criterion
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles_OneToMany(
    CMS3ObjectHelpers::ObjectMatchType type, double match_thr,
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map< T, typename std::vector<U> >& key_val_map
  );

  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles_OneToMany(
    bool(*comparator)(T, U),
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map< T, typename std::vector<U> >& key_val_map
  );

}

template<typename T> void CMS3ObjectHelpers::getObjectPointer(typename edm::View<T>::const_iterator const& it, T const*& ptr){ ptr = &(*it); }
template<typename T> void CMS3ObjectHelpers::getObjectPointer(typename std::vector<T>::const_iterator const& it, T const*& ptr){ ptr = &(*it); }
template<typename T> void CMS3ObjectHelpers::getObjectPointer(typename std::vector<T const*>::const_iterator const& it, T const*& ptr){ ptr = *it; }

template<typename T> bool CMS3ObjectHelpers::objHasGreaterPt(T const& earlier, T const& later){ return (earlier.pt() > later.pt()); }
template<typename T> bool CMS3ObjectHelpers::ptrHasGreaterPt(T const* earlier, T const* later){ return (earlier && (!later || (earlier->pt() > later->pt()))); }

template<typename T> bool CMS3ObjectHelpers::objHasGreaterPz(T const& earlier, T const& later){ return (earlier.pz() > later.pz()); }
template<typename T> bool CMS3ObjectHelpers::ptrHasGreaterPz(T const* earlier, T const* later){ return (earlier && (!later || (earlier->pz() > later->pz()))); }

template<typename T> void CMS3ObjectHelpers::sortByGreaterPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), CMS3ObjectHelpers::objHasGreaterPt<T>); }
template<typename T> void CMS3ObjectHelpers::sortByGreaterPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), CMS3ObjectHelpers::ptrHasGreaterPt<T>); }

template<typename T> void CMS3ObjectHelpers::sortByGreaterPz(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), CMS3ObjectHelpers::objHasGreaterPz<T>); }
template<typename T> void CMS3ObjectHelpers::sortByGreaterPz(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), CMS3ObjectHelpers::ptrHasGreaterPz<T>); }

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void CMS3ObjectHelpers::matchParticles(
  CMS3ObjectHelpers::ObjectMatchType type,
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map<T, U>& key_val_map
){
  bool matchByDeltaR = (type==kMatchBy_DeltaR || type==kMatchBy_DeltaR_PDGid);
  bool matchByFourMomentum = (type==kMatchBy_FourMomentum || type==kMatchBy_FourMomentum_PDGid);
  bool matchByPDGid = (type==kMatchBy_DeltaR_PDGid || type==kMatchBy_FourMomentum_PDGid);

  std::vector<T> remaining_keys; remaining_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_keys.push_back(ptr);
  }
  std::vector<U> remaining_vals; remaining_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_vals.push_back(ptr);
  }

  while (!remaining_keys.empty() && !remaining_vals.empty()){
    auto chosenKey = remaining_keys.end();
    auto chosenVal = remaining_vals.end();
    double minDeltaR=-1;
    double minDeltaEucDot=-1;
    for (auto it_key = remaining_keys.begin(); it_key != remaining_keys.end(); it_key++){
      T key;
      CMS3ObjectHelpers::getObjectPointer(it_key, key);
      auto const pKey = key->p4();

      for (auto it_val = remaining_vals.begin(); it_val != remaining_vals.end(); it_val++){
        U val;
        CMS3ObjectHelpers::getObjectPointer(it_val, val);

        if (matchByPDGid && key->pdgId()!=val->pdgId()) continue;
        auto const pVal = val->p4();

        double deltaR = std::abs(reco::deltaR(pVal, pKey));
        double euc_dot_prod = pKey.px()*pVal.px() + pKey.py()*pVal.py() + pKey.pz()*pVal.pz() + pKey.energy()*pVal.energy();
        double euc_dot_prod_ref = pKey.px()*pKey.px() + pKey.py()*pKey.py() + pKey.pz()*pKey.pz() + pKey.energy()*pKey.energy();
        double deltaEucDot = std::abs(euc_dot_prod - euc_dot_prod_ref);

        if (
          (matchByDeltaR && (minDeltaR==-1. || deltaR<minDeltaR))
          ||
          (matchByFourMomentum && (minDeltaEucDot==-1. || deltaEucDot<minDeltaEucDot))
          ){
          if (matchByDeltaR) minDeltaR=deltaR;
          else if (matchByFourMomentum) minDeltaEucDot=deltaEucDot;
          chosenKey=it_key;
          chosenVal=it_val;
        }
      }
    }

    if (chosenKey!=remaining_keys.end() && chosenVal!=remaining_vals.end()){
      T key; U val;
      CMS3ObjectHelpers::getObjectPointer(chosenKey, key);
      CMS3ObjectHelpers::getObjectPointer(chosenVal, val);
      key_val_map[key] = val;
    }
    for (auto it=remaining_keys.begin(); it!=remaining_keys.end(); it++){ if (it == chosenKey){ remaining_keys.erase(it); break; } }
    for (auto it=remaining_vals.begin(); it!=remaining_vals.end(); it++){ if (it == chosenVal){ remaining_vals.erase(it); break; } }
  }
}

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void CMS3ObjectHelpers::matchParticles(
  CMS3ObjectHelpers::ObjectMatchType type, double match_thr,
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map<T, U>& key_val_map
){
  bool matchByDeltaR = (type==kMatchBy_DeltaR || type==kMatchBy_DeltaR_PDGid);
  bool matchByFourMomentum = (type==kMatchBy_FourMomentum || type==kMatchBy_FourMomentum_PDGid);
  bool matchByPDGid = (type==kMatchBy_DeltaR_PDGid || type==kMatchBy_FourMomentum_PDGid);

  std::vector<T> remaining_keys; remaining_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_keys.push_back(ptr);
  }
  std::vector<U> remaining_vals; remaining_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_vals.push_back(ptr);
  }

  while (!remaining_keys.empty() && !remaining_vals.empty()){
    auto chosenKey = remaining_keys.end();
    auto chosenVal = remaining_vals.end();
    double minDeltaR=-1;
    double minDeltaEucDot=-1;
    for (auto it_key = remaining_keys.begin(); it_key != remaining_keys.end(); it_key++){
      T key;
      CMS3ObjectHelpers::getObjectPointer(it_key, key);
      auto const pKey = key->p4();

      for (auto it_val = remaining_vals.begin(); it_val != remaining_vals.end(); it_val++){
        U val;
        CMS3ObjectHelpers::getObjectPointer(it_val, val);

        if (matchByPDGid && key->pdgId()!=val->pdgId()) continue;
        auto const pVal = val->p4();

        double deltaR = std::abs(reco::deltaR(pVal, pKey));
        double euc_dot_prod = pKey.px()*pVal.px() + pKey.py()*pVal.py() + pKey.pz()*pVal.pz() + pKey.energy()*pVal.energy();
        double euc_dot_prod_ref = pKey.px()*pKey.px() + pKey.py()*pKey.py() + pKey.pz()*pKey.pz() + pKey.energy()*pKey.energy();
        double deltaEucDot = std::abs(euc_dot_prod - euc_dot_prod_ref);

        if (
          (matchByDeltaR && (minDeltaR==-1. || deltaR<minDeltaR))
          ||
          (matchByFourMomentum && (minDeltaEucDot==-1. || deltaEucDot<minDeltaEucDot))
          ){
          if (matchByDeltaR) minDeltaR=deltaR;
          else if (matchByFourMomentum) minDeltaEucDot=deltaEucDot;
          chosenKey=it_key;
          chosenVal=it_val;
        }
      }
    }

    if (chosenKey!=remaining_keys.end() && chosenVal!=remaining_vals.end()){
      T key; U val;
      CMS3ObjectHelpers::getObjectPointer(chosenKey, key);
      CMS3ObjectHelpers::getObjectPointer(chosenVal, val);
      key_val_map[key] = val;
      if (match_thr>0.){
        if (matchByDeltaR){
          double const deltaR = std::abs(reco::deltaR(val->p4(), key->p4()));
          if (deltaR>=match_thr){
            key_val_map[key] = nullptr;
            chosenVal = remaining_vals.end();
          }
        }
      }
    }
    for (auto it=remaining_keys.begin(); it!=remaining_keys.end(); it++){ if (it == chosenKey){ remaining_keys.erase(it); break; } }
    for (auto it=remaining_vals.begin(); it!=remaining_vals.end(); it++){ if (it == chosenVal){ remaining_vals.erase(it); break; } }
  }
}

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void CMS3ObjectHelpers::matchParticles_OneToMany(
  CMS3ObjectHelpers::ObjectMatchType type, double match_thr,
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map< T, typename std::vector<U> >& key_val_map
){
  bool matchByDeltaR = (type==kMatchBy_DeltaR || type==kMatchBy_DeltaR_PDGid);
  bool matchByPDGid = (type==kMatchBy_DeltaR_PDGid);

  std::vector<T> all_keys; all_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_keys.push_back(ptr);
  }

  std::vector<U> all_vals; all_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_vals.push_back(ptr);
  }

  {
    std::vector<U> empty_vals; if (!all_vals.empty()) empty_vals.reserve(all_vals.size());
    for (auto const& key:all_keys) key_val_map[key] = empty_vals;
  }
  if (!matchByDeltaR) return;

  for (auto it_key = all_keys.begin(); it_key != all_keys.end(); it_key++){
    T key;
    CMS3ObjectHelpers::getObjectPointer(it_key, key);
    auto const pKey = key->p4();
    auto map_iterator = key_val_map.find(key);
    auto& target_val_list = map_iterator->second;

    typename std::vector<std::pair<double, U>> criterion_val_pair_list;
    for (auto it_val = all_vals.begin(); it_val != all_vals.end(); it_val++){
      U val;
      CMS3ObjectHelpers::getObjectPointer(it_val, val);

      if (matchByPDGid && key->pdgId()!=val->pdgId()) continue;
      auto const pVal = val->p4();

      double deltaR = std::abs(reco::deltaR(pVal, pKey));
      if (matchByDeltaR && (match_thr<0. || deltaR<match_thr)) HelperFunctions::addByLowest(criterion_val_pair_list, deltaR, val);
    }

    for (auto const& pp:criterion_val_pair_list) target_val_list.push_back(pp.second);
  }
}

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void CMS3ObjectHelpers::matchParticles_OneToMany(
  bool(*comparator)(T, U),
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map< T, typename std::vector<U> >& key_val_map
){
  std::vector<T> all_keys; all_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_keys.push_back(ptr);
  }

  std::vector<U> all_vals; all_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    CMS3ObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_vals.push_back(ptr);
  }

  {
    std::vector<U> empty_vals; if (!all_vals.empty()) empty_vals.reserve(all_vals.size());
    for (auto const& key:all_keys) key_val_map[key] = empty_vals;
  }

  for (auto it_key = all_keys.begin(); it_key != all_keys.end(); it_key++){
    T key;
    CMS3ObjectHelpers::getObjectPointer(it_key, key);
    auto map_iterator = key_val_map.find(key);
    auto& target_val_list = map_iterator->second;

    for (auto it_val = all_vals.begin(); it_val != all_vals.end(); it_val++){
      U val;
      CMS3ObjectHelpers::getObjectPointer(it_val, val);

      if (comparator(key, val)) target_val_list.push_back(val);
    }

  }
}


#endif
