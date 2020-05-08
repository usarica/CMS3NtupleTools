#ifndef PARTICLEOBJECTHELPERS_H
#define PARTICLEOBJECTHELPERS_H

#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include "ParticleObject.h"
#include "HelperFunctions.h"
#include "TLorentzVector.h"


namespace ParticleObjectHelpers{
  enum ObjectMatchType{
    kMatchBy_DeltaR,
    kMatchBy_DeltaR_PDGid,
    kMatchBy_FourMomentum,
    kMatchBy_FourMomentum_PDGid
  };

  template<typename T> bool objHasGreaterPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterPt(T const* earlier, T const* later);

  template<typename T> bool objHasGreaterScalarSumPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterScalarSumPt(T const* earlier, T const* later);

  template<typename T> bool objHasGreaterScalarSumPt_ImmediateDaughters(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterScalarSumPt_ImmediateDaughters(T const* earlier, T const* later);

  template<typename T> void sortByGreaterPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterPt(std::vector<T*>& vec);

  template<typename T> void sortByGreaterScalarSumPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterScalarSumPt(std::vector<T*>& vec);

  template<typename T> void sortByGreaterScalarSumPt_ImmediateDaughters(std::vector<T>& vec);
  template<typename T> void sortByGreaterScalarSumPt_ImmediateDaughters(std::vector<T*>& vec);

  template<typename T> TLorentzVector convertCMSLorentzVectorToTLorentzVector(T const& p4);

  // Matching funtions
  template<typename T> void getObjectPointer(typename std::vector<T>::const_iterator const& it, T const*& ptr);
  template<typename T> void getObjectPointer(typename std::vector<T const*>::const_iterator const& it, T const*& ptr);

  // T and U are pointer types
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles(
    ParticleObjectHelpers::ObjectMatchType type,
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map<T, U>& key_val_map
  );
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles(
    ParticleObjectHelpers::ObjectMatchType type, double match_thr,
    Iterable_T const& keys_begin, Iterable_T const& keys_end,
    Iterable_U const& vals_begin, Iterable_U const& vals_end,
    typename std::unordered_map<T, U>& key_val_map
  );
  // The values in the map below are ordered in the smallest matching criterion
  template<typename T, typename U, typename Iterable_T, typename Iterable_U> void matchParticles_OneToMany(
    ParticleObjectHelpers::ObjectMatchType type, double match_thr,
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

template<typename T> bool ParticleObjectHelpers::objHasGreaterPt(T const& earlier, T const& later){ return (earlier.pt() > later.pt()); }
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterPt(T const* earlier, T const* later){ return (earlier && (!later || (earlier->pt() > later->pt()))); }

template<typename T> bool ParticleObjectHelpers::objHasGreaterScalarSumPt(T const& earlier, T const& later){
  std::vector<ParticleObject const*> deepDaus_earlier; earlier.getDeepDaughters(deepDaus_earlier);
  std::vector<ParticleObject const*> deepDaus_later; later.getDeepDaughters(deepDaus_later);
  float scsumpt_earlier=0; for (auto const& dau:deepDaus_earlier) scsumpt_earlier += dau->pt();
  float scsumpt_later=0; for (auto const& dau:deepDaus_later) scsumpt_later += dau->pt();
  return (scsumpt_earlier > scsumpt_later);
}
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterScalarSumPt(T const* earlier, T const* later){
  return (earlier && ((later && ParticleObjectHelpers::objHasGreaterScalarSumPt(*earlier, *later)) || !later));
}

template<typename T> bool ParticleObjectHelpers::objHasGreaterScalarSumPt_ImmediateDaughters(T const& earlier, T const& later){
  std::vector<ParticleObject*> const& shallowDaus_earlier = earlier.getDaughters();
  std::vector<ParticleObject*> const& shallowDaus_later = later.getDaughters();
  float scsumpt_earlier=0; for (auto const& dau:shallowDaus_earlier) scsumpt_earlier += dau->pt();
  float scsumpt_later=0; for (auto const& dau:shallowDaus_later) scsumpt_later += dau->pt();
  return (scsumpt_earlier > scsumpt_later);
}
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterScalarSumPt_ImmediateDaughters(T const* earlier, T const* later){
  return (earlier && ((later && ParticleObjectHelpers::objHasGreaterScalarSumPt_ImmediateDaughters(*earlier, *later)) || !later));
}

template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterPt<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterPt<T>); }

template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterScalarSumPt<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterScalarSumPt<T>); }

template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt_ImmediateDaughters(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterScalarSumPt_ImmediateDaughters<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt_ImmediateDaughters(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterScalarSumPt_ImmediateDaughters<T>); }

template<typename T> TLorentzVector ParticleObjectHelpers::convertCMSLorentzVectorToTLorentzVector(T const& p4){ return TLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.T()); }

template<typename T> void ParticleObjectHelpers::getObjectPointer(typename std::vector<T>::const_iterator const& it, T const*& ptr){ ptr = &(*it); }
template<typename T> void ParticleObjectHelpers::getObjectPointer(typename std::vector<T const*>::const_iterator const& it, T const*& ptr){ ptr = *it; }

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void ParticleObjectHelpers::matchParticles(
  ParticleObjectHelpers::ObjectMatchType type,
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
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_keys.push_back(ptr);
  }
  std::vector<U> remaining_vals; remaining_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_vals.push_back(ptr);
  }

  while (!remaining_keys.empty() && !remaining_vals.empty()){
    auto chosenKey = remaining_keys.end();
    auto chosenVal = remaining_vals.end();
    double minDeltaR=-1;
    double minDeltaEucDot=-1;
    for (auto it_key = remaining_keys.begin(); it_key != remaining_keys.end(); it_key++){
      T key;
      ParticleObjectHelpers::getObjectPointer(it_key, key);
      auto const pKey = key->p4();

      for (auto it_val = remaining_vals.begin(); it_val != remaining_vals.end(); it_val++){
        U val;
        ParticleObjectHelpers::getObjectPointer(it_val, val);

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
      ParticleObjectHelpers::getObjectPointer(chosenKey, key);
      ParticleObjectHelpers::getObjectPointer(chosenVal, val);
      key_val_map[key] = val;
    }
    for (auto it=remaining_keys.begin(); it!=remaining_keys.end(); it++){ if (it == chosenKey){ remaining_keys.erase(it); break; } }
    for (auto it=remaining_vals.begin(); it!=remaining_vals.end(); it++){ if (it == chosenVal){ remaining_vals.erase(it); break; } }
  }
}

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void ParticleObjectHelpers::matchParticles(
  ParticleObjectHelpers::ObjectMatchType type, double match_thr,
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
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_keys.push_back(ptr);
  }
  std::vector<U> remaining_vals; remaining_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) remaining_vals.push_back(ptr);
  }

  while (!remaining_keys.empty() && !remaining_vals.empty()){
    auto chosenKey = remaining_keys.end();
    auto chosenVal = remaining_vals.end();
    double minDeltaR=-1;
    double minDeltaEucDot=-1;
    for (auto it_key = remaining_keys.begin(); it_key != remaining_keys.end(); it_key++){
      T key;
      ParticleObjectHelpers::getObjectPointer(it_key, key);
      auto const pKey = key->p4();

      for (auto it_val = remaining_vals.begin(); it_val != remaining_vals.end(); it_val++){
        U val;
        ParticleObjectHelpers::getObjectPointer(it_val, val);

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
      ParticleObjectHelpers::getObjectPointer(chosenKey, key);
      ParticleObjectHelpers::getObjectPointer(chosenVal, val);
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

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void ParticleObjectHelpers::matchParticles_OneToMany(
  ParticleObjectHelpers::ObjectMatchType type, double match_thr,
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map< T, typename std::vector<U> >& key_val_map
){
  bool matchByDeltaR = (type==kMatchBy_DeltaR || type==kMatchBy_DeltaR_PDGid);
  bool matchByPDGid = (type==kMatchBy_DeltaR_PDGid);

  std::vector<T> all_keys; all_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_keys.push_back(ptr);
  }

  std::vector<U> all_vals; all_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_vals.push_back(ptr);
  }

  {
    std::vector<U> empty_vals; if (!all_vals.empty()) empty_vals.reserve(all_vals.size());
    for (auto const& key:all_keys) key_val_map[key] = empty_vals;
  }
  if (!matchByDeltaR) return;

  for (auto it_key = all_keys.begin(); it_key != all_keys.end(); it_key++){
    T key;
    ParticleObjectHelpers::getObjectPointer(it_key, key);
    auto const pKey = key->p4();
    auto map_iterator = key_val_map.find(key);
    auto& target_val_list = map_iterator->second;

    typename std::vector<std::pair<double, U>> criterion_val_pair_list;
    for (auto it_val = all_vals.begin(); it_val != all_vals.end(); it_val++){
      U val;
      ParticleObjectHelpers::getObjectPointer(it_val, val);

      if (matchByPDGid && key->pdgId()!=val->pdgId()) continue;
      auto const pVal = val->p4();

      double deltaR = std::abs(reco::deltaR(pVal, pKey));
      if (matchByDeltaR && (match_thr<0. || deltaR<match_thr)) HelperFunctions::addByLowest(criterion_val_pair_list, deltaR, val);
    }

    for (auto const& pp:criterion_val_pair_list) target_val_list.push_back(pp.second);
  }
}

template<typename T, typename U, typename Iterable_T, typename Iterable_U> void ParticleObjectHelpers::matchParticles_OneToMany(
  bool(*comparator)(T, U),
  Iterable_T const& keys_begin, Iterable_T const& keys_end,
  Iterable_U const& vals_begin, Iterable_U const& vals_end,
  typename std::unordered_map< T, typename std::vector<U> >& key_val_map
){
  std::vector<T> all_keys; all_keys.reserve(static_cast<size_t>(keys_end-keys_begin));
  for (Iterable_T it = keys_begin; it != keys_end; it++){
    T ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_keys.push_back(ptr);
  }

  std::vector<U> all_vals; all_vals.reserve(static_cast<size_t>(vals_end-vals_begin));
  for (Iterable_U it = vals_begin; it != vals_end; it++){
    U ptr;
    ParticleObjectHelpers::getObjectPointer(it, ptr);
    if (ptr) all_vals.push_back(ptr);
  }

  {
    std::vector<U> empty_vals; if (!all_vals.empty()) empty_vals.reserve(all_vals.size());
    for (auto const& key:all_keys) key_val_map[key] = empty_vals;
  }

  for (auto it_key = all_keys.begin(); it_key != all_keys.end(); it_key++){
    T key;
    ParticleObjectHelpers::getObjectPointer(it_key, key);
    auto map_iterator = key_val_map.find(key);
    auto& target_val_list = map_iterator->second;

    for (auto it_val = all_vals.begin(); it_val != all_vals.end(); it_val++){
      U val;
      ParticleObjectHelpers::getObjectPointer(it_val, val);

      if (comparator(key, val)) target_val_list.push_back(val);
    }

  }
}


#endif
