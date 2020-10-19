#ifndef CMS3_HLTOBJECTPROPERTIES_H
#define CMS3_HLTOBJECTPROPERTIES_H

#include <vector>
#include <utility>
#include <algorithm>
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include "ParticleObject.h"


class HLTObjectProperties{
public:
  enum CutType{
    kPt,
    kPtHigh,
    kEta,
    kEtaLow,
    kMass,
    nCutTypes
  };
  enum TriggerObjectType{
    kMuon=0,
    kElectron,
    kPhoton,
    kAK4Jet,
    kAK4DiJetSumWithDEtaDPhi,
    kAK8Jet,
    kHT,
    kHT_NoMu,
    kMET,
    kMET_NoMu,
    nTriggerObjectTypes
  };

protected:
  TriggerObjectType type;

  float pt_cut;
  float pt_high_cut;
  float eta_cut;
  float eta_low_cut;
  float mass_cut;

  bool has_pt_cut;
  bool has_pt_high_cut;
  bool has_eta_cut;
  bool has_eta_low_cut;
  bool has_mass_cut;

public:
  HLTObjectProperties();
  HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_);
  HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_, std::vector< std::pair<HLTObjectProperties::CutType, float> > const& cuts_);
  HLTObjectProperties(HLTObjectProperties const& other);

  void resetCuts();

  void setPtCut(float const& val){ pt_cut=val; has_pt_cut=(pt_cut>0.f); }
  void setPtLowCut(float const& val){ setPtCut(val); }
  void setPtHighCut(float const& val){ pt_high_cut=val; has_pt_high_cut=(pt_high_cut>0.f); }
  void setEtaCut(float const& val){ eta_cut=val; has_eta_cut=(eta_cut>0.f); }
  void setEtaLowCut(float const& val){ eta_low_cut=val; has_eta_low_cut=(eta_low_cut>0.f); }
  void setEtaHighCut(float const& val){ setEtaCut(val); }
  void setMassCut(float const& val){ mass_cut=val; has_mass_cut=(mass_cut>0.f); }

  TriggerObjectType const& getType() const{ return type; }
  bool hasSameType(HLTObjectProperties::TriggerObjectType const& type_) const{ return (type==type_); }
  bool testCuts(ParticleObject::LorentzVector_t const& p4, HLTObjectProperties::TriggerObjectType const& type_) const;

  bool operator == (HLTObjectProperties const& other) const;
  bool operator != (HLTObjectProperties const& other) const{ return !(*this == other); }
  bool operator > (HLTObjectProperties const& other) const;
  bool operator >= (HLTObjectProperties const& other) const{ return (*this == other || *this > other); }
  bool operator < (HLTObjectProperties const& other) const{ return (type == other.type && !(*this >= other)); }
  bool operator <= (HLTObjectProperties const& other) const{ return (*this == other || *this < other); }

  static bool isMoreRestrictive(HLTObjectProperties const& earlier, HLTObjectProperties const& later){ return (earlier > later); }
  static bool isMoreRestrictive_ptr(HLTObjectProperties const* earlier, HLTObjectProperties const* later){ return (earlier && (!later || (*earlier > *later))); }

  static void sortByMoreRestrictive(std::vector<HLTObjectProperties>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive); }
  static void sortByMoreRestrictive(std::vector<HLTObjectProperties*>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive_ptr); }
  static void sortByMoreRestrictive(std::vector<HLTObjectProperties const*>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive_ptr); }

};


#endif
