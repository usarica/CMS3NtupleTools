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
    kEta,
    kMass,
    nCutTypes
  };
  enum TriggerObjectType{
    kMuon=0,
    kElectron,
    kPhoton,
    kAK4Jet,
    kAK8Jet,
    kHT,
    kMHT,
    kMET,
    nTriggerObjectTypes
  };

protected:
  TriggerObjectType type;

  float pt_cut;
  float eta_cut;
  float mass_cut;

  bool has_pt_cut;
  bool has_eta_cut;
  bool has_mass_cut;

public:
  HLTObjectProperties();
  HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_);
  HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_, std::vector< std::pair<HLTObjectProperties::CutType, float> > const& cuts_);
  HLTObjectProperties(HLTObjectProperties const& other);

  void setPtCut(float const& val){ pt_cut=val; has_pt_cut=(pt_cut>0.f); }
  void setEtaCut(float const& val){ eta_cut=val; has_eta_cut=(eta_cut>0.f); }
  void setMassCut(float const& val){ mass_cut=val; has_mass_cut=(mass_cut>0.f); }

  bool hasSameType(HLTObjectProperties::TriggerObjectType const& type_) const{ return (type==type_); }
  bool testCuts(ParticleObject::LorentzVector_t const& p4, HLTObjectProperties::TriggerObjectType const& type_) const;

  bool operator == (HLTObjectProperties const& other) const;
  bool operator != (HLTObjectProperties const& other) const{ return !(*this == other); }
  bool operator > (HLTObjectProperties const& other) const;
  bool operator >= (HLTObjectProperties const& other) const{ return (*this == other || *this > other); }

  static bool isMoreRestrictive(HLTObjectProperties const& earlier, HLTObjectProperties const& later){ return (earlier > later); }
  static bool isMoreRestrictive_ptr(HLTObjectProperties const* earlier, HLTObjectProperties const* later){ return (earlier && (!later || (*earlier > *later))); }

  static void sortByMoreRestrictive(std::vector<HLTObjectProperties>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive); }
  static void sortByMoreRestrictive(std::vector<HLTObjectProperties*>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive_ptr); }
  static void sortByMoreRestrictive(std::vector<HLTObjectProperties const*>& vec){ std::sort(vec.begin(), vec.end(), HLTObjectProperties::isMoreRestrictive_ptr); }

};


#endif
