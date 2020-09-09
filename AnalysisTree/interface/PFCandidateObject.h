#ifndef PFCANDIDATEOBJECT_H
#define PFCANDIDATEOBJECT_H

#include <string>

#include "ParticleObject.h"


#define PFCANDIDATE_VARIABLES \
PFCANDIDATE_VARIABLE(cms3_pfcand_qualityflag_t, qualityFlag, 0) \
PFCANDIDATE_VARIABLE(bool, is_associated_firstPV, 0) \
PFCANDIDATE_VARIABLE(float, dxy_associatedPV, 0) \
PFCANDIDATE_VARIABLE(float, dz_associatedPV, 0) \
PFCANDIDATE_VARIABLE(float, dxy_firstPV, 0) \
PFCANDIDATE_VARIABLE(float, dz_firstPV, 0) \
PFCANDIDATE_VARIABLE(cms3_listIndex_signed_short_t, matched_FSRCandidate_index, -1) \
PFCANDIDATE_VARIABLE(std::vector<cms3_listIndex_short_t>, matched_muon_index_list, std::vector<cms3_listIndex_short_t>()) \
PFCANDIDATE_VARIABLE(std::vector<cms3_listIndex_short_t>, matched_electron_index_list, std::vector<cms3_listIndex_short_t>()) \
PFCANDIDATE_VARIABLE(std::vector<cms3_listIndex_short_t>, matched_photon_index_list, std::vector<cms3_listIndex_short_t>()) \
PFCANDIDATE_VARIABLE(std::vector<cms3_listIndex_short_t>, matched_ak4jet_index_list, std::vector<cms3_listIndex_short_t>()) \
PFCANDIDATE_VARIABLE(std::vector<cms3_listIndex_short_t>, matched_ak8jet_index_list, std::vector<cms3_listIndex_short_t>())


class PFCandidateVariables{
public:
#define PFCANDIDATE_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  PFCANDIDATE_VARIABLES;
#undef PFCANDIDATE_VARIABLE

  PFCandidateVariables();
  PFCandidateVariables(PFCandidateVariables const& other);
  PFCandidateVariables& operator=(const PFCandidateVariables& other);

  void swap(PFCandidateVariables& other);

};

class PFCandidateObject : public ParticleObject{
private:
  // These enums come from DataFormats/PatCandidates/interface/PackedCandidate.h (latest update: CMSSW_10_2_X)
  enum qualityFlagsShiftsAndMasks {
    assignmentQualityMask = 0x7, assignmentQualityShift = 0,
    trackHighPurityMask  = 0x8, trackHighPurityShift=3,
    lostInnerHitsMask = 0x30, lostInnerHitsShift=4,
    muonFlagsMask = 0x0600, muonFlagsShift=9,
    egammaFlagsMask = 0x0800, egammaFlagsShift=11
  };

public:
  PFCandidateVariables extras;

  PFCandidateObject();
  PFCandidateObject(cms3_id_t id_);
  PFCandidateObject(cms3_id_t id_, LorentzVector_t const& mom_);
  PFCandidateObject(const PFCandidateObject& other);
  PFCandidateObject& operator=(const PFCandidateObject& other);
  ~PFCandidateObject();

  void swap(PFCandidateObject& other);

  bool isGlobalMuon() const;
  bool isStandaloneMuon() const;

};

#endif
