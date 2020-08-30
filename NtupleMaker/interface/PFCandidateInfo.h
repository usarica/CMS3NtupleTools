#ifndef CMS3_PFCANDIDATEINFO_H
#define CMS3_PFCANDIDATEINFO_H


#include <vector>

#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>

#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>


struct PFCandidateInfo{
  typedef cms3_listIndex_short_t index_t;
  typedef std::vector<index_t> index_list_t;

  pat::PackedCandidate const* obj;

  cms3_listIndex_signed_short_t matchedFSRCandidate;

  index_list_t matched_muons;
  index_list_t matched_electrons;
  index_list_t matched_photons;

  index_list_t matched_ak4jets;
  index_list_t matched_ak8jets;

  PFCandidateInfo();
  PFCandidateInfo(pat::PackedCandidate const* obj_);
  PFCandidateInfo(PFCandidateInfo const&);
  ~PFCandidateInfo(){}

  CMSLorentzVector_d p4() const;
  CMSLorentzVector_d::Scalar pt() const;
  CMSLorentzVector_d::Scalar eta() const;
  CMSLorentzVector_d::Scalar phi() const;
  CMSLorentzVector_d::Scalar mass() const;
  cms3_id_t pdgId() const;
  float charge() const;

  index_list_t& getParticleMatchList(cms3_id_t const& id_);
  index_list_t const& getParticleMatchList(cms3_id_t const& id_) const;

  bool findParticleMatch(cms3_id_t const& id_, index_t const& idx) const;
  bool findAK4JetMatch(index_t const& idx) const;
  bool findAK8JetMatch(index_t const& idx) const;

  void addParticleMatch(cms3_id_t const& id_, index_t const& idx);
  void addAK4JetMatch(index_t const& idx);
  void addAK8JetMatch(index_t const& idx);

  void setFSRCandidateMatch(cms3_listIndex_short_t const& idx){ matchedFSRCandidate = idx; }

  unsigned short countRawParticleAssociations() const{ return matched_muons.size() + matched_electrons.size() + matched_photons.size(); }

  // Gives the actual count of particle overlaps
  void analyzeParticleOverlaps(
    std::vector<pat::Muon const*> const& muons, std::vector<pat::Electron const*> const& electrons, std::vector<pat::Photon const*> const& photons,
    cms3_listIndex_long_t& nImperfectOverlaps, cms3_listIndex_long_t& nPerfectOverlaps
  ) const;

  static PFCandidateInfo& make_and_get_PFCandidateInfo(std::vector<PFCandidateInfo>& pfcandInfos, pat::PackedCandidate const* pfcand);

  static void linkFSRCandidates(std::vector<FSRCandidateInfo>& fsrcandInfos, std::vector<PFCandidateInfo>& pfcandInfos);

};

#endif
