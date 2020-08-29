#ifndef CMS3_PFCANDIDATEINFO_H
#define CMS3_PFCANDIDATEINFO_H


#include <vector>

#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>

#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>

#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>


struct PFCandidateInfo{
  pat::PackedCandidate const* obj;

  short matchedFSRCandidate;

  std::vector<unsigned int> matched_muons;
  std::vector<unsigned int> matched_electrons;
  std::vector<unsigned int> matched_photons;

  std::vector<unsigned int> matched_ak4jets;
  std::vector<unsigned int> matched_ak8jets;

  PFCandidateInfo();
  PFCandidateInfo(pat::PackedCandidate const* obj_);
  PFCandidateInfo(PFCandidateInfo const&);
  ~PFCandidateInfo(){}

  CMSLorentzVector_d p4() const;
  CMSLorentzVector_d::Scalar pt() const;
  CMSLorentzVector_d::Scalar eta() const;
  CMSLorentzVector_d::Scalar phi() const;
  CMSLorentzVector_d::Scalar mass() const;
  int pdgId() const;
  float charge() const;

  std::vector<unsigned int>& getParticleMatchList(int const& id_);
  std::vector<unsigned int> const& getParticleMatchList(int const& id_) const;

  bool findParticleMatch(int const& id_, unsigned int const& idx) const;
  bool findAK4JetMatch(unsigned int const& idx) const;
  bool findAK8JetMatch(unsigned int const& idx) const;

  void addParticleMatch(int const& id_, unsigned int const& idx);
  void addAK4JetMatch(unsigned int const& idx);
  void addAK8JetMatch(unsigned int const& idx);

  void setFSRCandidateMatch(unsigned int const& idx){ matchedFSRCandidate = idx; }

  unsigned short countRawParticleAssociations() const{ return matched_muons.size() + matched_electrons.size() + matched_photons.size(); }

  // Gives the actual count of particle overlaps
  void analyzeParticleOverlaps(
    std::vector<pat::Muon const*> const& muons, std::vector<pat::Electron const*> const& electrons, std::vector<pat::Photon const*> const& photons,
    unsigned short& nImperfectOverlaps, unsigned short& nPerfectOverlaps
  ) const;

  static PFCandidateInfo& make_and_get_PFCandidateInfo(std::vector<PFCandidateInfo>& pfcandInfos, pat::PackedCandidate const* pfcand);

  static void linkFSRCandidates(std::vector<FSRCandidateInfo>& fsrcandInfos, std::vector<PFCandidateInfo>& pfcandInfos);

};

#endif
