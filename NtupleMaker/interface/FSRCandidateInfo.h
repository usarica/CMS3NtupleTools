#ifndef CMS3_FSRCANDIDATEINFO_H
#define CMS3_FSRCANDIDATEINFO_H


#include <vector>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>

#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>


struct FSRCandidateInfo{
  pat::PackedCandidate const* obj;
  double fsrIso;

  std::vector<pat::Electron const*> veto_electron_list;
  std::vector<pat::Photon const*> veto_photon_list;

  std::vector<pat::Muon const*> matched_muon_list;
  std::vector<pat::Electron const*> matched_electron_list;

  FSRCandidateInfo();
  FSRCandidateInfo(FSRCandidateInfo const&);
  ~FSRCandidateInfo(){}

  CMSLorentzVector_d p4() const;
  CMSLorentzVector_d::Scalar pt() const;
  CMSLorentzVector_d::Scalar eta() const;
  CMSLorentzVector_d::Scalar phi() const;
  CMSLorentzVector_d::Scalar mass() const;
  int pdgId() const;
  float charge() const;

};

#endif
