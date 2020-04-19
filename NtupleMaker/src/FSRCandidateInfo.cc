#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>


FSRCandidateInfo::FSRCandidateInfo() :
  obj(nullptr),
  fsrIso(0)
{}
FSRCandidateInfo::FSRCandidateInfo(FSRCandidateInfo const& other) :
  obj(other.obj),
  fsrIso(other.fsrIso),

  veto_electron_list(other.veto_electron_list),
  veto_photon_list(other.veto_photon_list),

  matched_muon_list(other.matched_muon_list),
  matched_electron_list(other.matched_electron_list)
{}

CMSLorentzVector_d FSRCandidateInfo::p4() const{ return obj->p4(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::pt() const{ return obj->pt(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::eta() const{ return obj->eta(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::phi() const{ return obj->phi(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::mass() const{ return obj->mass(); }
int FSRCandidateInfo::pdgId() const{ return obj->pdgId(); }
float FSRCandidateInfo::charge() const{ return obj->charge(); }
