#include <CMS3/NtupleMaker/interface/FSRCandidateInfo.h>


CMSLorentzVector_d FSRCandidateInfo::p4() const{ return obj->p4(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::pt() const{ return obj->pt(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::eta() const{ return obj->eta(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::phi() const{ return obj->phi(); }
CMSLorentzVector_d::Scalar FSRCandidateInfo::mass() const{ return obj->mass(); }
int FSRCandidateInfo::pdgId() const{ return obj->pdgId(); }
float FSRCandidateInfo::charge() const{ return obj->charge(); }
