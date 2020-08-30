#include <cmath>

#include <FWCore/Utilities/interface/EDMException.h>

#include <CMS3/NtupleMaker/interface/PFCandidateInfo.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>


PFCandidateInfo::PFCandidateInfo() :
  obj(nullptr),
  matchedFSRCandidate(-1)
{}
PFCandidateInfo::PFCandidateInfo(pat::PackedCandidate const* obj_) :
  obj(obj_),
  matchedFSRCandidate(-1)
{}
PFCandidateInfo::PFCandidateInfo(PFCandidateInfo const& other) :
  obj(other.obj),

  matchedFSRCandidate(other.matchedFSRCandidate),

  matched_muons(other.matched_muons),
  matched_electrons(other.matched_electrons),
  matched_photons(other.matched_photons),

  matched_ak4jets(other.matched_ak4jets),
  matched_ak8jets(other.matched_ak8jets)
{}

CMSLorentzVector_d PFCandidateInfo::p4() const{ return obj->p4(); }
CMSLorentzVector_d::Scalar PFCandidateInfo::pt() const{ return obj->pt(); }
CMSLorentzVector_d::Scalar PFCandidateInfo::eta() const{ return obj->eta(); }
CMSLorentzVector_d::Scalar PFCandidateInfo::phi() const{ return obj->phi(); }
CMSLorentzVector_d::Scalar PFCandidateInfo::mass() const{ return obj->mass(); }
cms3_id_t PFCandidateInfo::pdgId() const{ return obj->pdgId(); }
float PFCandidateInfo::charge() const{ return obj->charge(); }

PFCandidateInfo::index_list_t& PFCandidateInfo::getParticleMatchList(cms3_id_t const& id_){
  cms3_absid_t abs_id = std::abs(id_);
  index_list_t* idx_list = nullptr;
  switch (abs_id){
  case 13:
    idx_list = &(this->matched_muons);
    break;
  case 11:
    idx_list = &(this->matched_electrons);
    break;
  case 22:
    idx_list = &(this->matched_photons);
    break;
  default:
    throw cms::Exception("PFCandidateInfo") << "PFCandidateInfo::getParticleMatchList: Particle id " << id_ << " is not defined.";
    break;
  }
  return *idx_list;
}
PFCandidateInfo::index_list_t const& PFCandidateInfo::getParticleMatchList(cms3_id_t const& id_) const{
  cms3_absid_t abs_id = std::abs(id_);
  index_list_t const* idx_list = nullptr;
  switch (abs_id){
  case 13:
    idx_list = &(this->matched_muons);
    break;
  case 11:
    idx_list = &(this->matched_electrons);
    break;
  case 22:
    idx_list = &(this->matched_photons);
    break;
  default:
    throw cms::Exception("PFCandidateInfo") << "PFCandidateInfo::getParticleMatchList: Particle id " << id_ << " is not defined.";
    break;
  }
  return *idx_list;
}

bool PFCandidateInfo::findParticleMatch(cms3_id_t const& id_, PFCandidateInfo::index_t const& idx) const{
  index_list_t const& idx_list = getParticleMatchList(id_);
  return HelperFunctions::checkListVariable(idx_list, idx);
}
bool PFCandidateInfo::findAK4JetMatch(PFCandidateInfo::index_t const& idx) const{
  return HelperFunctions::checkListVariable(matched_ak4jets, idx);
}
bool PFCandidateInfo::findAK8JetMatch(PFCandidateInfo::index_t const& idx) const{
  return HelperFunctions::checkListVariable(matched_ak8jets, idx);
}

void PFCandidateInfo::addParticleMatch(cms3_id_t const& id_, PFCandidateInfo::index_t const& idx){
  if (!findParticleMatch(id_, idx)){
    index_list_t& idx_list = getParticleMatchList(id_);
    idx_list.push_back(idx);
  }
}
void PFCandidateInfo::addAK4JetMatch(PFCandidateInfo::index_t const& idx){
  if (!findAK4JetMatch(idx)) matched_ak4jets.push_back(idx);
}
void PFCandidateInfo::addAK8JetMatch(PFCandidateInfo::index_t const& idx){
  if (!findAK8JetMatch(idx)) matched_ak8jets.push_back(idx);
}

void PFCandidateInfo::analyzeParticleOverlaps(
  std::vector<pat::Muon const*> const& /*muons*/, std::vector<pat::Electron const*> const& electrons, std::vector<pat::Photon const*> const& photons,
  cms3_listIndex_long_t& nImperfectOverlaps, cms3_listIndex_long_t& nPerfectOverlaps
) const{
  // Count total possible matches
  cms3_listIndex_long_t n_matches_muon = matched_muons.size();
  cms3_listIndex_long_t n_matches_electron = matched_electrons.size();
  cms3_listIndex_long_t n_matches_photon = matched_photons.size();
  nImperfectOverlaps = (
    (n_matches_muon*(n_matches_muon-1))/2 + n_matches_muon*(n_matches_electron + n_matches_photon)
    +
    (n_matches_electron*(n_matches_electron-1))/2 + n_matches_electron*n_matches_photon
    +
    (n_matches_photon*(n_matches_photon-1))/2
    );

  std::vector< std::pair<index_t, index_t> > eg_perfect_matches; eg_perfect_matches.reserve(matched_electrons.size()*matched_photons.size());
  for (auto const& iele:matched_electrons){
    pat::Electron const* electron = electrons.at(iele);
    cms3_listSize_t electron_n_associated_pfcands = electron->userInt("n_associated_pfcands");
    float electron_associated_pfcands_sum_sc_pt = electron->userFloat("associated_pfcands_sum_sc_pt");
    for (auto const& ipho:matched_photons){
      pat::Photon const* photon = photons.at(ipho);
      cms3_listSize_t photon_n_associated_pfcands = photon->userInt("n_associated_pfcands");
      float photon_associated_pfcands_sum_sc_pt = photon->userFloat("associated_pfcands_sum_sc_pt");
      if (
        electron_n_associated_pfcands == photon_n_associated_pfcands
        &&
        std::abs(electron_associated_pfcands_sum_sc_pt/photon_associated_pfcands_sum_sc_pt-1.f)<1e-3
        ) eg_perfect_matches.emplace_back(iele, ipho);
    }
  }
  nPerfectOverlaps = eg_perfect_matches.size();
  nImperfectOverlaps -= nPerfectOverlaps;
}

PFCandidateInfo& PFCandidateInfo::make_and_get_PFCandidateInfo(std::vector<PFCandidateInfo>& pfcandInfos, pat::PackedCandidate const* pfcand){
  if (!pfcand) throw cms::Exception("PFCandidateInfo") << "PFCandidateInfo::make_and_get_PFCandidateInfo: PF candidate is null.";
  PFCandidateInfo* res = nullptr;
  for (auto& pfcandInfo:pfcandInfos){
    if (pfcandInfo.obj == pfcand){
      res = &pfcandInfo;
      break;
    }
  }
  if (!res){
    pfcandInfos.emplace_back(pfcand);
    res = &(pfcandInfos.back());
  }
  return *res;
}

void PFCandidateInfo::linkFSRCandidates(std::vector<FSRCandidateInfo>& fsrcandInfos, std::vector<PFCandidateInfo>& pfcandInfos){
  for (auto& pfcandInfo:pfcandInfos){
    unsigned int ifsr = 0;
    for (auto const& fsrcandInfo:fsrcandInfos){
      if (pfcandInfo.obj == fsrcandInfo.obj){
        pfcandInfo.setFSRCandidateMatch(ifsr);
        break;
      }
      ifsr++;
    }
  }
}
