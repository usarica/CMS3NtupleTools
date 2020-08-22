#ifndef CMS3_PHOTONSELECTIONHELPERS_H
#define CMS3_PHOTONSELECTIONHELPERS_H

#include <cmath>
#include <unordered_map>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>


namespace PhotonSelectionHelpers{
  enum EffectiveAreaType{
    PhotonEA_ch,
    PhotonEA_nh,
    PhotonEA_em
  };

  // Skim selection
  constexpr double selection_skim_pt = 5.;
  constexpr double selection_skim_eta = 5.0;
  constexpr double selection_iso_deltaR = 0.3; // This is a constant since the only isolation available uses dR<0.3

  // FSR skim
  constexpr double selection_match_fsr_deltaR = 0.5;
  constexpr double selection_skim_fsr_pt = 2.;
  constexpr double selection_skim_fsr_min_reldr = 0.012; // min. DR/pt**2 threshold
  constexpr double selection_skim_fsr_reliso = 1.8; // fsrIso / pT threshold

  float photonEffArea(pat::Photon const& obj, int const& year, PhotonSelectionHelpers::EffectiveAreaType const& eatype); // For PF and mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/photons_cff.py EAFile_MiniIso entries

  template<typename PFCandIterable> float photonPFIsoWorstChargedHadron(pat::Photon const& obj, int const& year, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end);

  float photonPFIsoChargedHadron(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoNeutralHadron(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoEM(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoComb(pat::Photon const& obj, int const& year, double const& rho); // Absolute PF iso. value, uses rho instead of delta beta

  // Same as PFEGammaFilters::passPhotonSelection
  bool testEGammaPFPhotonSelection(pat::Photon const& obj, int const& year);
  // Same as PFEGammaFilters::isPhotonSafeForJetMET
  bool testEGammaPFPhotonMETSafetySelection(pat::Photon const& obj, pat::PackedCandidate const* pfCand, int const& year);

  bool testSkimPhoton(pat::Photon const& obj, int const& year);

}

template<typename PFCandIterable> float PhotonSelectionHelpers::photonPFIsoWorstChargedHadron(pat::Photon const& obj, int const& /*year*/, PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end){
  constexpr double cut_deltaR = selection_iso_deltaR;

  constexpr double cut_deltaRself = 1e-5;
  constexpr double cut_pt = 0.;

  auto associated_pfcands_refvec = obj.associatedPackedPFCandidates();
  std::vector<pat::PackedCandidate const*> associated_pfcands; associated_pfcands.reserve(associated_pfcands_refvec.size());
  for (auto const& it_pfcand:associated_pfcands_refvec) associated_pfcands.push_back(&(*it_pfcand));

  double res = 0;
  std::unordered_map<reco::Vertex const*, double> vtx_iso_map;
  for (PFCandIterable it_pfcands = pfcands_begin; it_pfcands!=pfcands_end; it_pfcands++){
    pat::PackedCandidate const* pfcand;
    CMS3ObjectHelpers::getObjectPointer(it_pfcands, pfcand);
    if (!pfcand) continue;

    // Omit candidates in veto cone
    if (HelperFunctions::checkListVariable(associated_pfcands, pfcand)) continue;

    unsigned int abs_id = std::abs(pfcand->pdgId());
    if (abs_id!=211) continue;

    int charge = pfcand->charge();
    if (charge==0) continue;

    double pt = pfcand->pt();
    if (pt<=cut_pt) continue;

    double dr = reco::deltaR(obj.p4(), pfcand->p4());
    if (dr>=cut_deltaR || dr<cut_deltaRself) continue;

    reco::VertexRef vtxref = pfcand->vertexRef();
    if (!vtxref.isNonnull()) continue;

    if (pfcand->fromPV(vtxref.key()) < pat::PackedCandidate::PVTight) continue;

    reco::Vertex const* vtx = vtxref.get();
    auto it_tmp = vtx_iso_map.find(vtx);
    if (it_tmp == vtx_iso_map.end()){
      vtx_iso_map[vtx] = 0;
      it_tmp = vtx_iso_map.find(vtx);
    }
    it_tmp->second += pt;
  }

  for (auto& it:vtx_iso_map) res = std::max(res, it.second);
  return res;
}


#endif
