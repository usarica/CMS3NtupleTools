#ifndef CMS3_PHOTONSELECTIONHELPERS_H
#define CMS3_PHOTONSELECTIONHELPERS_H

#include <iostream>
#include <cmath>
#include <unordered_map>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>

#include "MELAStreamHelpers.hh"


namespace reco{
  typedef std::vector<reco::VertexRef> VertexRefCollection;
}


namespace PhotonSelectionHelpers{
  enum EffectiveAreaType{
    PhotonEA_ch,
    PhotonEA_nh,
    PhotonEA_em
  };

  // Skim selection
  constexpr double selection_skim_pt = 20.;
  constexpr double selection_skim_eta = 5.0;
  constexpr double selection_iso_deltaR = 0.3; // This is a constant since the only isolation available uses dR<0.3

  float photonEffArea(pat::Photon const& obj, int const& year, PhotonSelectionHelpers::EffectiveAreaType const& eatype); // For PF and mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/photons_cff.py EAFile_MiniIso entries

  template<typename PFCandIterable> float photonPFIsoWorstChargedHadron(
    pat::Photon const& obj, int const& /*year*/, reco::VertexRefCollection const& vtxrefs,
    PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end,
    bool const& useLooseVtxFit, double const& thr_dz
  );

  float photonPFIsoChargedHadron(pat::Photon const& obj, int const& year, double const& rho, double const* iso_override=nullptr); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoNeutralHadron(pat::Photon const& obj, int const& year, double const& rho, double const* iso_override=nullptr); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoEM(pat::Photon const& obj, int const& year, double const& rho, double const* iso_override=nullptr); // Absolute PF iso. value, uses rho instead of delta beta
  float photonPFIsoComb(pat::Photon const& obj, int const& year, double const& rho, double const* iso_override_ch=nullptr, double const* iso_override_nh=nullptr, double const* iso_override_em=nullptr); // Absolute PF iso. value, uses rho instead of delta beta

  // Same as PFEGammaFilters::passPhotonSelection
  bool testEGammaPFPhotonSelection(pat::Photon const& obj, int const& year);
  // Same as PFEGammaFilters::isPhotonSafeForJetMET
  bool testEGammaPFPhotonMETSafetySelection(pat::Photon const& obj, pat::PackedCandidate const* pfCand, int const& year);

  bool testSkimPhoton(pat::Photon const& obj, int const& year);

}

template<typename PFCandIterable> float PhotonSelectionHelpers::photonPFIsoWorstChargedHadron(
  pat::Photon const& obj, int const& /*year*/, reco::VertexRefCollection const& vtxrefs,
  PFCandIterable const& pfcands_begin, PFCandIterable const& pfcands_end,
  bool const& useLooseVtxFit, double const& thr_dz
){
  constexpr double cut_deltaR = selection_iso_deltaR;
  constexpr double cut_deltaRself = 1e-5;
  constexpr double cut_pt = 0.;
  bool const addDZor = (thr_dz>0.);

  auto associated_pfcands_refvec = obj.associatedPackedPFCandidates();
  std::vector<pat::PackedCandidate const*> associated_pfcands; associated_pfcands.reserve(associated_pfcands_refvec.size());
  for (auto const& it_pfcand:associated_pfcands_refvec) associated_pfcands.push_back(&(*it_pfcand));

  double res = -1;
  std::unordered_map<int, double> vtx_iso_map;
  for (auto const& vtxref:vtxrefs) vtx_iso_map[vtxref.key()] = 0;
  for (PFCandIterable it_pfcands = pfcands_begin; it_pfcands!=pfcands_end; it_pfcands++){
    pat::PackedCandidate const* pfcand;
    CMS3ObjectHelpers::getObjectPointer(it_pfcands, pfcand);
    if (!pfcand) continue;

    unsigned int abs_id = std::abs(pfcand->pdgId());
    if (abs_id!=211) continue;

    // Omit candidates in veto cone
    if (HelperFunctions::checkListVariable(associated_pfcands, pfcand)) continue;

    double pt = pfcand->pt();
    if (pt<=cut_pt) continue;

    for (auto const& vtxref:vtxrefs){
      math::XYZVector phoWrtVtx = obj.superCluster()->position() - vtxref.get()->position();

      double dr = 0;
      HelperFunctions::deltaR<double>(phoWrtVtx.Eta(), phoWrtVtx.Phi(), pfcand->eta(), pfcand->phi(), dr);
      if (dr>=cut_deltaR || dr<cut_deltaRself) continue;

      int const vtxkey = vtxref.key();
      bool isFromPV = (
        //(!useLooseVtxFit && pfcand->fromPV(vtxkey)>=pat::PackedCandidate::PVTight)
        (!useLooseVtxFit && pfcand->vertexRef().key()==vtxkey && pfcand->pvAssociationQuality()>=pat::PackedCandidate::UsedInFitTight)
        ||
        (useLooseVtxFit && pfcand->vertexRef().key()==vtxkey && pfcand->pvAssociationQuality()>=pat::PackedCandidate::UsedInFitLoose)
        ||
        (addDZor && std::abs(pfcand->dz(vtxkey))<thr_dz)
        );
      if (!isFromPV) continue;

      auto it_tmp = vtx_iso_map.find(vtxkey);
      it_tmp->second += pt;
    }
  }

  auto it_found = vtx_iso_map.end();
  for (auto it = vtx_iso_map.begin(); it != vtx_iso_map.end(); it++){
    if (it_found==vtx_iso_map.end() || it->second>it_found->second) it_found = it;
  }
  if (it_found!=vtx_iso_map.end()) res = it_found->second;

  return res;
}


#endif
