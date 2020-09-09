#include "MuonHandler.h"
#include "ElectronHandler.h"
#include "PhotonHandler.h"
#include "FSRHandler.h"
#include "ParticleDisambiguator.h"
#include "ParticleSelectionHelpers.h"



void ParticleDisambiguator::disambiguateParticles(
  std::vector<MuonObject*>*& muons,
  std::vector<ElectronObject*>*& electrons,
  std::vector<PhotonObject*>*& photons,
  bool doDeleteObjects
){
  // First disambiguate electrons from muons
  if (electrons){
    std::vector<ElectronObject*> electrons_new; electrons_new.reserve(electrons->size());
    for (auto*& product:(*electrons)){
      bool isTightProduct = ParticleSelectionHelpers::isTightParticle(product);
      bool isLooseProduct = ParticleSelectionHelpers::isLooseParticle(product);
      bool doRemove=false;
      if (!doRemove && muons){
        for (auto const* part:(*muons)){
          bool isTightPart = ParticleSelectionHelpers::isTightParticle(part);
          bool isLoosePart = ParticleSelectionHelpers::isLooseParticle(part);
          bool isVetoPart = ParticleSelectionHelpers::isVetoParticle(part);
          if (
            !(
              isTightPart
              ||
              (isLoosePart && !isTightProduct)
              ||
              (isVetoPart && !isLooseProduct)
              )
            ) continue;
          if (reco::deltaR(product->p4(), part->p4())<0.05){ doRemove=true; break; }
        }
      }
      if (!doRemove) electrons_new.push_back(product);
      else if (doDeleteObjects) delete product;
    }
    *electrons = electrons_new;
  }
  // Disambiguate photons from muons and *new* electrons
  if (photons){
    std::vector<PhotonObject*> photons_new; photons_new.reserve(photons->size());
    for (auto*& product:(*photons)){
      bool isTightProduct = ParticleSelectionHelpers::isTightParticle(product);
      bool isLooseProduct = ParticleSelectionHelpers::isLooseParticle(product);
      bool doRemove=false;
      if (!doRemove && muons){
        for (auto const* part:(*muons)){
          bool isTightPart = ParticleSelectionHelpers::isTightParticle(part);
          bool isLoosePart = ParticleSelectionHelpers::isLooseParticle(part);
          bool isVetoPart = ParticleSelectionHelpers::isVetoParticle(part);
          if (
            !(
              isTightPart
              ||
              (isLoosePart && !isTightProduct)
              ||
              (isVetoPart && !isLooseProduct)
              )
            ) continue;
          if (reco::deltaR(product->p4(), part->p4())<0.1){ doRemove=true; break; }
        }
      }
      if (!doRemove && electrons){
        for (auto const* part:(*electrons)){
          bool isTightPart = ParticleSelectionHelpers::isTightParticle(part);
          bool isLoosePart = ParticleSelectionHelpers::isLooseParticle(part);
          bool isVetoPart = ParticleSelectionHelpers::isVetoParticle(part);
          bool shareAllPFCands = (
            part->extras.n_associated_pfcands == product->extras.n_associated_pfcands
            &&
            part->extras.associated_pfcands_sum_px == product->extras.associated_pfcands_sum_px
            &&
            part->extras.associated_pfcands_sum_py == product->extras.associated_pfcands_sum_py
            );
          if (
            !(
              isTightPart
              ||
              (isLoosePart && !isTightProduct)
              ||
              (isVetoPart && !isLooseProduct)
              )
            ) continue;
          if (shareAllPFCands || reco::deltaR(product->p4(), part->p4())<0.1){ doRemove=true; break; }
        }
      }
      if (!doRemove) photons_new.push_back(product);
      else if (doDeleteObjects) delete product;
    }
    *photons = photons_new;
  }
}

void ParticleDisambiguator::disambiguateParticles(
  MuonHandler* muonHandle,
  ElectronHandler* electronHandle,
  PhotonHandler* photonHandle
){
  std::vector<MuonObject*>* muons = (muonHandle ? &(muonHandle->productList) : nullptr);
  std::vector<ElectronObject*>* electrons = (electronHandle ? &(electronHandle->productList) : nullptr);
  std::vector<PhotonObject*>* photons = (photonHandle ? &(photonHandle->productList) : nullptr);

  disambiguateParticles(muons, electrons, photons, true);
}


void ParticleDisambiguator::disambiguateParticles(FSRHandler* fsrHandle){
  if (!fsrHandle) return;

  std::vector<MuonObject*>* muons = &(fsrHandle->muons_postFSR);
  std::vector<ElectronObject*>* electrons = &(fsrHandle->electrons_postFSR);
  std::vector<PhotonObject*>* photons = &(fsrHandle->photons_postFSR);

  // Objects in the post-FSR containers will not be inaccessible if these containers drop them, so no need to delete the objects themselves.
  disambiguateParticles(muons, electrons, photons, false);
}
