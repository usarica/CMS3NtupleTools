#ifndef NTUPLEMAKER_PHOTONMAKER_H
#define NTUPLEMAKER_PHOTONMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <CommonTools/Utils/interface/StringCutObjectSelector.h>

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>


class PhotonMaker : public edm::stream::EDProducer<>{
public:
  explicit PhotonMaker(const edm::ParameterSet&);
  ~PhotonMaker();

protected:
  std::string aliasprefix_;
  int year_;

  edm::VParameterSet MVACuts_;

  std::unordered_map< std::string, std::vector< StringCutObjectSelector<pat::Photon, true> > > MVACutObjects;

  edm::EDGetTokenT< edm::View<pat::Photon> > photonsToken;

  edm::EDGetTokenT< edm::View<pat::PackedCandidate> > pfcandsToken;

  edm::EDGetTokenT< double > rhoToken;

  edm::EDGetTokenT< EcalRecHitCollection > ebhitsToken;
  edm::EDGetTokenT< EcalRecHitCollection > eehitsToken;


private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

  void setupMVACuts();
  void setMVAIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;
  void setCutBasedIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;

  void setCutBasedHGGIdSelectionBits(edm::View<pat::Photon>::const_iterator const&, pat::Photon&) const;
  void setEGammaPFPhotonIdSelectionBits(edm::View<pat::Photon>::const_iterator const&, pat::PackedCandidate const*, pat::Photon&) const;

  static float getRecHitEnergyTime(DetId const&, EcalRecHitCollection const*, EcalRecHitCollection const*, unsigned short, unsigned short, float* outtime=nullptr);

  void get_photon_pfphoton_matchMap(
    edm::Event const& iEvent,
    edm::Handle< edm::View<pat::Photon> > const& photonsHandle, edm::Handle< edm::View<pat::PackedCandidate> > const& pfcandsHandle,
    std::unordered_map<pat::Photon const*, pat::PackedCandidate const*>& res
  ) const;

};


#endif
