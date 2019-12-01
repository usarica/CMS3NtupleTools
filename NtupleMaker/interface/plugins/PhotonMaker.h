#ifndef NTUPLEMAKER_PHOTONMAKER_H
#define NTUPLEMAKER_PHOTONMAKER_H

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <CommonTools/Utils/interface/StringCutObjectSelector.h>

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"


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
  edm::EDGetTokenT< double > rhoToken;

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

  void setupMVACuts();
  void setMVAIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;
  void setCutBasedIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;

};


#endif