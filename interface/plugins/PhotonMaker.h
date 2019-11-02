// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PhotonMaker.h,v 1.6 2010/04/24 03:58:59 warren Exp $
//
//
#ifndef NTUPLEMAKER_PHOTONMAKER_H
#define NTUPLEMAKER_PHOTONMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

  edm::EDGetTokenT< edm::View<pat::Photon> > photonsToken;

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void produce(edm::Event&, const edm::EventSetup&);

  void setMVAIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;
  void setCutBasedIdUserVariables(edm::View<pat::Photon>::const_iterator const&, pat::Photon&, std::string const&, std::string const&) const;

};


#endif
