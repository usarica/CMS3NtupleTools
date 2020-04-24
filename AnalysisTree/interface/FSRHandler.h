#ifndef FSRHANDLER_H
#define FSRHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "FSRObject.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "ParticleDisambiguator.h"


class FSRHandler : public IvyBase{
public:
  static const std::string colName;

protected:
  friend class ParticleDisambiguator;

  // Owned products
  std::vector<FSRObject*> fsrCandidates;
  std::vector<MuonObject*> muons_owned;
  std::vector<ElectronObject*> electrons_owned;

  // List of products to return, owned or not
  std::vector<MuonObject*> muons_postFSR;
  std::vector<ElectronObject*> electrons_postFSR;
  std::vector<PhotonObject*> photons_postFSR;

  void clear();

  bool constructFSRObjects();
  bool reconstructPostFSRObjects(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

public:
  // Constructors
  FSRHandler();

  // Destructors
  ~FSRHandler(){ clear(); }

  bool constructPostFSRParticles(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons);

  std::vector<FSRObject*> const& getFSRCandidates() const{ return fsrCandidates; }
  std::vector<MuonObject*> const& getMuons() const{ return muons_postFSR; }
  std::vector<ElectronObject*> const& getElectrons() const{ return electrons_postFSR; }
  std::vector<PhotonObject*> const& getPhotons() const{ return photons_postFSR; }

  static void bookBranches(BaseTree* tree);

};


#endif
