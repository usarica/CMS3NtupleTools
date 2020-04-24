#ifndef PARTICLEDISAMBIGUATOR_H
#define PARTICLEDISAMBIGUATOR_H

#include <vector>
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"


class MuonHandler;
class ElectronHandler;
class PhotonHandler;
class FSRHandler;


class ParticleDisambiguator{
protected:
  void disambiguateParticles(
    std::vector<MuonObject*>*& muons,
    std::vector<ElectronObject*>*& electrons,
    std::vector<PhotonObject*>*& photons,
    bool doDeleteObjects
  );

public:
  ParticleDisambiguator(){};

  void disambiguateParticles(
    MuonHandler* muonHandle,
    ElectronHandler* electronHandle,
    PhotonHandler* photonHandle
  );

  void disambiguateParticles(FSRHandler* fsrHandle);

};


#endif
