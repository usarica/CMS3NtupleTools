#ifndef PARTICLEDISAMBIGUATOR_H
#define PARTICLEDISAMBIGUATOR_H


class MuonHandler;
class ElectronHandler;
class PhotonHandler;


class ParticleDisambiguator{
public:
  ParticleDisambiguator(){};

  void disambiguateParticles(
    MuonHandler* muonHandle,
    ElectronHandler* electronHandle,
    PhotonHandler* photonHandle
  );

};


#endif
