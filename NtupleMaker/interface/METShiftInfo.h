#ifndef METSHIFTINFO_H
#define METSHIFTINFO_H

#include <string>
#include <vector>

#include <DataFormats/Candidate/interface/Particle.h>


struct METShiftInfo{
  std::vector<reco::Particle::LorentzVector> metshifts;

  METShiftInfo();
  METShiftInfo(const METShiftInfo&);

};


#endif
