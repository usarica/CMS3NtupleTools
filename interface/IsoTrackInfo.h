#ifndef ISOTRACKINFO_H
#define ISOTRACKINFO_H

#include <string>
#include <vector>
#include <unordered_map>
#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>


struct IsoTrackInfo{
  CMSLorentzVector_d p4;

  int charge;
  int pdgId;

  float pfIso03_ch;
  float pfIso03_nh;
  float pfIso03_em;
  float pfIso03_db;
  float pfRelIso_comb;
  float miniIso_ch;
  float miniIso_nh;
  float miniIso_em;
  float miniIso_db;
  float miniRelIso_comb;

  bool fromPV;
  float dxy;
  float dz;
  float dxyerr;
  float dzerr;

  float deltaEta;
  float deltaPhi;

  bool is_pfCand;
  bool is_lostTrack;

  int nearestPFcand_id;
  float nearestPFcand_deltaR;
  CMSLorentzVector_d nearestPFcand_p4;

  float pterr;
  float normChi2;

  unsigned int tracker_hit_signature;
  int lastSubdet;
  int lastLayer;

  unsigned int n_layers_with_measurement;
  unsigned int n_layers_pixel_with_measurement;
  unsigned int n_valid_pixel_hits;
  unsigned int n_lost_pixel_inner_hits;
  unsigned int n_missing_inner_hits;
  unsigned int n_missing_outer_hits;

  float dEdxStrip;
  float dEdxPixel;

  float track_pt;
  float track_eta;
  float track_phi;

  bool is_highPurityTrack;
  bool is_tightTrack;

  float lepOverlap;
  float pfNeutralSum;

  float matchedCaloJetEmEnergy;
  float matchedCaloJetHadEnergy;

  //std::vector<uint16_t> crossedEcalStatus;
  //std::vector<uint32_t> crossedHcalStatus;

  IsoTrackInfo();
  IsoTrackInfo(IsoTrackInfo const&);

};


#endif