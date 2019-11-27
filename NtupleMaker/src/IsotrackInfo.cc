#include <CMS3/NtupleMaker/interface/IsotrackInfo.h>


IsotrackInfo::IsotrackInfo() :
  p4(0, 0, 0, 0),

  charge(0),
  id(0),

  pfIso03_ch(-1),
  pfIso03_nh(-1),
  pfIso03_em(-1),
  pfIso03_db(-1),
  pfIso03_comb_nofsr(-1),
  miniIso_ch(-1),
  miniIso_nh(-1),
  miniIso_em(-1),
  miniIso_db(-1),
  miniIso_comb_nofsr(-1),

  fromPV(false),
  dxy(0),
  dz(0),
  dxyerr(-1),
  dzerr(-1),

  deltaEta(0),
  deltaPhi(0),

  is_pfCand(false),
  is_lostTrack(false),

  nearestPFcand_id(0),
  nearestPFcand_deltaR(-1),
  nearestPFcand_p4(0, 0, 0, 0),

  pterr(-1),
  normChi2(-1),

  tracker_hit_signature(0),
  lastSubdet(-1),
  lastLayer(-1),

  n_layers_with_measurement(0),
  n_layers_pixel_with_measurement(0),
  n_valid_pixel_hits(0),
  n_lost_pixel_inner_hits(0),
  n_missing_inner_hits(0),
  n_missing_outer_hits(0),

  dEdxStrip(0),
  dEdxPixel(0),

  track_pt(0),
  track_eta(0),
  track_phi(0),

  is_highPurityTrack(0),
  is_tightTrack(0),

  lepOverlap(0),
  pfNeutralSum(0),

  matchedCaloJetEmEnergy(0),
  matchedCaloJetHadEnergy(0)
{}

IsotrackInfo::IsotrackInfo(IsotrackInfo const& other) : 
  p4(other.p4),

  charge(other.charge),
  id(other.id),

  pfIso03_ch(other.pfIso03_ch),
  pfIso03_nh(other.pfIso03_nh),
  pfIso03_em(other.pfIso03_em),
  pfIso03_db(other.pfIso03_db),
  pfIso03_comb_nofsr(other.pfIso03_comb_nofsr),
  miniIso_ch(other.miniIso_ch),
  miniIso_nh(other.miniIso_nh),
  miniIso_em(other.miniIso_em),
  miniIso_db(other.miniIso_db),
  miniIso_comb_nofsr(other.miniIso_comb_nofsr),

  fromPV(other.fromPV),
  dxy(other.dxy),
  dz(other.dz),
  dxyerr(other.dxyerr),
  dzerr(other.dzerr),

  deltaEta(other.deltaEta),
  deltaPhi(other.deltaPhi),

  is_pfCand(other.is_pfCand),
  is_lostTrack(other.is_lostTrack),

  nearestPFcand_id(other.nearestPFcand_id),
  nearestPFcand_deltaR(other.nearestPFcand_deltaR),
  nearestPFcand_p4(other.nearestPFcand_p4),

  pterr(other.pterr),
  normChi2(other.normChi2),

  tracker_hit_signature(other.tracker_hit_signature),
  lastSubdet(other.lastSubdet),
  lastLayer(other.lastLayer),

  n_layers_with_measurement(other.n_layers_with_measurement),
  n_layers_pixel_with_measurement(other.n_layers_pixel_with_measurement),
  n_valid_pixel_hits(other.n_valid_pixel_hits),
  n_lost_pixel_inner_hits(other.n_lost_pixel_inner_hits),
  n_missing_inner_hits(other.n_missing_inner_hits),
  n_missing_outer_hits(other.n_missing_outer_hits),

  dEdxStrip(other.dEdxStrip),
  dEdxPixel(other.dEdxPixel),

  track_pt(other.track_pt),
  track_eta(other.track_eta),
  track_phi(other.track_phi),

  is_highPurityTrack(other.is_highPurityTrack),
  is_tightTrack(other.is_tightTrack),

  lepOverlap(other.lepOverlap),
  pfNeutralSum(other.pfNeutralSum),

  matchedCaloJetEmEnergy(other.matchedCaloJetEmEnergy),
  matchedCaloJetHadEnergy(other.matchedCaloJetHadEnergy)

  //crossedEcalStatus(other.crossedEcalStatus),
  //crossedHcalStatus(other.crossedHcalStatus)
{}
