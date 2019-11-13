#include "CMS3/NtupleMaker/interface/METInfo.h"


METInfo::METInfo() :
  met(-1),
  metPhi(0),
  sumEt(-1),
  metSignificance(0),
  met_over_sqrtSumEt(-1),

  met_raw(-1),
  metPhi_raw(0),
  sumEt_raw(-1),

  met_JERUp(-1),
  metPhi_JERUp(0),
  met_JERDn(-1),
  metPhi_JERDn(0),

  met_JECUp(-1),
  metPhi_JECUp(0),
  met_JECDn(-1),
  metPhi_JECDn(0),

  met_MuonEnUp(-1),
  metPhi_MuonEnUp(0),
  met_MuonEnDn(-1),
  metPhi_MuonEnDn(0),

  met_ElectronEnUp(-1),
  metPhi_ElectronEnUp(0),
  met_ElectronEnDn(-1),
  metPhi_ElectronEnDn(0),

  met_TauEnUp(-1),
  metPhi_TauEnUp(0),
  met_TauEnDn(-1),
  metPhi_TauEnDn(0),

  met_UnclusteredEnUp(-1),
  metPhi_UnclusteredEnUp(0),
  met_UnclusteredEnDn(-1),
  metPhi_UnclusteredEnDn(0),

  met_PhotonEnUp(-1),
  metPhi_PhotonEnUp(0),
  met_PhotonEnDn(-1),
  metPhi_PhotonEnDn(0),

  calo_met(-1),
  calo_metPhi(0),

  gen_met(-1),
  gen_metPhi(0)
{}

METInfo::METInfo(const METInfo& other) :
  met(other.met),
  metPhi(other.metPhi),
  sumEt(other.sumEt),
  metSignificance(other.metSignificance),
  met_over_sqrtSumEt(other.met_over_sqrtSumEt),

  met_raw(other.met_raw),
  metPhi_raw(other.metPhi_raw),
  sumEt_raw(other.sumEt_raw),

  met_JERUp(other.met_JERUp),
  metPhi_JERUp(other.metPhi_JERUp),
  met_JERDn(other.met_JERDn),
  metPhi_JERDn(other.metPhi_JERDn),

  met_JECUp(other.met_JECUp),
  metPhi_JECUp(other.metPhi_JECUp),
  met_JECDn(other.met_JECDn),
  metPhi_JECDn(other.metPhi_JECDn),

  met_MuonEnUp(other.met_MuonEnUp),
  metPhi_MuonEnUp(other.metPhi_MuonEnUp),
  met_MuonEnDn(other.met_MuonEnDn),
  metPhi_MuonEnDn(other.metPhi_MuonEnDn),

  met_ElectronEnUp(other.met_ElectronEnUp),
  metPhi_ElectronEnUp(other.metPhi_ElectronEnUp),
  met_ElectronEnDn(other.met_ElectronEnDn),
  metPhi_ElectronEnDn(other.metPhi_ElectronEnDn),

  met_TauEnUp(other.met_TauEnUp),
  metPhi_TauEnUp(other.metPhi_TauEnUp),
  met_TauEnDn(other.met_TauEnDn),
  metPhi_TauEnDn(other.metPhi_TauEnDn),

  met_UnclusteredEnUp(other.met_UnclusteredEnUp),
  metPhi_UnclusteredEnUp(other.metPhi_UnclusteredEnUp),
  met_UnclusteredEnDn(other.met_UnclusteredEnDn),
  metPhi_UnclusteredEnDn(other.metPhi_UnclusteredEnDn),

  met_PhotonEnUp(other.met_PhotonEnUp),
  metPhi_PhotonEnUp(other.metPhi_PhotonEnUp),
  met_PhotonEnDn(other.met_PhotonEnDn),
  metPhi_PhotonEnDn(other.metPhi_PhotonEnDn),

  calo_met(other.calo_met),
  calo_metPhi(other.calo_metPhi),

  gen_met(other.gen_met),
  gen_metPhi(other.gen_metPhi)
{}
