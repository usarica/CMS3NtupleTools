import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDFilter("MuonMaker",
    # muon collection
    muonsInputTag = cms.InputTag("muons"),                         
    #track collection
    tracksInputTag = cms.InputTag("generalTracks"),
    #beamSpot (from CMS2)
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    #isoDeposits tags						 
    isoDepTrackInputTag = cms.InputTag("muIsoDepositTk"),
	isoDepEcalInputTag  = cms.InputTag("muIsoDepositCalByAssociatorTowers","ecal"),
	isoDepHcalInputTag  = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
)


