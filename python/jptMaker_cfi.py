import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.ZSPJetCorrections219_cff import ZSPJetCorrectorIcone5

ZSPJetCorrectorSisCone5 = ZSPJetCorrectorIcone5.clone()

ZSPJetCorrectorSisCone5.label = cms.string('ZSPJetCorrectorSisCone5')

ZSPJetCorJetIcone5 = cms.EDProducer("CaloJetCorrectionProducer"                              ,
                                    src        = cms.InputTag("sisCone5CaloJets"            ),
                                    correctors = cms.vstring('ZSPJetCorrectorSisCone5'      ),
                                    alias      = cms.untracked.string('ZSPJetCorJetSisCone5')
)

from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import JetPlusTrackZSPCorrectorIcone5
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import ZSPiterativeCone5JetTracksAssociatorAtVertex
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import ZSPiterativeCone5JetTracksAssociatorAtCaloFace
from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import ZSPiterativeCone5JetExtender

JetPlusTrackZSPCorrectorSisCone5 = JetPlusTrackZSPCorrectorIcone5.clone()

JetPlusTrackZSPCorrectorSisCone5.JetTrackCollectionAtCalo   = cms.InputTag("ZSPsisCone5JetTracksAssociatorAtCaloFace")
JetPlusTrackZSPCorrectorSisCone5.JetTrackCollectionAtVertex = cms.InputTag("ZSPsisCone5JetTracksAssociatorAtVertex"  )
JetPlusTrackZSPCorrectorSisCone5.label                      = cms.string('JetPlusTrackZSPCorrectorSisCone5'          )

JetPlusTrackZSPCorJetSisCone5 = cms.EDProducer("CaloJetCorrectionProducer"                                       ,
                                               src        = cms.InputTag("ZSPJetCorJetSisCone5"                 ),
                                               correctors = cms.vstring('JetPlusTrackZSPCorrectorSisCone5'      ),
                                               alias      = cms.untracked.string('JetPlusTrackZSPCorJetSisCone5')
)

ZSPSisCone5JetTracksAssociatorAtVertex      = ZSPiterativeCone5JetTracksAssociatorAtVertex.clone()
ZSPSisCone5JetTracksAssociatorAtVertex.jets = cms.InputTag("ZSPJetCorJetSiscone5")

ZSPSisCone5JetTracksAssociatorAtCaloFace      = ZSPiterativeCone5JetTracksAssociatorAtCaloFace.clone()
ZSPSisCone5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ZSPJetCorJetSiscone5")

ZSPSisCone5JetExtender                  = ZSPiterativeCone5JetExtender.clone()
ZSPSisCone5JetExtender.jets             = cms.InputTag("ZSPJetCorJetSiscone5"                    )
ZSPSisCone5JetExtender.jet2TracksAtCALO = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtCaloFace")
ZSPSisCone5JetExtender.jet2TracksAtVX   = cms.InputTag("ZSPSisCone5JetTracksAssociatorAtVertex"  )

ZSPrecoJetAssociations = cms.Sequence(ZSPSisCone5JetTracksAssociatorAtVertex*ZSPSisCone5JetTracksAssociatorAtCaloFace*ZSPSisCone5JetExtender)

JetPlusTrackCorrections = cms.Sequence(ZSPrecoJetAssociations*JetPlusTrackZSPCorJetSisCone5)

jptMaker = cms.EDFilter("JPTMaker",
                        jptInputTag      = cms.InputTag("JetPlusTrackZSPCorJetSisCone5"),
                        caloJetsInputTag = cms.InputTag("sisCone5CaloJets"             )
)

JPTCorrections = cms.Sequence( JetPlusTrackCorrections * jptMaker )
