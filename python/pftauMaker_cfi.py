import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDProducer("PFTauMaker",
                            aliasPrefix = cms.untracked.string("taus_pf"),
                            # cms2PFJetsTag = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
                            # referencePFJetsTag = cms.InputTag("ak5PFJets"),
                            # particleFlowTag = cms.InputTag("particleFlow"),

                            # PFTau collection
                            # pftausInputTag = cms.InputTag("hpsPFTauProducer"),
                            pftausInputTag = cms.InputTag("slimmedTaus"),

                            # add discriminators to the list below
                            # details: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID
                            # check available discriminators: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py
                            tauIDCollection                     = cms.untracked.vstring(
        "againstElectronDeadECAL",
        "againstElectronLoose",
        "againstElectronMedium",
        "againstElectronTight",
        "againstMuonLoose",     
        "againstMuonLoose2",    
        "againstMuonLoose3",    
        "againstMuonMedium",    
        "againstMuonMedium2",   
        "againstMuonTight",     
        "againstMuonTight2",    
        "againstMuonTight3",   
        "byCombinedIsolationDeltaBetaCorrRaw",
        "byCombinedIsolationDeltaBetaCorrRaw3Hits",
        "byLooseCombinedIsolationDeltaBetaCorr",  
        "byLooseCombinedIsolationDeltaBetaCorr3Hits",
        "byMediumCombinedIsolationDeltaBetaCorr", 
        "byMediumCombinedIsolationDeltaBetaCorr3Hits", 
        "byTightCombinedIsolationDeltaBetaCorr",  
        "byTightCombinedIsolationDeltaBetaCorr3Hits", 
        "byVLooseCombinedIsolationDeltaBetaCorr", 
        "byDecayModeFinding",
        ),

                            # ** remove this later after things work ** #
        # discriminator
        # available in miniAOD
        # againstElectronDeadECAL                     = cms.InputTag("hpsPFTauDiscriminationByDeadECALElectronRejection"),
        # againstElectronLoose                   = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"),
        # againstElectronMedium                  = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
        # againstElectronTight                   = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection"),
        # againstMuonLoose                       = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),
        # againstMuonLoose2                      = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection2"),
        # againstMuonLoose3                      = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),
        # againstMuonMedium                      = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection"),
        # againstMuonMedium2                     = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection2"),
        # againstMuonTight                       = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection"),
        # againstMuonTight2                      = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection2"),
        # againstMuonTight3                      = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
        # byCombinedIsolationDeltaBetaCorrRaw    = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr"),
        # byCombinedIsolationDeltaBetaCorrRaw3Hits    = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
        # byLooseCombinedIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
        # byLooseCombinedIsolationDeltaBetaCorr3Hits  = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
        # byMediumCombinedIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"),
        # byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
        # byTightCombinedIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"),
        # byTightCombinedIsolationDeltaBetaCorr3Hits  = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
        # byVLooseCombinedIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"),
        # byDecayModeFinding                     = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),

        # # not available in miniAOD
        # byIsolationMVAraw                      = cms.InputTag("hpsPFTauDiscriminationByIsolationMVAraw"),
        # byLooseIsolationMVA                    = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA"),
        # byMediumIsolationMVA                   = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA"),
        # byTightIsolationMVA                    = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA"),
        # byIsolationMVA2raw                     = cms.InputTag("hpsPFTauDiscriminationByIsolationMVA2raw"),
        # byLooseIsolationMVA2                   = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA2"),
        # byMediumIsolationMVA2                  = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA2"),
        # byTightIsolationMVA2                   = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA2"),
        # againstElectronMVA                     = cms.InputTag("hpsPFTauDiscriminationByMVAElectronRejection"),
        # againstElectronMVA2raw                 = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection"),
        # againstElectronMVA2category            = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection:category"),
        # againstElectronVLooseMVA2              = cms.InputTag("hpsPFTauDiscriminationByMVA2VLooseElectronRejection"),
        # againstElectronLooseMVA2               = cms.InputTag("hpsPFTauDiscriminationByMVA2LooseElectronRejection"),
        # againstElectronMediumMVA2              = cms.InputTag("hpsPFTauDiscriminationByMVA2MediumElectronRejection"),
        # againstElectronTightMVA2               = cms.InputTag("hpsPFTauDiscriminationByMVA2TightElectronRejection"),
        # againstElectronMVA3raw                      = cms.InputTag("hpsPFTauDiscriminationByMVA3rawElectronRejection"),
        # againstElectronMVA3category                 = cms.InputTag("hpsPFTauDiscriminationByMVA3rawElectronRejection:category"),
        # againstElectronLooseMVA3                    = cms.InputTag("hpsPFTauDiscriminationByMVA3LooseElectronRejection"),
        # againstElectronMediumMVA3                   = cms.InputTag("hpsPFTauDiscriminationByMVA3MediumElectronRejection"),
        # againstElectronTightMVA3                    = cms.InputTag("hpsPFTauDiscriminationByMVA3TightElectronRejection"),
        # againstElectronVTightMVA3                   = cms.InputTag("hpsPFTauDiscriminationByMVA3VTightElectronRejection"),

)




