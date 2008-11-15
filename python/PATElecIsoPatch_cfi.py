import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *

# input pat sequences
from PhysicsTools.PatAlgos.patLayer0_cff import *
from PhysicsTools.PatAlgos.patLayer1_cff import  *


eleIsoDepositHcalFromTowersDepth1 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet( EleIsoHcalFromTowersExtractorBlock )
)
eleIsoDepositHcalFromTowersDepth1.ExtractorPSet.hcalDepth = cms.int32(1)

eleIsoDepositHcalFromTowersDepth2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet( EleIsoHcalFromTowersExtractorBlock )
)
eleIsoDepositHcalFromTowersDepth2.ExtractorPSet.hcalDepth = cms.int32(2)


#this is a list of all the isolation deposits we need to convert to value maps
#so pat can use them
#step 1) get a list of isolations for which we need to do this
#step 2) convert the isoDeposits to value maps
#step 3) associate the value maps to the pat layer 0 output

patPatchElectronIsolationLabels = cms.VInputTag(
        cms.InputTag("eleIsoDepositEcalFromHits"),
        cms.InputTag("eleIsoDepositHcalFromTowersDepth1"),
        cms.InputTag("eleIsoDepositHcalFromTowersDepth2")
)

# read and convert to ValueMap<IsoDeposit> keyed to Candidate
patPatchElectronIsolationMaps = cms.EDFilter("MultipleIsoDepositsToValueMaps",
    collection   = cms.InputTag("pixelMatchGsfElectrons"),
    associations = patPatchElectronIsolationLabels,
)

#pat currently (219) needs the isolation deposits in a value map keyed to the layer0 output
#this does this
patLayer0PatchElectronIsolations = cms.EDFilter("CandManyValueMapsSkimmerIsoDeposits",
    collection   = cms.InputTag("allLayer0Electrons"),
    backrefs     = cms.InputTag("allLayer0Electrons"),
    commonLabel  = cms.InputTag("patPatchElectronIsolationMaps"),
    associations = patPatchElectronIsolationLabels,
)

#now we need a sequence to actually run
patchIsolSequence = cms.Sequence(eleIsoDepositEcalFromHits*eleIsoDepositHcalFromTowersDepth1*
                                         eleIsoDepositHcalFromTowersDepth2)
#patchIsolSequence = cms.Sequence(eleIsoDepositHcalFromTowersDepth1*
#                                         eleIsoDepositHcalFromTowersDepth2)

patLayer0PatchIsolSequence = cms.Sequence(patchIsolSequence*
                                         patPatchElectronIsolationMaps*
                                         patLayer0PatchElectronIsolations)

#now in theory if we run this sequence before the pat electron producer, we should be able
#to just load these isolations into  the pat electron





#now we need to run the patch isolation modules

#adding in patch isolations to the PAT electron
#this bit is the iso deposits
allLayer1Electrons.isoDeposits = cms.PSet(
        tracker = cms.InputTag("layer0ElectronIsolations","eleIsoDepositTk"),
        ecal    = cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositEcalFromHits"),
        hcal     = cms.InputTag("layer0ElectronIsolations","eleIsoDepositHcalFromTowers"),
        user    = cms.VInputTag( cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositHcalFromTowersDepth1"), 
                                 cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositHcalFromTowersDepth2"))
)

#now calculate the isolations and store them
allLayer1Electrons.isolation = cms.PSet(
    tracker = cms.PSet(
    # source IsoDeposit
    src = cms.InputTag("layer0ElectronIsolations","eleIsoDepositTk"),
    # parameters to compute isolation (Egamma POG defaults)
    deltaR = cms.double(0.3),
    vetos = cms.vstring('0.015', # inner radius veto cone
                        'Threshold(1.0)'),       # threshold on individual track pt
    skipDefaultVeto = cms.bool(True),
    ),
    ecal = cms.PSet(
    # source IsoDeposit
    src = cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositEcalFromHits"),
    # parameters to compute isolation (Egamma POG defaults)
    deltaR = cms.double(0.42),
    vetos = cms.vstring('EcalBarrel:0.040', 'EcalBarrel:RectangularEtaPhiVeto(-0.01,0.01,-0.5,0.5)',  # Barrel (|eta| < 1.479)
                        'EcalEndcaps:0.070','EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'),
    skipDefaultVeto = cms.bool(True),
    ),
    
    hcal = cms.PSet(
    src = cms.InputTag("layer0ElectronIsolations","eleIsoDepositHcalFromTowers"),
    # parameters to compute isolation (Egamma POG defaults)
    deltaR = cms.double(0.4),
    vetos = cms.vstring('0.1'),
    skipDefaultVeto = cms.bool(True),
    ),
    user = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositHcalFromTowersDepth1"),
    # parameters to compute isolation (Egamma POG defaults)
    deltaR = cms.double(0.4),
    vetos = cms.vstring('0.1'),
    skipDefaultVeto = cms.bool(True),
    ),

    cms.PSet(
    src = cms.InputTag("patLayer0PatchElectronIsolations","eleIsoDepositHcalFromTowersDepth2"),
    # parameters to compute isolation (Egamma POG defaults)
    deltaR = cms.double(0.4),
    vetos = cms.vstring(),
    skipDefaultVeto = cms.bool(True),
    )
    
    )#end user VPSet
)



patchPATSequence = cms.Sequence(patLayer0*
                               patLayer0PatchIsolSequence* #has to be before layer1
                               patLayer1)









