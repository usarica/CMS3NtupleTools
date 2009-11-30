import FWCore.ParameterSet.Config as cms

beamFlaggerMaker = cms.EDProducer( "HLTLevel1GTSeed",

#
# see 
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HLTrigger/HLTfilters/python/hltLevel1GTSeed_cfi.py?revision=1.4&view=markup
#
#process.beamFlagger.L1TechTriggerSeeding = cms.bool(True)
#process.beamFlagger.L1SeedsLogicalExpression = cms.string('40 OR 41')


    # default: true
    #    seeding done via L1 trigger object maps, with objects that fired 
    #    only objects from the central BxInEvent (L1A) are used
    # if false:
    #    seeding is done ignoring if a L1 object fired or not, 
    #    adding all L1EXtra objects corresponding to the object types 
    #    used in all conditions from the algorithms in logical expression 
    #    for a given number of BxInEvent
    L1UseL1TriggerObjectMaps = cms.bool(True),
    #
    # option used forL1UseL1TriggerObjectMaps = False only
    # number of BxInEvent: 1: L1A=0; 3: -1, L1A=0, 1; 5: -2, -1, L1A=0, 1, 2
    L1NrBxInEvent = cms.int32(3),
    #
    # seeding done via technical trigger bits, if value is "true";
    # default: false (seeding via physics algorithms)
    #
    L1TechTriggerSeeding = cms.bool(True),
    #
    # seeding done with aliases for physics algorithms
    L1UseAliasesForSeeding = cms.bool(True),
    #
    # logical expression for the required L1 algorithms;
    # the algorithms are specified by name
    # allowed operators: "AND", "OR", "NOT", "(", ")"
    #
    # by convention, "L1GlobalDecision" logical expression means global decision
    # 
    L1SeedsLogicalExpression = cms.string('40 OR 41'),
    #
    # InputTag for the L1 Global Trigger DAQ readout record
    #   GT Emulator = gtDigis
    #   GT Unpacker = l1GtUnpack
    #
    #   cloned GT unpacker in HLT = gtDigis
    #
    L1GtReadoutRecordTag = cms.InputTag("gtDigis"),
    #
    # InputTag for L1 Global Trigger object maps
    #   only the emulator produces the object maps
    #   GT Emulator = gtDigis
    #
    #   cloned GT emulator in HLT = l1GtObjectMap
    #
    L1GtObjectMapTag = cms.InputTag("l1GtObjectMap"),
    #
    # InputTag for L1 particle collections (except muon)
    #   L1 Extra = l1extraParticles
    #
    L1CollectionsTag = cms.InputTag("l1extraParticles"),
    #
    # InputTag for L1 muon collection
    L1MuonCollectionTag = cms.InputTag("l1extraParticles"),
    #
    # saveTags for AOD book-keeping
    saveTags = cms.untracked.bool(True)


)


