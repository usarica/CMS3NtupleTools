#Contains a reduced set of the core CMS2 makers for met template studies. Does not contain Gen or PAT makers
import FWCore.ParameterSet.Config as cms

#from CMS2.NtupleMaker.aSkimFilter_cfi import *
from CMS2.NtupleMaker.beamSpotMaker_cfi import *
from CMS2.NtupleMaker.beamHaloSequence_cff import *
#from CMS2.NtupleMaker.bTaggingSequence_cfi import *
#from CMS2.NtupleMaker.bTaggingTrkSequence_cfi import *
#from CMS2.NtupleMaker.bTagMaker_cfi import *
#from CMS2.NtupleMaker.bTagTrkMaker_cfi import *
from CMS2.NtupleMaker.electronMaker_cfi import *
from CMS2.NtupleMaker.electronSequence_cfi import *
from CMS2.NtupleMaker.eventMaker_cfi import *
from CMS2.NtupleMaker.eventSelectionMaker_cfi import *
from CMS2.NtupleMaker.hltMaker_cff import *
from CMS2.NtupleMaker.jetSequence_cff import *
from CMS2.NtupleMaker.jetMaker_cfi import *
from CMS2.NtupleMaker.jetToElAssMaker_cfi import *
from CMS2.NtupleMaker.jetToMuAssMaker_cfi import *
from CMS2.NtupleMaker.jptSequence_cff import *
from CMS2.NtupleMaker.jptMaker_cfi import *
from CMS2.NtupleMaker.jptToCaloJetAssMaker_cfi import *
from CMS2.NtupleMaker.l1Maker_cfi import *
from CMS2.NtupleMaker.metSequence_cff import *
from CMS2.NtupleMaker.metMaker_cfi import *
from CMS2.NtupleMaker.muonMaker_cfi import *
from CMS2.NtupleMaker.pfJetMaker_cfi import *
from CMS2.NtupleMaker.pfmetMaker_cfi import *
from CMS2.NtupleMaker.photonMaker_cfi import *
from CMS2.NtupleMaker.scMaker_cfi import *
from CMS2.NtupleMaker.tcmetMaker_cfi import *
from CMS2.NtupleMaker.trkJetMaker_cfi import *
from CMS2.NtupleMaker.trkJetSequence_cfi import *


CMS2Reco      = cms.Sequence(egammaElectronIDCMS2 * cms2CaloJetSequence * metCorSequence * cms2TrkJetSequence * cms2beamHaloSequence)

eventmakers   = cms.Sequence(beamSpotMaker * eventMaker * eventSelectionMaker * cms2beamHaloSequence)

trigmakers   = cms.Sequence(l1Maker * hltMakerSequence)

makers        = cms.Sequence(muonMaker * scMaker * electronMaker * photonMaker * jetMaker  * JPTCorrections * metMaker * tcmetMaker)

pflowmakers   = cms.Sequence(pfmetMaker * pfJetMaker )

cms2CoreSequence = cms.Sequence(CMS2Reco * eventmakers * trigmakers * makers * pflowmakers)
