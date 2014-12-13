#Contains the core CMS2 makers. Does not contain Gen or PAT makers
import FWCore.ParameterSet.Config as cms

#from CMS3.NtupleMaker.aSkimFilter_cfi              import *
from CMS3.NtupleMaker.beamSpotMaker_cfi            import *
from CMS3.NtupleMaker.bTaggingSequence_cfi         import *
from CMS3.NtupleMaker.bTagJPTSequence_cfi          import *
from CMS3.NtupleMaker.bTaggingTrkSequence_cfi      import *
from CMS3.NtupleMaker.bTagMaker_cfi                import *
from CMS3.NtupleMaker.bTagJPTJetMaker_cfi          import *
from CMS3.NtupleMaker.bTagTrkMaker_cfi             import *
from CMS3.NtupleMaker.tauTaggingSequence_cfi       import *
from CMS3.NtupleMaker.pftauMaker_cfi               import *
from CMS3.NtupleMaker.elCaloIsoSequence_cff        import *
from CMS3.NtupleMaker.elTkJuraIsoMaker_cfi         import *
from CMS3.NtupleMaker.electronMaker_cfi            import *
from CMS3.NtupleMaker.electronSequence_cfi         import *
from CMS3.NtupleMaker.elToJetAssMaker_cfi          import *
from CMS3.NtupleMaker.elToMuAssMaker_cfi           import *
from CMS3.NtupleMaker.eventMaker_cfi               import *
from CMS3.NtupleMaker.gsfTrackMaker_cfi            import *
from CMS3.NtupleMaker.hcalNoiseSummaryMaker_cfi    import *
from CMS3.NtupleMaker.hltMaker_cff                 import *
from CMS3.NtupleMaker.hypDilepMaker_cfi            import *
from CMS3.NtupleMaker.hypDilepVertexMaker_cfi      import *
from CMS3.NtupleMaker.hypTrilepMaker_cfi           import *
from CMS3.NtupleMaker.hypQuadlepMaker_cfi          import *
from CMS3.NtupleMaker.hypIsoMaker_cfi              import *
from CMS3.NtupleMaker.jetSequence_cff              import *
from CMS3.NtupleMaker.jetMaker_cfi                 import *
from CMS3.NtupleMaker.jetToElAssMaker_cfi          import *
from CMS3.NtupleMaker.jetToMuAssMaker_cfi          import *
from CMS3.NtupleMaker.jptMaker_cfi                 import *
from CMS3.NtupleMaker.l1Maker_cfi                  import *
from CMS3.NtupleMaker.luminosityMaker_cfi          import *
#from CMS3.NtupleMaker.metSequence_cff              import *
from CMS3.NtupleMaker.metMaker_cfi                 import *
from CMS3.NtupleMaker.muonMaker_cfi                import *
from CMS3.NtupleMaker.muToElsAssMaker_cfi          import *
from CMS3.NtupleMaker.muToJetAssMaker_cfi          import *
from CMS3.NtupleMaker.photonMaker_cfi              import *
from CMS3.NtupleMaker.recoErrorLogMaker_cfi        import *
from CMS3.NtupleMaker.recoConversionMaker_cfi      import *
from CMS3.NtupleMaker.scMaker_cfi                  import *
from CMS3.NtupleMaker.secVertexMaker_cfi           import *
from CMS3.NtupleMaker.tcmetSequence_cff            import *
from CMS3.NtupleMaker.trackMaker_cfi               import *
from CMS3.NtupleMaker.trackToElsAssMaker_cfi       import *
from CMS3.NtupleMaker.trackToMuonAssMaker_cfi      import *
from CMS3.NtupleMaker.trkJetMaker_cfi              import *
from CMS3.NtupleMaker.trkToVtxAssMaker_cfi         import *
from CMS3.NtupleMaker.vertexMaker_cfi              import *
from CMS3.NtupleMaker.beamHaloMaker_cfi            import *
from CMS3.NtupleMaker.fastJetSequence_cff          import *
from CMS3.NtupleMaker.pfJetMaker_cfi               import *
from CMS3.NtupleMaker.subJetMaker_cfi               import *
from CMS3.NtupleMaker.muToTrigAssMaker_cfi         import *
from CMS3.NtupleMaker.elToTrigAssMaker_cfi         import *
#from CMS3.NtupleMaker.ecalDRFilterMaker_cff        import *
#from CMS3.NtupleMaker.ecalTPFilterMaker_cff        import *
#from CMS3.NtupleMaker.eeBadRecovMaker_cff          import *
from CMS3.NtupleMaker.metFilterSequence_cff        import *
from CMS3.NtupleMaker.metFilterMaker_cfi           import *

from CMS3.NtupleMaker.cms2PFSequence_cff           import *

from CMS3.NtupleMaker.sParmMaker_cff               import * # doesn't always get loaded

#CMS2Reco         = cms.Sequence( cms2JetSequence * metCorSequence * CMS2Btagging * CMS2TrkBtagging * CMS2JPTBtagging )
CMS2Reco         = cms.Sequence( cms2JetSequence * CMS2Btagging * CMS2TrkBtagging * CMS2JPTBtagging )
eventmakers      = cms.Sequence( beamSpotMaker * vertexMaker * vertexMakerWithBS * eventMaker * hcalNoiseSummaryMaker * cms2EgammaElectronID )
eventmakerswsparm= cms.Sequence( eventmakers * sParmMaker ) # build up alternate sequence to be swapped in main config file

trigmakers       = cms.Sequence( l1Maker * hltMakerSequence )

makers           = cms.Sequence( trackMaker * gsfTrackMaker * muonMaker * scMaker * fastJetSequence * electronMaker * photonMaker * jetMaker * jptMaker * trkJetMaker * pfJetMaker * subJetMaker * metMaker * 
                                 tcmetMaker * luminosityMaker * recoErrorLogMaker * beamHaloMaker * recoConversionMaker * cms2MetFilterSequence * metFilterMaker )

assmakers        = cms.Sequence( jetToMuAssMaker * jetToElAssMaker * muToElsAssMaker * muToJetAssMaker * elToMuAssMaker * elToJetAssMaker * trackToMuonAssMaker * trackToElsAssMaker * trkToVtxAssMaker * muToTrigAssMaker * elToTrigAssMaker)
hypmakers        = cms.Sequence( hypDilepMaker * hypDilepVertexMaker * hypTrilepMaker * hypQuadlepMaker )
othermakers      = cms.Sequence( elCaloIsoSequence * elTkJuraIsoMaker * bTagMaker *  bTagTrkMaker * bTagJPTJetMaker )#* pftauMaker )
#cms2CoreSequence = cms.Sequence( CMS2Reco * eventmakers * trigmakers * makers * assmakers * hypmakers * CMS2PFtautagging * othermakers )
cms2CoreSequence = cms.Sequence( CMS2Reco *  eventmakers * trigmakers * makers * assmakers * hypmakers * othermakers )

## the CMS2PFtautagging need the fastJetSequence before so the order is makers / CMS2PFtautagging / othermakers 
