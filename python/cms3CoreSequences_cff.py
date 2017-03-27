#Contains the core CMS2 makers. Does not contain Gen or PAT makers
import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

from CMS3.NtupleMaker.pftauMaker_cfi               import *
from CMS3.NtupleMaker.electronMaker_cfi            import *
from CMS3.NtupleMaker.eventMaker_cfi               import *
from CMS3.NtupleMaker.hltMaker_cff                 import *
from CMS3.NtupleMaker.hypDilepMaker_cfi            import *
from CMS3.NtupleMaker.muonMaker_cfi                import *
from CMS3.NtupleMaker.photonMaker_cfi              import *
from CMS3.NtupleMaker.secVertexMaker_cfi           import *
from CMS3.NtupleMaker.vertexMaker_cfi              import *
from CMS3.NtupleMaker.pfJetMaker_cfi               import *
from CMS3.NtupleMaker.subJetMaker_cfi              import *
from CMS3.NtupleMaker.fastJetSequence_cff          import *
if configProcessName.isFastSim == False:
    from CMS3.NtupleMaker.metFilterMaker_cfi           import *
else:
    from CMS3.NtupleMaker.sParmMaker_cff               import * # doesn't always get loaded
