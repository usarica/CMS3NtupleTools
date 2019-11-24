#Contains the core CMS3 makers. Does not contain Gen or PAT makers
import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

from CMS3.NtupleMaker.electronMaker_cfi import *
from CMS3.NtupleMaker.lumiFilter_cfi import *
from CMS3.NtupleMaker.hltMaker_cfi import *
from CMS3.NtupleMaker.muonMomentumCorrector_cfi import *
from CMS3.NtupleMaker.muonMaker_cfi  import *
from CMS3.NtupleMaker.photonMaker_cfi import *
from CMS3.NtupleMaker.secVertexMaker_cfi import *
from CMS3.NtupleMaker.vertexMaker_cfi import *
from CMS3.NtupleMaker.pfJetMaker_cfi import *
from CMS3.NtupleMaker.metFilterMaker_cfi import *
from CMS3.NtupleMaker.pfmetMaker_cfi import *
from CMS3.NtupleMaker.pftauMaker_cfi import *
from CMS3.NtupleMaker.pfCandidateMaker_cfi import *
from CMS3.NtupleMaker.isoTrackMaker_cfi import *
