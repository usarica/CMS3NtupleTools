import FWCore.ParameterSet.Config as cms

#from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *
from CMS2.NtupleMaker.electronDuplicateRemover_cfi import *
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

###
### Electron ID
###

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustLooseCMS2.src = "gsfElectrons"

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustTightCMS2.src = "gsfElectrons"
#eidRobustTightCMS2.robustEleIDCuts.barrel = [0.015, 0.0092, 0.020, 0.0025]
#eidRobustTightCMS2.robustEleIDCuts.endcap = [0.018, 0.025, 0.020, 0.0040]

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustHighEnergyCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustHighEnergyCMS2.src = "gsfElectrons"
#eidRobustHighEnergyCMS2.robustEleIDCuts.barrel = [0.050, 0.011, 0.090, 0.005]
#eidRobustHighEnergyCMS2.robustEleIDCuts.endcap = [0.100, 0.0275, 0.090, 0.007]

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidLooseCMS2.src = "gsfElectrons"
eidLooseCMS2.electronQuality = 'loose'

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidTightCMS2.src = "gsfElectrons"
eidTightCMS2.electronQuality = 'tight'

###
### Sequences
###

### Electron ID sequences
egammaElectronIDCMS2 = cms.Sequence(eidRobustLooseCMS2+eidRobustTightCMS2+eidRobustHighEnergyCMS2+eidLooseCMS2+eidTightCMS2)


