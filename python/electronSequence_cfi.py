
import FWCore.ParameterSet.Config as cms

###
### Electron ID
###

from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

###
### Sequences
###

### Electron ID sequences
cms2EgammaElectronID = cms.Sequence(egammaIDLikelihood)

