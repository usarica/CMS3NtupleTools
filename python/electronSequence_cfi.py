
import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

###
### Electron ID
###

###
### Here version V01 is used for all IDs.
### This version was tuned in 2_2_X releases and represents 
### the latest official Egamma electron ID.
###

### Robust loose
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustLooseCMS2.electronIDType = 'robust'
eidRobustLooseCMS2.electronQuality = 'loose'
eidRobustLooseCMS2.electronVersion = 'V01'

### Robust tight
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustTightCMS2.electronIDType = 'robust'
eidRobustTightCMS2.electronQuality = 'tight'
eidRobustTightCMS2.electronVersion = 'V01'

### Robust high energy (HEEP)
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustHighEnergyCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustHighEnergyCMS2.electronIDType = 'robust'
eidRobustHighEnergyCMS2.electronQuality = 'highenergy'
eidRobustHighEnergyCMS2.electronVersion = 'V01'

### Loose catagorised (Sani)
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidLooseCMS2.electronIDType 	= 'classbased'
eidLooseCMS2.electronQuality 	= 'loose'
eidLooseCMS2.electronVersion 	= 'V01'

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidTightCMS2.electronIDType     = 'classbased'
eidTightCMS2.electronQuality    = 'tight'
eidTightCMS2.electronVersion   	= 'V01'

###
### Sequences
###

### Electron ID sequences
egammaElectronIDCMS2 = cms.Sequence(eidRobustLooseCMS2+eidRobustTightCMS2+eidRobustHighEnergyCMS2+eidLooseCMS2+eidTightCMS2)

