import FWCore.ParameterSet.Config as cms


from MitEdm.Producers.conversionRejection_cff import *
from CMS2.NtupleMaker.mitConversionMaker_cfi import *
from CMS2.NtupleMaker.elToMITConvAssMaker_cfi import * 

mvfConversionRemoval.useRhoMin = False

mitConversionSequence = cms.Sequence(conversionRejection * mitConversionMaker * elToMITConvAssMaker)
