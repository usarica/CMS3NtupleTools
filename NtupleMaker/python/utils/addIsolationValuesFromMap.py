import FWCore.ParameterSet.Config as cms

def addIsolationValuesFromMap(obj,isPhoton):
   if not hasattr(obj,"modifierConfig"): return
   if not hasattr(obj.modifierConfig,"modifications"): return
   for par in obj.modifierConfig.modifications:
      cfgpar = None

      if not hasattr(par, "modifierName"):
         continue
      if par.modifierName != "EGExtraInfoModifierFromFloatValueMaps":
         continue

      if not isPhoton and hasattr(par,"electron_config"):
         cfgpar = par.electron_config
      elif isPhoton and hasattr(par,"photon_config"):
         cfgpar = par.photon_config
      else:
         continue

      # Determine what to add first
      isWrongPset = False
      add_phoWorstChargedIsolation = True
      add_phoWorstChargedIsolationConeVeto = True
      add_phoWorstChargedIsolationConeVetoPVConstr = True
      for attribute, value in cfgpar.__dict__.items():
         if "ecalEnergy" in attribute:
            isWrongPset = True
            break
         if attribute == "phoWorstChargedIsolation":
            add_phoWorstChargedIsolation = False
         if attribute == "phoWorstChargedIsolationConeVeto":
            add_phoWorstChargedIsolationConeVeto = False
         if attribute == "phoWorstChargedIsolationConeVetoPVConstr":
            add_phoWorstChargedIsolationConeVetoPVConstr = False

      if isWrongPset:
         continue

      if add_phoWorstChargedIsolation:
         print "addIsolationValuesFromMap: Adding phoWorstChargedIsolation"
         setattr(cfgpar,"phoWorstChargedIsolation",cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"))
      if add_phoWorstChargedIsolationConeVeto:
         print "addIsolationValuesFromMap: Adding phoWorstChargedIsolationConeVeto"
         setattr(cfgpar,"phoWorstChargedIsolationConeVeto",cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationConeVeto"))
      if add_phoWorstChargedIsolationConeVetoPVConstr:
         print "addIsolationValuesFromMap: Adding phoWorstChargedIsolationConeVetoPVConstr"
         setattr(cfgpar,"phoWorstChargedIsolationConeVetoPVConstr",cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationConeVetoPVConstr"))
