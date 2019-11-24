import FWCore.ParameterSet.Config as cms

def replaceMVAValuesByRaw(obj, vpset):
   # Get values that are NOT supposed to be replaced
   noreplace=[]
   for pset in vpset:
      noreplacestr=pset.mvaLabel.value()
      if "RawValues" in noreplacestr:
         continue
      noreplacestrremove = noreplacestr.split("_")[-1]
      noreplacestr = noreplacestr.replace("_{}".format(noreplacestrremove),"")
      print "replaceMVAValuesByRaw: Will avoid replacing {}".format(noreplacestr)
      noreplace.append(noreplacestr)
   if not hasattr(obj,"modifierConfig"): return
   if not hasattr(obj.modifierConfig,"modifications"): return
   for par in obj.modifierConfig.modifications:
      cfgpars = []
      if hasattr(par,"electron_config"):
         cfgpars.append(par.electron_config)
      if hasattr(par,"photon_config"):
         cfgpars.append(par.photon_config)
      for cfgpar in cfgpars:
         for attribute, value in cfgpar.__dict__.items():
            if ("Values" in attribute) and not ("RawValues" in attribute):
               for attr,val in value.__dict__.items():
                  if ("productInstance" in attr) and not ("RawValues" in val):
                     newval=val
                     if not (val in noreplace):
                        newval = val.replace("Values", "RawValues")
                        print "Replacing {} with {}".format(val,newval)
                     else:
                        print "Keeping {} same as before".format(val)
                     setattr(value,attr,newval)
               setattr(cfgpar,attribute,value)

def getMVACutflowDictionary(obj):
   res=cms.VPSet()
   for physobj in obj.physicsObjectIDs:
      if hasattr(physobj,"idDefinition"):
         if hasattr(physobj.idDefinition,"cutFlow"):
            for par in physobj.idDefinition.cutFlow:
               if hasattr(par,"mvaValueMapName"):
                  mvaname=""
                  for attr,val in par.mvaValueMapName.__dict__.items():
                     if ("productInstance" in attr) and ("Values" in val):
                        mvaname=val
                        break
                  mvaname = cms.string("{}_{}".format(mvaname,physobj.idDefinition.idName.value().split('-')[-1]))
                  parres = cms.PSet(
                     mvaLabel = mvaname,
                     mvaCuts = par.mvaCuts
                     )
                  print "getMVACutflowDictionary: Adding label {}:".format(parres.mvaLabel.value())
                  print parres.mvaCuts
                  res.append(parres)
   return res
