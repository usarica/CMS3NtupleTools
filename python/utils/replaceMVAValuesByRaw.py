def replaceMVAValuesByRaw(obj):
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
                     newval = val.replace("Values", "RawValues")
                     print "Replacing {} with {}".format(val,newval)
                     setattr(value,attr,newval)
               setattr(cfgpar,attribute,value)
