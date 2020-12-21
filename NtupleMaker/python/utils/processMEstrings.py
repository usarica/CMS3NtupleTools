import re

def getSampleMassValue(samplestr):
   matchres = re.search("_M[0-9]*_",samplestr)
   massval = None
   if matchres is not None:
      matchstr = matchres.group()
      massval = matchstr
      massval = massval.replace("_M","")
      massval = massval.replace("_","")
   if massval is not None:
      print "Mass value: {}".format(massval)
   else:
      raise RuntimeError("processMEstrings::getSampleMassValue: Cannot find sample mass value for {}".format(samplestr))
   return massval

#def getCMSBestMassReweightingLine(strme):
#   if "Name:SamplePropagator" in strme or "Name:SampleHypothesis" in strme or "Name:SampleDecayHypothesis" in strme or "Name:SampleProductionHypothesis" in strme or \
#      "Alias:SamplePropagator" in strme or "Alias:SampleHypothesis" in strme or "Alias:SampleDecayHypothesis" in strme or "Alias:SampleProductionHypothesis" in strme:
#      return strme
#   else:


def processMEstrings(samplestr, MElist):
   newMElist=[]
   matchres = re.search("_M[0-9]*_",samplestr)
   # Acquire mass value
   massval = getSampleMassValue(samplestr)
   for strme in MElist:
      newstrme = strme
      newstrme = newstrme.replace("<HMASS>",massval)

      print newstrme
      newMElist.append(newstrme)
   return newMElist

