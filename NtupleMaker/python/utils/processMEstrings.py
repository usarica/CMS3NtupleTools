import re

def processMEstrings(samplestr, MElist):
   newMElist=[]
   matchres = re.search("_M[0-9]*_",samplestr)
   massval=None
   if matchres is not None:
      matchstr = matchres.group()
      massval = matchstr
      massval = massval.replace("_M","")
      massval = massval.replace("_","")
      print "Mass value: {}".format(massval)
   for strme in MElist:
      newstrme = strme
      newstrme = newstrme.replace("<HMASS>",massval)
      print newstrme
      newMElist.append(newstrme)
   return newMElist

