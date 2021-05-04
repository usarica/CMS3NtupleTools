#!/bin/bash

date=$1
period=$2
prodVersion=$3
fitVersion=$4
strFilter=""
if [[ "$5" != "" ]]; then
  strFilter="$5"
fi
script=produceLeptonEfficiencies_RooFit.cc
functions=( summarizeFits plotAllFits )
jobdate="${date}_LeptonTnP_FitSummary"

for function in "${functions[@]}"; do
  arguments=''
  if [[ "$function" == "summarizeFits" ]]; then
    arguments='"<period>","<strdate>","<fitVersion>",<is_ee>,<eeGapCode>,<resPdgId>'
  else
    arguments='"<period>","<prodVersion>","<fitVersion>",<is_ee>,<eeGapCode>,<resPdgId>,"<strFilter>"'
  fi
  arguments="${arguments/<period>/$period}"
  arguments="${arguments/<strdate>/$date}"
  arguments="${arguments/<fitVersion>/$fitVersion}"
  arguments="${arguments/<prodVersion>/$prodVersion}"
  arguments="${arguments/<resPdgId>/23}"
  arguments="${arguments/<strFilter>/${strFilter}}"

  # Do ee events first
  for eeGapCode in {-1..1}; do
    strargs="${arguments}"
    strargs="${strargs/<is_ee>/true}"
    strargs="${strargs/<eeGapCode>/${eeGapCode}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"

  done

  # Do mumu events next
  for eeGapCode in -1; do
    strargs="${arguments}"
    strargs="${strargs/<is_ee>/false}"
    strargs="${strargs/<eeGapCode>/${eeGapCode}}"

    submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"

  done

done
