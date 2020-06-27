#!/bin/bash

date=$1
period=$2
prodVersion=$3
script=produceLeptonEfficiencies_RooFit.cc
function=getEfficiencies
jobdate="${date}_LeptonTnP_RooFit"
arguments='"<period>","<prodVersion>","<strdate>",<is_ee>,<eeGapCode>,<resPdgId>,"<systOptions>",<minPt_tag>,<fit_low>,<fit_high>'
arguments="${arguments/<period>/$period}"
arguments="${arguments/<strdate>/$date}"
arguments="${arguments/<prodVersion>/$prodVersion}"

declare -a csvList

declare -a minPtTags_ee
declare -a minPtTags_mumu
declare -a systOptions

systOptions=( "" "ALTBkg" "ALTBkg2" "TightTag" "TightTag.ALTBkg" "TightTag.ALTBkg2" "MC_2l2nu" "MC_4l" "PUDn" "PUUp" )

if [[ "$period" == "2016" ]]; then
  minPtTags_ee=( 28 30 )
  minPtTags_mumu=( 27 29 )
elif [[ "$period" == "2017" ]]; then
  minPtTags_ee=( 38 40 )
  minPtTags_mumu=( 30 32 )
elif [[ "$period" == "2018" ]]; then
  minPtTags_ee=( 35 37 )
  minPtTags_mumu=( 27 29 )
fi


# Do ee events first
for eeGapCode in -1 0 1; do
for minPtTag in "${minPtTags_ee[@]}"; do
for syst in "${systOptions[@]}"; do

  strargscore="${arguments}"
  strargscore="${strargscore/<is_ee>/true}"
  strargscore="${strargscore/<eeGapCode>/${eeGapCode}}"
  strargscore="${strargscore/<minPt_tag>/${minPtTag}}"
  strargscore="${strargscore/<systOptions>/${syst}}"

  # Do Z mass windows
  for iw in {1..3}; do
    if [[ "$syst" == "MC"* ]] && [[ $iw -gt 1 ]]; then
      continue
    fi
    if [[ "$syst" == "PU"* ]] && [[ $iw -gt 1 ]]; then
      continue
    fi
    if [[ "$syst" == *"ALTBkg"* ]] && [[ $iw -ne 3 ]]; then
      continue
    fi

    let fit_low=60
    let fit_high=120
    if [[ $iw -eq 2 ]]; then
      let fit_low=65
      let fit_high=115
    elif [[ $iw -eq 3 ]]; then
      let fit_low=70
      let fit_high=110
    fi

    strargs="${strargscore}"
    strargs="${strargs/<resPdgId>/23}"
    strargs="${strargs/<fit_low>/${fit_low}}"
    strargs="${strargs/<fit_high>/${fit_high}}"

    echo submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done

done
done
done

# Do mumu events first
for eeGapCode in -1; do
for minPtTag in "${minPtTags_mumu[@]}"; do
for syst in "${systOptions[@]}"; do

  strargscore="${arguments}"
  strargscore="${strargscore/<is_ee>/false}"
  strargscore="${strargscore/<eeGapCode>/${eeGapCode}}"
  strargscore="${strargscore/<minPt_tag>/${minPtTag}}"
  strargscore="${strargscore/<systOptions>/${syst}}"

  # Do Z mass windows
  for iw in {1..3}; do
    if [[ "$syst" == "MC"* ]] && [[ $iw -gt 1 ]]; then
      continue
    fi
    if [[ "$syst" == "PU"* ]] && [[ $iw -gt 1 ]]; then
      continue
    fi
    if [[ "$syst" == *"ALTBkg"* ]] && [[ $iw -ne 3 ]]; then
      continue
    fi

    let fit_low=60
    let fit_high=120
    if [[ $iw -eq 2 ]]; then
      let fit_low=65
      let fit_high=115
    elif [[ $iw -eq 3 ]]; then
      let fit_low=70
      let fit_high=110
    fi

    strargs="${strargscore}"
    strargs="${strargs/<resPdgId>/23}"
    strargs="${strargs/<fit_low>/${fit_low}}"
    strargs="${strargs/<fit_high>/${fit_high}}"

    echo submitCMS3AnalysisProduction.sh script="${script}" function="${function}" arguments="${strargs}" date="${jobdate}"
  done

done
done
done
