#!/bin/bash

chkdir=$1

for f in $(find $chkdir -name condor.sub); do
  d=${f//\/condor.sub}
  cd $d
  
  arguments=""
  initialdir=""
  while IFS='' read -r line || [[ -n "$line" ]]; do
    llow="$(echo $line | awk '{print tolower($0)}')"
    if [[ "$llow" == "initialdir"*"="* ]];then
      fcnargname="$line"
      fcnargname="${fcnargname#*=}"
      fcnargname="${fcnargname#* }"
      initialdir="$fcnargname"
    elif [[ "$llow" == "arguments"*"="* ]];then
      fcnargname="$line"
      fcnargname="${fcnargname#*=}"
      fcnargname="${fcnargname#* }"
      arguments="$fcnargname"
    fi
  done < "condor.sub"

  RUNFILE=main_pset.py
  RUN_CMD=$(runGenericExecutable.py --executable="$RUNFILE" --command="$arguments" --dry)
  if [[ "$RUN_CMD" == "Running "* ]];then
    RUN_CMD=${RUN_CMD//"Running "}
  else
    echo "Run command ${RUN_CMD} is invalid."
    exit 1
  fi

  INFILENAMES=""
  OUTFILENAME=""
  fcnarglist=($(echo $RUN_CMD))
  for fargo in "${fcnarglist[@]}";do
    farg="${fargo//\"}"
    fargl="$(echo $farg | awk '{print tolower($0)}')"
    if [[ "$fargl" == "output="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      OUTFILENAME="${fcnargname}"
    elif [[ "$fargl" == "inputs="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      INFILENAMES="${fcnargname}"
    fi
  done

  cd - &> /dev/null

  infilelist=($(echo ${INFILENAMES//,/ }))
  if [[ ${#infilelist[@]} -eq 1 ]];then
    continue
  fi

  let ctr=0
  for infile in "${infilelist[@]}";do
    newd="${d}_recovery${ctr}"
    newinitialdir="${initialdir}_recovery${ctr}"
    newoutfile="${OUTFILENAME/.root/_recovery${ctr}.root}"
    mkdir -p $newd/Logs
    cp ${d}/condor.sub ${newd}/condor.sub
    ln -sf $(readlink -f ${d}/cms3ntuplemaker.tar) ${newd}/cms3ntuplemaker.tar
    sed -i "s|${OUTFILENAME}|${newoutfile}|g" ${newd}/condor.sub
    sed -i "s|${initialdir}|${newinitialdir}|g" ${newd}/condor.sub
    sed -i "s|${INFILENAMES}|${infile}|g" ${newd}/condor.sub
    echo "Created new submission directory $newd with the following replacements:"
    echo "- initialdir = $newinitialdir"
    echo "- inputs = $infile"
    echo "- output = $newoutfile"
    let ctr=${ctr}+1
  done

  rm -rf $d
  
done
