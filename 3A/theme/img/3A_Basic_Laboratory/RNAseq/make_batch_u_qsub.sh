#!/bin/bash

Ssymbol="[$"
Esymbol="]"
function AssignParams () {
  : '
  @params params.json     {"qsub_q": "u-debug"}
  @params templates.sh    PBS -q [${n_node}]
  @return batch_u_qsub.sh PBS -q u-debug
  '

  line="${1}"
  while true
  do
    Spos=`echo "${line}" | awk -v pattern=${Ssymbol} '{print index($0, pattern)}'`
    Epos=`echo "${line}" | awk -v pattern=${Esymbol} '{print index($0, pattern)}'`
    if [ ${Spos} -ne 0 -a ${Epos} -ne 0 -a ${Spos} -lt ${Epos} ]; then
      PREFIX=`echo "${line:0:$(($Spos-1))}"`
      VARIABLE=`echo "${line}" | cut -c "$(($Spos+1))-$(($Epos-1))"`
      CONTENT=`eval echo ${VARIABLE}`
      SUFFIX=`echo "${line:$Epos}"`
      line="${PREFIX}${CONTENT}${SUFFIX}"
    else
      break
    fi
  done
  echo "${line}" >> $OUTPUT_FILE
}

#=== START ===
echo -n "Path to parameter file (.json) : "
read PARAMS_FILE
# Convert each element of json to a variable.
KEYS=`jq -r 'keys[]' $PARAMS_FILE`
for key in $KEYS; do
  eval $key=`jq -r .${key} $PARAMS_FILE`
done

OUTPUT_FILE="batch_u_qsub.sh"
if [ -e ${OUTPUT_FILE} ]; then
  rm ${OUTPUT_FILE}
fi
if [ ! -d ${ID} ]; then
  mkdir ${ID}
fi
echo -n "Path to templates file (.sh)   : "
read TEMPLATES_FILE

PRE_IFS=$IFS
IFS=$'\n'
for line in `cat ${TEMPLATES_FILE}`
do
  AssignParams ${line}
done
IFS=$PRE_IFS
