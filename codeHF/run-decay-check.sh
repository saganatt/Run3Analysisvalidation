#!/bin/bash

source "/home/mkabus/Run3Analysisvalidation/exec/utilities.sh"

DIR="lc-dl-checks"
CONFIG_INPUT="config_input.sh"

for dl in `seq 0. 0.002 0.002` ; do
  echo "dl ${dl}"
  CUR_JSON="dpl-config_run3_edit_dl${dl}.json"
  CUR_DIR="${DIR}/dl${dl}"
  rm -rf "${CUR_DIR}"
  sed -i "s/JSON=\".*\"/JSON=\"${CUR_JSON}\"/g" "${CONFIG_INPUT}" || ErrExit "Could not edit JSON"
  bash runtest.sh > debug_dl.txt 2>&1 || ErrExit "Runtest failed"
  ./move-o2-results.sh "${CUR_DIR}" || ErrExit "Move results failed"
  cp "${CUR_JSON}" "${CUR_DIR}" || ErrExit "Move JSON failed"
done
