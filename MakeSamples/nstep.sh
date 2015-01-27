#!/bin/bash

#echo TWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEPTWOSTEP
#echo twostep.sh called with args:
#echo $*
#echo NJob is $NJob
#echo MaxEvents is $MaxEvents
#echo RUNTIME_AREA is $RUNTIME_AREA
#echo
#echo PRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENV
#echo printenv
#printenv
#echo PRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENVPRINTENV
#echo
echo start step1 at `date`
echo
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml PSet.py $1
exit_code=$?
if [ $exit_code -ne 0 ]; then
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  echo @@@@ cmsRun exited step1 with error code $exit_code
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  exit $exit_code
fi
echo
echo done with step1 at `date`, starting step2
echo
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml step2_cfg.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  echo @@@@ cmsRun exited step2 with error code $exit_code
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  exit $exit_code
fi
echo
echo done with step2 at `date`, starting step3
echo
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml step3_cfg.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  echo @@@@ cmsRun exited step3 with error code $exit_code
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  exit $exit_code
fi
echo
echo done with step3 at `date`, starting step4
echo
cmsRun -j $RUNTIME_AREA/crab_fjr_$NJob.xml step4_cfg.py
exit_code=$?
if [ $exit_code -ne 0 ]; then
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  echo @@@@ cmsRun exited step4 with error code $exit_code
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  exit $exit_code
fi
echo
echo done with step4 at `date`
