#!/bin/bash

sleeptime=0.1

if [ $# -eq 0 ]; then
  echo "Usage: $0 <command> [args]"
  echo
  echo "The following information is polled each $sleeptime seconds"
  echo "from /proc/[pid]/status:"
  echo
  echo "  VmPeak: Peak virtual memory size."
  echo "  VmSize: Virtual memory size."
  echo "  VmLck: Locked memory size."
  echo "  VmHWM: Peak resident set size (\"high water mark\")."
  echo "  VmRSS: Resident set size."
  echo "  VmData, VmStk, VmExe: Size of data, stack, and text segments."
  echo "  VmLib: Shared library code size."
  echo "  VmPTE: Page table entries size (since Linux 2.6.10)."
  echo
  echo "The command is run under /usr/bin/time."
  exit
fi

# code inspired by:
#   http://stackoverflow.com/questions/1080461/
#         /peak-memory-measurement-of-long-running-process-in-linux
function __measure_space_peak {
  types="Peak Size Lck HWM RSS Data Stk Exe Lib PTE"
  declare -A maxVm
  for vm in $types; do maxVm[$vm]=0; done
  ppid=$$
  /usr/bin/time $@ &
  tpid=`pgrep -P ${ppid} -n -f time`
  if [[ ${tpid} -ne "" ]]; then
    pid=`pgrep -P ${tpid} -n -f $1` # $! may work here but not later
  fi
  declare -A Vm
  while [[ ${tpid} -ne "" ]]; do
    for vm in $types; do
      if [[ ${pid} -ne "" ]]; then
        Vm[$vm]=`cat /proc/${pid}/status 2> /dev/null \
        | grep Vm${vm} | awk '{print $2}'`
        if [[ ${Vm[$vm]} -gt ${maxVm[$vm]} ]]; then
          maxVm[$vm]=${Vm[$vm]}
        fi
      fi
    done
    sleep $sleeptime
    savedtpid=${tpid}
    tpid=`pgrep -P ${ppid} -n -f time`
  done
  wait ${savedtpid} # don't wait, job is finished
  exitstatus=$?   # catch the exit status of wait, the same of $@
  echo "Memory usage for $@:"
  for vm in $types; do echo "  Vm$vm: ${maxVm[$vm]} kB"; done
  echo "Exit status: ${exitstatus}"
}
__measure_space_peak $*
