#!/usr/bin/env bash
#
# Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
# Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

#
# This script is used to run benchmarks for readjoiner and some related
# programs and save the session transcrips.
#

#
# Directories and filenames: (<TOP>= main working dir)
# [auto] means automatically created/mantained from this script
#
TOP="$HOME/local"
LOGSDIR="$HOME/logs"
READSETSDIR="$TOP/readsets"
REFSEQS="$TOP/refseqs"
REFRESULTS="$TOP/refresults"
#
# <TOP>/readsets ($READSETSDIR)     readsets, fasta format, single file
#                                   filename: <READSET>.reads.fas
#
# <TOP>/refseqs                     reference sequences for alignment test
#                                   filename: <READSET>.reference.fas
#
# <TOP>/refresults                  reference results files for diff test
#
# <TOP>/<READSET>/<PROGRAM>         [auto] results/transcripts
#
# <TOP>/run[1-9]                    benchmark scripts, should source this file
#
# <TOP>/run.include.local           sourced by this file
#
# ~/logs ($LOGSDIR)                 [auto] transcripts are collected here
# ~/logs/.NEXT                      [auto] next free log number
# ~/logs/RUNNING                    [auto] list of currently running benchmarks
# ~/logs/DONE                       [auto] list of completed benchmarks
#

#
# default parameters:
# (override these in run[1-9|include.local])
#

MINLEN=45                      # minimal SPM length

# binaries:
GTDIR="$TOP/gt64"              # genometools git repository
GT="$GTDIR/bin/gt"             # genometools binary
GTOPTS=""                      # this is use e.g. to pass the -debug opt
SGA="sga"                      # sga binary
EDENA="edena64"                # edena binary
LEAP="leap"                    # leap binary

# scripts:
SPACEPEAK="rdj-spacepeak.sh"   # spacepeak benchmark script
PARSETRANSCRIPT="rdj-parsetranscript.rb"        # transcript parser script
CLEANTRANSCRIPT="rdj-parsetranscript.rb -clean" # transcript cleaner script
CHECKCONTIGS="rdj-checkcontigs.rb"              # checkcontigs script
INITRAM=initRAM                # program to call to ensure a cold start

# scripts parameters:
RAM=0                          # argument of $INITRAM
SLEEP=0                        # amount of time to sleep after initRAM
MAXERR=2                       # maximal error accepted by checkcontigs script

# load local parameters
if [[ -e $TOP/run.include.local ]]; then source $TOP/run.include.local; fi

# -- private functions --

# $1: filename under $LOGSDIR (DONE|RUNNING)
# _LOGNUM/_PROG/_READSET must be set by the caller
function __add_to_list {
  DESC="$DESC "
  printf '%-6s%-13s%-15s%-50s%-10s%s\n' \
         "$_LOGNUM" "$_STEP" "$_READSET" "$DESC" `date +%y%m%d` `hostname` \
         >> $LOGSDIR/$1
}

# $1: filename under $LOGSDIR (DONE|RUNNING)
function __rm_from_list {
  ruby -e "\
  file = File.open(\"$LOGSDIR/$1\", \"r+\");\
  file.flock(File::LOCK_EX);\
  out = String.new;\
  file.each {|line| (out << line) if line !~ /$_LOGNUM/ };\
  file.pos = 0;\
  file.print out;\
  file.truncate(file.pos);\
  file.close"
}

# saves the next free lognumber in $_LOGNUM
# and updates $LOGSDIR/.NEXT
function __get_lognumber {
  _LOGNUM=`cat $LOGSDIR/.NEXT`
  __NEXTFREE=`echo $_LOGNUM + 1 | bc`
  echo $__NEXTFREE > $LOGSDIR/.NEXT
}

#
# dividing the run in __start, __run, __end allows one to
# inject additional code, such as remove files (run_leap)
# or output assembly statistics (_assembly functions)
#

# set $_PIPELINE, $_READSET, $_STEP, $_CMD, $_CMD_ARGS before calling this;
# after this you get a $_LOG filename where to append additional information
function __start {
  _PREVDIR=$PWD
  mkdir -p $TOP/$_READSET/$_PIPELINE
  _LOG=$TOP/$_READSET/$_PIPELINE/$_STEP.log
  cd $TOP/$_READSET/$_PIPELINE > /dev/null
  rm -f $_LOG
  echo "# benchmark transcript"    | tee -a $_LOG
  echo "# date: `date`"            | tee -a $_LOG
  echo "# hostname: `hostname`"    | tee -a $_LOG
  echo "# pipeline: $_PIPELINE"    | tee -a $_LOG
  echo "# step: $_STEP"            | tee -a $_LOG
  echo "# command: $_CMD"          | tee -a $_LOG
  echo "# arguments: $_CMD_ARGS"   | tee -a $_LOG
  echo -n "# working directory: "  | tee -a $_LOG
  echo "$PWD"                      | tee -a $_LOG
  __get_lognumber
  echo "Transcript number: $_LOGNUM"
  __add_to_list RUNNING
}

function __run {
  if [ $RAM   -gt 0 ]; then initRAM $RAM; fi
  if [ $SLEEP -gt 0 ]; then sleep $SLEEP; fi
  export GT_ENV_OPTIONS='-spacepeak -showtime'
  __RUN="$_CMD $_CMD_ARGS"
  echo Run $__RUN
  script -c "$SPACEPEAK $__RUN" -a $_LOG
  unset GT_ENV_OPTIONS
}

function __run_without_benchmark {
  __SPACEPEAK=$SPACEPEAK
  SPACEPEAK=""
  __run
  SPACEPEAK=$__SPACEPEAK
}

function __end {
  __TO=$LOGSDIR/$_LOGNUM.$_STEP.$_READSET
  cp $_LOG $__TO
  $CLEANTRANSCRIPT $__TO
  echo "Trascript file: $__TO"
  $PARSETRANSCRIPT $__TO
  __add_to_list  DONE
  __rm_from_list RUNNING
  cd $_PREVDIR
}

function __end_without_parsetranscript {
  __PARSETRANSCRIPT=$PARSETRANSCRIPT
  PARSETRANSCRIPT=true
  __end
  PARSETRANSCRIPT=$__PARSETRANSCRIPT
}

# $1: contigs filename
# requires $_LOG
function __append_seqstat_contigs {
  $GT seqstat -contigs $1 | tee -a $_LOG
}

# $1: git repository path
function __append_git_info {
  cd $1 > /dev/null
  echo "# git repository: $PWD"                 | tee -a $_LOG
  echo -n "# git HEAD is "                      | tee -a $_LOG
  git branch --contains HEAD -v | cut -c 3-     | tee -a $_LOG
  echo -n "# git master is "                    | tee -a $_LOG
  git log master --oneline -n 1                 | tee -a $_LOG
  cd - > /dev/null
}

# requires $GT, $GT, $GTDIR, $_LOG, $_CMD_ARGS
function __append_version_genometools {
  echo -n "# version: "                                | tee -a $_LOG
  $GT --version | head -n 1                            | tee -a $_LOG
  # compiler:
  echo -n "# "                                         | tee -a $_LOG
  $GT --version | tail -n 2 | head -n 1 | cut -c 6-    | tee -a $_LOG
  # compile flags:
  echo -n "# c"                                        | tee -a $_LOG
  $GT --version | tail -n 1 | cut -c 2-                | tee -a $_LOG
  __append_git_info $GTDIR
}

# requires $SGA, $_LOG, $_CMD_ARGS
function __append_version_sga {
  echo "# version: `$SGA --version | head -n 1`"       | tee -a $_LOG
}

# requires $EDENA, $_LOG, $_CMD_ARGS
function __append_version_edena {
  echo "# version: `$EDENA -v | head -n 1`"            | tee -a $_LOG
}

# requires $LEAP, $_LOG, $_CMD_ARGS
function __append_version_leap {
  echo "# version: n.a."                               | tee -a $_LOG
}

# requires _CONTIGS, _READSET and requirements of __start/__run/__end
function __checkcontigs {
  __REFERENCE=$REFSEQS/$_READSET.reference.fas
  _CMD=$CHECKCONTIGS
  _CMD_ARGS="$__REFERENCE $_CONTIGS $MAXERR"
  __start
  rm -rf $_CONTIGS.inexact $_CONTIGS.invalid
  __run_without_benchmark
  if [[ -e $_CONTIGS.invalid ]]; then
    $GT seqstat -distlen $_CONTIGS.invalid | tee -a $_LOG
  fi
  __end_without_parsetranscript
}

# -- public --

# run_* functions take one mandatory argument (readset name)
# and optionally a string to pass to the command;
# the minimal SPM length must be specified setting the $MINLEN variable

function run_suffixerator {
  _PIPELINE=sfx
  _STEP=sfx
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="suffixerator -mirrored -indexname $_READSET $2 -db $__READS "
  _CMD_ARGS+="-lcp -suf -ssp -v"
  __start
  __append_version_genometools
  __run
  __end
}

function run_readjoiner_correct {
  _PIPELINE=sfx
  _STEP=rdjC
  _READSET=$1
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="readjoiner correct $2 -ii $_READSET"
  __start
  __append_version_genometools
  __run
  __end
}

function run_readjoiner_prefilter {
  _PIPELINE=rdj
  _STEP=rdjP
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="readjoiner prefilter -v -readset $_READSET $2 -db $__READS"
  __start
  __append_version_genometools
  __run
  __end
}

function run_readjoiner_overlap {
  _PIPELINE=rdj
  _STEP=rdjO
  _READSET=$1
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="readjoiner overlap -l $MINLEN -v -readset $_READSET $2"
  __start
  __append_version_genometools
  __run
  __end
}

function run_readjoiner_assembly {
  _PIPELINE=rdj
  _STEP=rdjA
  _READSET=$1
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="readjoiner assembly -l $MINLEN -v -readset $_READSET $2"
  __start
  __append_version_genometools
  __run
  __end
}

function run_diff_encseq {
  _PIPELINE=rdj
  _STEP=diff_encseq
  _READSET=$1
  _REFRESULTSDIR=$HOME/data/refresults
  _CMD=diff
  _CMD_ARGS="-q -s $_READSET.esq "
  _CMD_ARGS+="$_REFRESULTSDIR/$_READSET/$_READSET.esq"
  __start
  __run_without_benchmark
  __end_without_parsetranscript
}

function run_diff_spm {
  _PIPELINE=rdj
  _STEP=diff_spm
  _READSET=$1
  _NOFSPMFILES=$2
  _REFRESULTSDIR=$HOME/data/refresults
  _CMD=diff
  _CMD_ARGS="-q -s"
  _CMD_ARGS+=" $_READSET.spm.txt.sorted"
  _CMD_ARGS+=" $_REFRESULTSDIR/$_READSET/$_READSET.spm.txt.sorted"
  __start
  rm -f $_READSET.spm.txt
  touch $_READSET.spm.txt
  for (( i=0; i<=$[$_NOFSPMFILES -1]; i++ ))
  do
    $GT readjoiner spmtest -test showlist \
      -readset $_READSET.$i >> $_READSET.spm.txt
  done
  sort $_READSET.spm.txt > $_READSET.spm.txt.sorted
  __run_without_benchmark
  __end_without_parsetranscript
}

function run_diff_contigs {
  _PIPELINE=rdj
  _STEP=diff_contigs
  _READSET=$1
  _REFRESULTSDIR=$HOME/data/refresults
  _CMD=diff
  _CMD_ARGS="-q -s $_READSET.contigs.raw "
  _CMD_ARGS+="$_REFRESULTSDIR/$_READSET/$_READSET.contigs.raw"
  __start
  $GT seq -recreate -showfasta -width 0 $_READSET.contigs.fas | \
    grep -v '^>' | sort > $_READSET.contigs.raw
  $GT clean
  __run_without_benchmark
  __end_without_parsetranscript
}

function run_encseq2spm {
  _PIPELINE=rdj
  _STEP=e2s
  _READSET=$1
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="encseq2spm -l $MINLEN -v -ii $_READSET $2"
  __start
  __append_version_genometools
  __run
  __end
}

function run_seqcorrect {
  _PIPELINE=rdj
  _STEP=seqc
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD="$GT $GTOPTS"
  _CMD_ARGS="dev seqcorrect -k $K -v -db $__READS -indexname $_READSET $2"
  __start
  __append_version_genometools
  __run
  __end
}

function run_edena_overlap {
  _PIPELINE=edena
  _STEP=ednO
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$EDENA
  _CMD_ARGS="-p $_READSET -r $__READS -M $MINLEN $2"
  __start
  __append_version_edena
  __run
  __end
}

function run_edena_assembly {
  _PIPELINE=edena
  _STEP=ednA
  _READSET=$1
  _CMD=$EDENA
  _CMD_ARGS="-p $_READSET -e $_READSET.ovl -m $MINLEN $2"
  __start
  __append_version_edena
  __run
  __append_seqstat_contigs ${_READSET}_contigs.fasta
  __end
}

function run_leap {
  _PIPELINE=leap
  _STEP=leap
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$LEAP
  _CMD_ARGS="$__READS $MINLEN output"
  __start
  __append_version_leap
  rm -rf output
  __run
  # first line of output of leap is blank, fix it:
  __firstline=`head -n 1 output/leap-all-contigs.fna`
  if [ "$__firstline" = "" ]; then
    mv output/leap-all-contigs.fna output/leap-all-contigs.fna.original
    tail +2 output/leap-all-contigs.fna.original > output/leap-all-contigs.fna
  fi
  __append_seqstat_contigs output/leap-all-contigs.fna
  __end
}

function run_sga_index {
  _PIPELINE=sga
  _STEP=sgaI
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$SGA
  _CMD_ARGS="index -v $2 $__READS"
  __start
  __append_version_sga
  __run
  __end
}

function run_sga_rmdup {
  _PIPELINE=sga
  _STEP=sgaR
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$SGA
  _CMD_ARGS="rmdup $2 -v $__READS"
  __start
  __append_version_sga
  __run
  __end
}

function run_sga_fm_merge {
  _PIPELINE=sga
  _STEP=sgaM
  _READSET=$1
  __READS=$_READSET.reads.rmdup.fa
  _CMD=$SGA
  _CMD_ARGS="fm-merge $2 -m $MINLEN -v $__READS"
  __start
  __append_version_sga
  __run
  __append_seqstat_contigs $_READSET.reads.rmdup.merged.fa
  __end
}

function run_sga_overlap {
  _PIPELINE=sga
  _STEP=sgaO
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$SGA
  _CMD_ARGS="overlap $2 -m $MINLEN -v $__READS"
  __start
  __append_version_sga
  __run
  __end
}

function run_sga_assembly {
  _PIPELINE=sga
  _STEP=sgaA
  _READSET=$1
  __READS=$READSETSDIR/$_READSET.reads.fas
  _CMD=$SGA
  _CMD_ARGS="assemble $2 -v $1.reads.asqg.gz"
  __start
  __append_version_sga
  __run
  __append_seqstat_contigs default-contigs.fa
  __end
}

function run_readjoiner_check {
  _PIPELINE=rdj
  _STEP=rdjC
  _READSET=$1
  _CONTIGS=$_READSET.contigs.fas
  __checkcontigs
}

function run_sga_check {
  _PIPELINE=sga
  _STEP=sgaC
  _READSET=$1
  _CONTIGS=default-contigs.fa
  __checkcontigs
}

function run_sga_fmmerge_check {
  _PIPELINE=sga
  _STEP=sgaCM
  _READSET=$1
  _CONTIGS=$_READSET.reads.rmdup.merged.fa
  __checkcontigs
}

function run_edena_check {
  _PIPELINE=edena
  _STEP=ednC
  _READSET=$1
  _CONTIGS=$_READSET
  _CONTIGS+=_contigs.fasta
  __checkcontigs
}

function run_leap_check {
  _PIPELINE=leap
  _STEP=leapC
  _READSET=$1
  _CONTIGS=output/leap-all-contigs.fna
  __checkcontigs
}

# run the next run[1-9] script (if exists)
function runnext {
  crr=`echo $0 | grep "run(\d)" -P -o | cut -c4`
  nxt=`echo $crr + 1 | bc`
  if [[ -e ./run$nxt ]]; then ./run$nxt; fi
}
