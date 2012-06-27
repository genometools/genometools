#!/bin/bash
#
# Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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

# default values:
BWA=${BWA:=bwa}
GT=${GT:=gt}
SAMTOOLS=${SAMTOOLS:=samtools}
BWA_T=${BWA_T:=4}
PARAMS=${PARAMS:=""}

# command line parameters:
USG="<action> <genome.fas> <CDS.gff> <reads.fastq> [<mates.fastq>]"
ACTION=$1
GENOME=$2
ANNOTATION=$3
READS=$4
MATES=$5
MAP="${READS}.${GENOME}"
MAPM="${MATES}.${GENOME}"

# hardcoded values:
NOGFF='--'

function f_help {
  echo "Reference-based correction of homopolymer length in sequencing reads."
  echo
  echo "  Usage: $0 $USG"
  echo
  echo " - <action> is one of the following:"
  echo
  echo "     prepare       prepare the data necessary for the correction"
  echo "                   (index, mapping, sort gff and mapping results)"
  echo
  echo "     correct       run the homopolymer error correction tool;"
  echo "                   (for single-end reads, reads are output in a "
  echo "                   single file and are usually sorted differently "
  echo "                   than in the original file;"
  echo "                   if paired-end reads are used, the sorting "
  echo "                   order is conserved, at cost of a longer running "
  echo "                   time and an higher memory requirement)"
  echo
  echo "     clean         remove index, mapping results and intermediate files"
  echo "                   (BE CAREFUL, no warranty it works correctly)"
  echo
  echo "     run           prepare, correct and clean"
  echo
  echo "     stats         collect statistics for each correction position"
  echo "                   (useful for debugging and evaluation,"
  echo "                   it may be slow and require lot of memory)"
  echo
  echo " - <genome.fas> is a single sequence Fasta file which"
  echo "   contains the reference against which the correction is done"
  echo
  echo " - CDS.gff contains the annotation of genome.fas;"
  echo "   instead of a filename '$NOGFF' can be specified; in this case,"
  echo "   homopolymer are corrected in the whole sequence instead of only "
  echo "   in coding regions; this is however not reccomended"
  echo
  echo " - <reads.fastq> contains the reads to correct"
  echo
  echo " - <reads2.fastq> contains the mate pairs (if available)"
  echo
  echo "Further parameters can be optionally passed using env variables:"
  echo " - BWA          path to bwa binary             (default: $BWA)"
  echo " - SAMTOOLS     path to samtools binary        (default: $SAMTOOLS)"
  echo " - GT           path to GenomeTools binary     (default: $GT)"
  echo " - BWA_T        number of bwa threads          (default: $BWA_T)"
  echo " - PARAMS       further params for hopcorrect  (default: $PARAMS)"
}

function f_die {
  f_help
  exit $E_BADARGS
}

function f_index {
  echo "==== Index the genome..."
  echo
  $BWA index $GENOME
  echo
}

function f_map {
  echo "==== Map reads to the genome..."
  echo
  $BWA aln -t ${BWA_T} $GENOME $READS > $MAP.sai
  if [ "$MATES" != "" ]; then
    $BWA aln -t ${BWA_T} $GENOME $MATES > $MAPM.sai
    $BWA sampe $GENOME $MAP.sai $MAPM.sai $READS $MATES > $MAP.sam
  else
    $BWA samse $GENOME $MAP.sai $READS > $MAP.sam
  fi
  echo
}

function f_sortmap {
  echo "==== Sort the mapping results..."
  echo
  $SAMTOOLS view -Shu $MAP.sam > $MAP.bam
  $SAMTOOLS sort $MAP.bam sorted.$MAP
  echo
}

function f_sortgff {
  echo "==== Sort the annotation..."
  echo
  $GT gff3 -force -o sorted.$ANNOTATION -sort -retainids $ANNOTATION
  echo
}

function f_encode {
  echo
  echo "==== Encode genome in GtEncseq format..."
  echo
  $GT encseq encode -v $GENOME
  echo
}

function f_prepare {
  f_index
  f_map
  f_sortmap
  f_encode
  if [ "$ANNOTATION" != "$NOGFF" ]; then f_sortgff; fi
}

function f_correct {
  echo "==== Correct homopolymers..."
  echo
  if [ "$ANNOTATION" != "$NOGFF" ]; then
    HOP_PARAMS="-a sorted.$ANNOTATION $HOP_PARAMS"
  fi
  if [ "$MATES" != "" ]; then
    HOP_PARAMS="-outorder $READS $MATES $HOP_PARAMS"
  fi
  $GT dev hopcorrect -v \
      -r $GENOME \
      -m sorted.$MAP.bam \
      -o hop_$READS \
      $HOP_PARAMS
  echo
  echo "==== done "
  echo
  echo "corrected reads: hop_$READS"
  echo
}

function f_clean {
  rm -f $GENOME.amb $GENOME.ann $GENOME.bwt $GENOME.pac $GENOME.sa
  rm -f $MAP.sai $MAP.sam
  rm -f $MAP.bam sorted.$MAP.bam
  rm -f $GENOME.esq $GENOME.des $GENOME.md5 $GENOME.sds
  if [ "$ANNOTATION" != "$NOGFF" ]; then
    rm -f sorted.$ANNOTATION
  fi
  if [ "$MATES" != "" ]; then
    rm -f $MAPM.sai
  fi
}

function f_run {
  f_prepare
  f_correct
  f_clean
}

S_END='$5'
S_CHAR='$7'
count_all () {
  RETVAL=`cat $1 | grep '^'$2 | wc -l`
}
count_all_pos () {
  RETVAL=`cat $1 | awk '/^'$2'/ && ('$S_END' == "'$3'")' | wc -l`
}
count () {
  RETVAL=`cat $1 | awk '/^'$2'/ && ('$S_CHAR' == "'$3'")' | wc -l`
}
count_pos () {
  RETVAL=`cat $1 | \
  awk '/^'$2'/ && ('$S_END' == "'$3'") && ('$S_CHAR' == "'$4'")' | wc -l`
}
nof_seq_long () {
  RETVAL=`grep "^$2--$2" $1 | awk '{ print $2 }'`
  if [ "$RETVAL" == "" ]; then RETVAL=0; fi
}

function f_pos_stats {
  OUTFILE=$2
  MAX_READ_LEN=`grep 'maximum length' $3 | awk '{ print $4 }'`
  NOFSEQ=`grep 'sequences' $3 | awk '{ print $2 }'`
  nof_seq_long $3 1
  S=$RETVAL # number of sequence which are shorter than the current length
  P=`echo "$S / $NOFSEQ" | bc -l`
  rm -f $OUTFILE
  touch $OUTFILE
  echo \
    "# pos: last position in original read of corrected homopolymer (1-based)"\
      >> $OUTFILE
  echo \
    "# %reads: percent of reads which are at least as long as pos" >> $OUTFILE
  echo "# (ins|del): insertions, deletions in read" >> $OUTFILE
  echo "# (ins|del)[ACGT]: insertions, deletions in read of symbol" >> $OUTFILE
  H="# pos\t%reads\tedits\tins\tdel"
  H="$H\tinsA\tinsC\tinsG\tinsT\tdelA\tdelC\tdelG\tdelT"
  echo -e "$H" >> $OUTFILE
  count_all $1 I
  I=$RETVAL
  count_all $1 D
  D=$RETVAL
  A=$[$I+ $D]
  count $1 I a
  Ia=$RETVAL
  count $1 I c
  Ic=$RETVAL
  count $1 I g
  Ig=$RETVAL
  count $1 I t
  It=$RETVAL
  count $1 D a
  Da=$RETVAL
  count $1 D c
  Dc=$RETVAL
  count $1 D g
  Dg=$RETVAL
  count $1 D t
  Dt=$RETVAL
  echo -e "all\t100.00\t$A\t$I\t$D\t$Ia\t$Ic\t$Ig\t$It\t$Da\t$Dc\t$Dg\t$Dt" \
    >> $OUTFILE
  for ((i = 2; i < $MAX_READ_LEN; i++)); do
    nof_seq_long $3 $i
    S=$[$RETVAL+$S]
    P=`echo "(1.0 - ($S / $NOFSEQ)) * 100" | bc -l`
    count_all_pos $1 I $i
    I=$RETVAL
    count_all_pos $1 D $i
    D=$RETVAL
    A=$[$I+ $D]
    count_pos $1 I $i a
    Ia=$RETVAL
    count_pos $1 D $i a
    Da=$RETVAL
    count_pos $1 I $i c
    Ic=$RETVAL
    count_pos $1 D $i c
    Dc=$RETVAL
    count_pos $1 I $i g
    Ig=$RETVAL
    count_pos $1 D $i g
    Dg=$RETVAL
    count_pos $1 I $i t
    It=$RETVAL
    count_pos $1 D $i t
    Dt=$RETVAL
    printf "$i\t%.2f\t$A\t$I\t$D\t$Ia\t$Ic\t$Ig\t$It\t$Da\t$Dc\t$Dg\t$Dt\n" $P \
      >> $OUTFILE
  done
}

function f_corrs_per_read {
  ID='$11'
  grep -v '^#' $1 | awk '{ print '$ID' }' | uniq -c > $MAP~tmpfile
  echo "# corrected reads: `cat $MAP~tmpfile | wc -l`" >> $1
  t=0
  for ((i=1;i<10;i++)); do
    c=`awk ' $1 == '$i $MAP~tmpfile | wc -l`
    if [ $c == 0  ]; then
      break
    else
      t=$[$t+($c*$i)]
      echo -e "# reads corrected $i time(s): $c" >> $1
    fi
  done
  echo -e "# total corrections: $t" >> $1
  rm $MAP~tmpfile
}

function f_stats {
  echo "==== Correcting and collecting statistics..."
  echo
  if [ "$ANNOTATION" != "$NOGFF" ]; then
    AOPT="-a sorted.$ANNOTATION"
  else
    AOPT=""
  fi
  $GT dev hopcorrect -v -r $GENOME -m sorted.$MAP.bam $AOPT \
      -o hop_$READS -stats $HOP_PARAMS >| $MAP.hop_stats
  echo
  f_corrs_per_read $MAP.hop_stats
  echo "==== done "
  echo
  echo "correction statistics:                $MAP.hop_stats"
  echo
}

function f_dists  {
  echo "==== Computing distributions..."
  echo
  $GT seqstat -distlen -b 1 $READS $MATES >| $READS.lendist
  $GT seqstat -distlen -b 1 hop_$READS >| hop_$READS.lendist
  f_pos_stats $MAP.hop_stats $MAP.hop_pos_stats $READS.lendist
  echo
  echo "==== done "
  echo
  echo "original reads length distribution:   $READS.lendist"
  echo "corrected reads length distribution:  hop_$READS.lendist"
  echo "corrections per read position:        $MAP.hop_pos_stats"
  echo
}

if [ $# -lt 4 -o $# -gt 5 ]; then f_die; fi
case "$ACTION" in
  'prepare')   f_prepare   ;;
  'correct')   f_correct   ;;
  'clean')     f_clean     ;;
  'run')       f_run       ;;
  'stats')     f_stats     ;;
  'dists')     f_dists     ;;
  *)           f_die       ;;
esac
