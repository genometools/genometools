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

# hardcoded values:
NOGFF='--'

# default values:
BWA=${BWA:=bwa}
GT=${GT:=gt}
SAMTOOLS=${SAMTOOLS:=samtools}
HOP_PARAMS=${HOP_PARAMS:=""}
BWA_PARAMS=${BWA_PARAMS:="-t 4"}
USESW=${USESW:=""}

# command line parameters:
USG="<action> <genome.fas> (<CDS.gff>|$NOGFF) <reads.fastq> [<mates.fastq>]"
ACTION=$1
if [[ "$ACTION" =~ "eval" ]]; then
  USG="eval <sequenced_genome.fas> <reference_genome.fas>"
  USG="$USG (<CDS.gff>|$NOGFF) <reads.fastq> [<mates.fastq>]"
  SGENOME=$2
  shift
fi
GENOME=$2
ANNOTATION=$3
READS=$4
MATES=$5
MAP="${READS}.${GENOME}"
MAPM="${MATES}.${GENOME}"

function show_globals {
  echo "USG=$USG"
  echo "ACTION=$ACTION"
  echo "SGENOME=$SGENOME"
  echo "GENOME=$GENOME"
  echo "ANNOTATION=$ANNOTATION"
  echo "READS=$READS"
  echo "MATES=$MATES"
  echo "MAP=$MAP"
  echo "MAPM=$MAPM"
}

function f_help_header {
  echo "Reference-based correction of homopolymer length in sequencing reads."
  echo
  echo "  Usage: $0 $USG"
  echo
}

function f_help_footer {
  echo
  echo "Important: this script only works correctly if all files are located"
  echo "           in the working directory where the script is called"
  echo "           (you can use symbolic links, of course)"
  echo
  echo "Further parameters can be optionally passed using env variables:"
  echo " - BWA          path to bwa binary             (default: $BWA)"
  echo " - SAMTOOLS     path to samtools binary        (default: $SAMTOOLS)"
  echo " - GT           path to GenomeTools binary     (default: $GT)"
  echo " - BWA_PARAMS   parameters for bwa aln         (default: $BWA_PARAMS)"
  echo " - HOP_PARAMS   parameters for hop             (default: $HOP_PARAMS)"
  echo " - USESW        set to \"TRUE\" to use bwa bwasw instead of bwa aln"
  echo " - SHOWONLY     set to \"TRUE\" to echo commands, without executing"
}

function f_shorthelp {
  f_help_header
  echo "Use \"$0 help\" for more information."
}

function f_help {
  f_help_header
  echo " - <action> is one of the following:"
  echo
  echo "     prepare       prepare the data necessary for the correction"
  echo "                   (index preparation, mapping, sorting, etc)"
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
  echo "                   (DISCLAIMER: this may have unintended effects, as"
  echo "                   it deletes files solely based on their filename)"
  echo
  echo "     run           prepare, correct and clean"
  echo
  echo "     stats         collect statistics for each correction position"
  echo "                   (useful for debugging and evaluation,"
  echo "                   it may be slow and require lot of memory)"
  echo
  echo "     eval-help     more information about evaluation actions;"
  echo
  echo "     help          shows this message"
  echo
  echo " - <genome.fas> is a single sequence Fasta file which"
  echo "   contains the reference against which the correction is done"
  echo
  echo " - CDS.gff contains the annotation of <genome.fas>;"
  echo "   instead of a filename '$NOGFF' can be specified; in this case,"
  echo "   homopolymer are corrected in the whole sequence instead of only "
  echo "   in coding regions"
  echo
  echo " - <reads.fastq> contains the reads to correct"
  echo
  echo " - <reads2.fastq> contains the mate pairs (if available)"
  f_help_footer
}

function f_eval_help {
  f_help_header
  echo " - <action> is one of the following:"
  echo
  echo "     eval              evaluate corrections to reads of a known genome"
  echo
  echo "     eval-make-true    prepare a set of \"true\" corrections"
  echo "                       for reads which have been sequenced from "
  echo "                       <sequenced_genome.fas>"
  echo
  echo "     eval-prepare      prepare for the evaluation of corrections "
  echo "                       against <reference_genome.fas>"
  echo
  echo "     eval-help         show this message"
  echo
  echo " - <sequenced_genome.fas> is a single sequence Fasta file which"
  echo "   contains the genome that has been sequenced"
  echo
  echo " - <reference_genome.fas> is a single sequence Fasta file which"
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
  echo "Before using eval, both eval-make-true and eval-prepare must be called."
  f_help_footer
}

function f_die {
  f_shorthelp
  exit $E_BADARGS
}

function f_index {
  echo "==== Create bwa genome index..."
  CMD="$BWA index $GENOME"
  echo
  echo $CMD
  if [ "$SHOWONLY" != "TRUE" ]; then $CMD; fi
  echo
}

function f_map {
  if [ "$USESW" == "TRUE" ]; then
    echo "==== Map reads to the genome using bwasw..."
    echo
    CMD="$BWA bwasw ${BWA_PARAMS} $GENOME $READS"
    echo $CMD '>' $MAP.sam
    if [ "$SHOWONLY" != "TRUE" ]; then $CMD > $MAP.sam; fi
  else
    echo "==== Map reads to the genome using bwa..."
    echo
    CMD="$BWA aln ${BWA_PARAMS} $GENOME $READS"
    echo $CMD '>' $MAP.sai
    if [ "$SHOWONLY" != "TRUE" ]; then $CMD > $MAP.sai; fi
    if [ "$MATES" != "" ]; then
      CMD1="$BWA aln ${BWA_PARAMS} $GENOME $MATES"
      CMD2="$BWA sampe $GENOME $MAP.sai $MAPM.sai $READS $MATES"
      echo $CMD1 '>' $MAPM.sai
      echo $CMD2 '>' $MAP.sam
      if [ "$SHOWONLY" != "TRUE" ]; then
        $CMD1 > $MAPM.sai
        $CMD2 > $MAP.sam
      fi
    else
      CMD="$BWA samse $GENOME $MAP.sai $READS"
      echo $CMD '>' $MAP.sam
      if [ "$SHOWONLY" != "TRUE" ]; then $CMD > $MAP.sam; fi
    fi
  fi
  echo
}

function f_sortmap {
  echo "==== Sort the mapping results using samtools..."
  echo
  CMD1="$SAMTOOLS view -Shu $MAP.sam"
  CMD2="$SAMTOOLS sort $MAP.bam sorted.$MAP"
  echo $CMD1 '>' $MAP.bam
  echo $CMD2
  if [ "$SHOWONLY" != "TRUE" ]; then
    $CMD1 > $MAP.bam
    $CMD2
  fi
  echo
}

function f_sortgff {
  echo "==== Process and sort the annotation..."
  echo
  awk '$3 == "CDS" \
       {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t." }' \
      $ANNOTATION > tmp_1.$ANNOTATION
  echo "##gff-version 3" > tmp_2.$ANNOTATION
  sort -k4 -n tmp_1.$ANNOTATION >> tmp_2.$ANNOTATION
  $GT gff3 -tidy -force -o sorted.$ANNOTATION tmp_2.$ANNOTATION
  rm tmp_1.$ANNOTATION
  rm tmp_2.$ANNOTATION
  echo
}

function f_encode {
  echo
  echo "==== Encode genome in GtEncseq format..."
  echo
  CMD="$GT encseq encode -v $GENOME"
  echo $CMD
  if [ "$SHOWONLY" != "TRUE" ]; then $CMD; fi
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
    HOP_PARAMS="-ann sorted.$ANNOTATION $HOP_PARAMS"
  fi
  if [ "$MATES" != "" -o "$USESW" == "TRUE" ]; then
    HOP_PARAMS="-reads $READS $MATES $HOP_PARAMS"
  else
    HOP_PARAMS="-o hop_$READS $HOP_PARAMS"
  fi
  CMD="$GT hop -v -ref $GENOME -map sorted.$MAP.bam $HOP_PARAMS"
  echo $CMD
  if [ "$SHOWONLY" != "TRUE" ]; then $CMD; fi
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

# position of fields in correction stats output table
declare -A FIELD
FIELD[R_HPOS]='$1'
FIELD[EDIT]='$2'
FIELD[S_HPOS]='$3'
FIELD[S_HEND]='$4'
FIELD[S_CHAR]='$5'
FIELD[S_OR]='$6'
FIELD[C_LEN]='$7'
FIELD[COVERAGE]='$8'
FIELD[R_HLEN]='$9'
FIELD[R_SUPP]='$10'
FIELD[S_HLEN]='$11'
FIELD[A_HLEN]='$12'
FIELD[A_SUPP]='$13'
FIELD[S_MAPQ]='$14'
FIELD[S_Q_BEF]='$15'
FIELD[S_Q_FIRST]='$16'
FIELD[S_Q_MIN]='$17'
FIELD[S_Q_AVE]='$18'
FIELD[S_Q_MAX]='$19'
FIELD[S_Q_RANGE]='$20'
FIELD[S_Q_LAST]='$21'
FIELD[S_Q_AFT]='$22'
FIELD[S_QUAL]='$23'
FIELD[S_ID]='$24'
FIELD[C_CLASS]='$25'

# combined fields
FIELD[EDIT_OP]="${FIELD[EDIT]}\"-\"${FIELD[S_CHAR]}"
FIELD[MULTIEDIT_OP]="${FIELD[EDIT]}\"-\"${FIELD[S_CHAR]}\"-\"${FIELD[C_LEN]}"
FIELD[A_HLEN_DIFF]="${FIELD[A_HLEN]}-${FIELD[R_HLEN]}"
FIELD[S_HLEN_DIFF]="${FIELD[S_HLEN]}-${FIELD[R_HLEN]}"
FIELD[NON_C_LEN]="(${FIELD[S_HLEN]}>${FIELD[R_HLEN]}?${FIELD[S_HLEN]}-${FIELD[R_HLEN]}:${FIELD[R_HLEN]}-${FIELD[S_HLEN]})-${FIELD[C_LEN]}"
FIELD[A_SUPP_DIFF]="${FIELD[A_SUPP]}-${FIELD[R_SUPP]}"

function f_pos_stats {
  count_all () { RETVAL=`grep '^'$2 $1 | wc -l`; }
  count_all_pos () { RETVAL=`awk '/^'$2'/ && ('${FIELD[S_HEND]}' == "'$3'")' $1 \
    | wc -l`; }
  count () { RETVAL=`awk '/^'$2'/ && ('${FIELD[S_CHAR]}' == "'$3'")' $1 | wc -l`; }
  count_pos () { RETVAL=`awk '/^'$2'/ && ('${FIELD[S_HEND]}' == "'$3'") \
    && ('${FIELD[S_CHAR]}' == "'$4'")' $1 | wc -l`; }
  nof_seq_long () { RETVAL=`grep "^$2--$2" $1 | awk '{ print $2 }'`
    if [ "$RETVAL" == "" ]; then RETVAL=0; fi; }
  OUTFILE=$2
  MAX_READ_LEN=`grep 'maximum length' $3 | awk '{ print $4 }'`
  NOFSEQ=`grep 'sequences' $3 | awk '{ print $2 }'`
  nof_seq_long $3 1; S=$RETVAL
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
  count_all $1 I; I=$RETVAL
  count_all $1 D; D=$RETVAL
  A=$[$I+ $D]
  count $1 I a; Ia=$RETVAL
  count $1 I c; Ic=$RETVAL
  count $1 I g; Ig=$RETVAL
  count $1 I t; It=$RETVAL
  count $1 D a; Da=$RETVAL
  count $1 D c; Dc=$RETVAL
  count $1 D g; Dg=$RETVAL
  count $1 D t; Dt=$RETVAL
  echo -e "all\t100.00\t$A\t$I\t$D\t$Ia\t$Ic\t$Ig\t$It\t$Da\t$Dc\t$Dg\t$Dt" \
    >> $OUTFILE
  for ((i = 2; i < $MAX_READ_LEN; i++)); do
    nof_seq_long $3 $i
    S=$[$RETVAL+$S]
    P=`echo "(1.0 - ($S / $NOFSEQ)) * 100" | bc -l`
    count_all_pos $1 I $i; I=$RETVAL
    count_all_pos $1 D $i; D=$RETVAL
    A=$[$I+ $D]
    count_pos $1 I $i a; Ia=$RETVAL
    count_pos $1 D $i a; Da=$RETVAL
    count_pos $1 I $i c; Ic=$RETVAL
    count_pos $1 D $i c; Dc=$RETVAL
    count_pos $1 I $i g; Ig=$RETVAL
    count_pos $1 D $i g; Dg=$RETVAL
    count_pos $1 I $i t; It=$RETVAL
    count_pos $1 D $i t; Dt=$RETVAL
    printf "$i\t%.2f\t$A\t$I\t$D\t$Ia\t$Ic\t$Ig\t$It\t$Da\t$Dc\t$Dg\t$Dt\n" $P \
      >> $OUTFILE
  done
}

function f_edits_per_read {
  grep -v '^#' $1 | awk '{ print '${FIELD[S_ID]}' }' | sort | uniq -c > $MAP~tmpfile
  echo "# edited reads: `cat $MAP~tmpfile | wc -l`" >> $1
  t=0
  for ((i=1;i<10;i++)); do
    c=`awk ' $1 == '$i $MAP~tmpfile | wc -l`
    if [ $c == 0  ]; then
      break
    else
      t=$[$t+($c*$i)]
      echo -e "# reads edited $i time(s): $c" >> $1
    fi
  done
  echo -e "# total edits: $t" >> $1
  rm $MAP~tmpfile
}

function f_stats {
  echo "==== Collecting correction statistics..."
  if [ "$ANNOTATION" != "$NOGFF" ]; then
    HOP_PARAMS="-ann sorted.$ANNOTATION $HOP_PARAMS"
  fi
  if [ "$MAKETRUE" == "TRUE" ]; then
    HOP_PARAMS="-state-of-truth $HOP_PARAMS"
  else
    HOP_PARAMS="-stats $HOP_PARAMS"
  fi
  if [ "$MATES" != "" -o "$USESW" == "TRUE" ]; then
    HOP_PARAMS="-reads $READS $MATES $HOP_PARAMS"
  else
    HOP_PARAMS="-o hop_$READS $HOP_PARAMS"
  fi
  CMD="$GT hop -v -ref $GENOME -map sorted.$MAP.bam $HOP_PARAMS"
  echo "$CMD" '>' $MAP.hop_stats
  if [ "$SHOWONLY" != "TRUE" ]; then
    $CMD > $MAP.hop_stats
    echo
    f_edits_per_read $MAP.hop_stats
  fi
  echo "==== done "
  echo
  echo "correction statistics:                $MAP.hop_stats"
  echo
}

stats2eds () {
  grep -v '^#' $MAP.hop_stats | awk '{ printf("%-30s\t%-5s\t%s\t%s\n", \
      '${FIELD[S_ID]}', '${FIELD[S_HEND]}', \
      '${FIELD[EDIT]}', '${FIELD[C_LEN]}') }' | sort > $MAP.eds
}

function f_eval_make_true {
  GENOME=$SGENOME
  MAP="${READS}.${GENOME}"
  MAPM="${MATES}.${GENOME}"
  ANNOTATION=$NOGFF
  if [ "$NOMAP" != "TRUE" ]; then f_prepare; fi
  MAKETRUE="TRUE"
  f_stats
  if [ "$SHOWONLY" != "TRUE" ]; then stats2eds; fi
}

function f_eval_prepare {
  f_prepare
}

function f_eval_mark {
  echo "==== Mark hop_stats lines using evaluation results"
  SMAP="${READS}.${SGENOME}"
  diff $SMAP.eds $MAP.eds | grep '^>' | sort > $MAP.fp
  ruby -e "
    fp = IO.read('$MAP.fp').split(%Q[\n])
    S_ID = %Q[${FIELD[S_ID]}][1..-1].to_i - 1
    S_HEND = %Q[${FIELD[S_HEND]}][1..-1].to_i - 1
    fp = fp.map{|x| x = x.split; %Q[#{x[1]} #{x[2]}]}
    infile = File.new('$MAP.hop_stats')
    outfile = File.new('$MAP.hop_stats_eval', 'w')
    infile.each do |line|
      line.chomp!
      if line =~ /^# coordinates/
        outfile.puts %Q[# c_class = classification of correction ]+
          %Q[(TP=true positive; FP=false positive)]
        outfile.puts line
      elsif line =~ /^# edit\sr_hpos/
        outfile.puts %Q[#{line}\tc_class]
      elsif line =~ /^#/
        outfile.puts line
      else
        splitted = line.split
        if (fp.include?(%Q[#{splitted[S_ID]} #{splitted[S_HEND]}]))
          outfile.puts %Q[#{line}\tFP]
        else
          outfile.puts %Q[#{line}\tTP]
        end
      end
    end
  "
  echo "marked hop_stats:   $MAP.hop_stats_eval"
}

function f_eval_proc_marked {
  awk ${FIELD[C_CLASS]}' == "FP"' $MAP.hop_stats_eval > $MAP.hop_stats_FP
  awk ${FIELD[C_CLASS]}' == "TP"' $MAP.hop_stats_eval > $MAP.hop_stats_TP
  echo "false positive only stats: $MAP.hop_stats_FP"
  echo "true positive only stats: $MAP.hop_stats_TP"
  QUANT_MEASURES="S_HPOS S_HEND C_LEN COVERAGE R_HLEN R_SUPP S_HLEN"
  QUANT_MEASURES="${QUANT_MEASURES} A_HLEN A_SUPP S_MAPQ S_Q_BEF S_Q_FIRST"
  QUANT_MEASURES="${QUANT_MEASURES} S_Q_MIN S_Q_AVE S_Q_MAX S_Q_RANGE S_Q_LAST"
  QUANT_MEASURES="${QUANT_MEASURES} S_Q_AFT A_HLEN_DIFF S_HLEN_DIFF NON_C_LEN"
  QUANT_MEASURES="${QUANT_MEASURES} A_SUPP_DIFF"
  for MEASURE in ${QUANT_MEASURES}; do
    awk '{print int('${FIELD[$MEASURE]}')}' $MAP.hop_stats_FP \
      | sort -n | uniq -c > $MAP.hop_FP_${MEASURE}_distri
    awk '{print int('${FIELD[$MEASURE]}')}' $MAP.hop_stats_TP \
      | sort -n | uniq -c > $MAP.hop_TP_${MEASURE}_distri
  done
  QUAL_MEASURES="S_CHAR EDIT S_OR EDIT_OP MULTIEDIT_OP"
  for MEASURE in ${QUAL_MEASURES}; do
    awk '{print '${FIELD[$MEASURE]}'}' $MAP.hop_stats_FP \
      | sort | uniq -c > $MAP.hop_FP_${MEASURE}_distri
    awk '{print '${FIELD[$MEASURE]}'}' $MAP.hop_stats_TP \
      | sort | uniq -c > $MAP.hop_TP_${MEASURE}_distri
  done
  MEASURES="${QUAL_MEASURES} ${QUANT_MEASURES}"
  for MEASURE in $MEASURES; do
    ruby -e "
      f = File.open('$MAP.hop_${MEASURE}_distri', 'w')
      tp={}
      fp={}
      IO.read(%Q[$MAP.hop_TP_${MEASURE}_distri]).split(%Q[\n]).
         each{|x| x=x.split; tp[x[1]]=x[0]; fp[x[1]]=0}
      IO.read(%Q[$MAP.hop_FP_${MEASURE}_distri]).split(%Q[\n]).
         each{|x| x=x.split; fp[x[1]]=x[0]; tp[x[1]]||=0}
      f.printf(%Q[${MEASURE}\tALL\tTP\tFP\tTPR\tFPR\n])
      keys = tp.keys
      if (keys[0].to_i.to_s == keys[0])
        keys = keys.sort_by {|x| x.to_i}
      else
        keys.sort!
      end
      keys.each do |k|
        f.printf(%Q[#{k}\t%s\t%s\t%s\t%.2f\t%.2f\n],
          tp[k].to_i + fp[k].to_i, tp[k], fp[k],
          tp[k].to_f * 100 / (tp[k].to_f+fp[k].to_f),
          fp[k].to_f * 100 / (tp[k].to_f+fp[k].to_f))
      end
      f.close
    "
    echo "${MEASURE} distribution: $MAP.hop_${MEASURE}_distri"
  done
}

function f_eval_stats {
  SMAP="${READS}.${SGENOME}"
  EDSDIFF="$SMAP.$GENOME.eds_diff"
  INFO="$SMAP.$GENOME.eval"
  diff $SMAP.eds $MAP.eds > $EDSDIFF
  CORR=`cat $MAP.eds | wc -l`
  ERR=`cat $SMAP.eds | wc -l`
  FN=`grep '^<' $EDSDIFF | wc -l`
  FP=`grep '^>' $EDSDIFF | wc -l`
  TP=$[$CORR-$FP]
  GENOMEDESC=`head -n 1 $GENOME`
  SGENOMEDESC=`head -n 1 $SGENOME`
  echo "==== Computing evaluation statistics..."
  echo | tee $INFO
  echo "Statistics of the correction of:" | tee -a $INFO
  if [ "$MATES" == "" ];then
    echo "- reads:                      $READS" | tee -a $INFO
  else
    echo "- reads:                      $READS + $MATES" | tee -a $INFO
  fi
  echo "- state of truth based upon: $SGENOMEDESC" | tee -a $INFO
  echo | tee -a $INFO
  MIN_HLEN=`grep -P -o "(?<=Distribution of homopolymers of length >= )\d+" \
    $MAP.hop_stats`
  if [ "$MIN_HLEN" == "" ]; then MIN_HLEN="n.a."; fi
  NOF_READS=`grep -P -o "(?<=segments in SAM file:)\s+\d+" \
    $MAP.hop_stats | grep -P -o "\d+"`
  NOF_NMAP=`grep -P -o "(?<=not mapping:)\s+\d+" \
    $MAP.hop_stats | grep -P -o "\d+"`
  if [ "$NOF_READS" != "" -a "$NOF_NMAP" != "" ]; then
    NOF_MAP=$[$NOF_READS - $NOF_NMAP]
  else
    NOF_MAP="n.a."
  fi
  if [ "$NOF_READS" == "" ]; then NOF_READS="n.a."; fi
  if [ "$NOF_NMAP" == "" ]; then NOF_NMAP="n.a."; fi
  echo "Correction parameters:"      | tee -a $INFO
  echo "- reference for correction:   $GENOMEDESC" | tee -a $INFO
  if [ "$ANNOTATION" != "$NOGFF" ]; then
    echo "- annotation:                 $ANNOTATION" | tee -a $INFO
  else
    echo "- annotation:                 none" | tee -a $INFO
  fi
  echo "- minimal homopol. length:    ${MIN_HLEN}" | tee -a $INFO
  echo | tee -a $INFO
  echo "Mapping results:" | tee -a $INFO
  echo "- total reads:        $NOF_READS" | tee -a $INFO
  echo "- mapping:            $NOF_MAP"   | tee -a $INFO
  echo "- not mapping:        $NOF_NMAP"  | tee -a $INFO
  echo | tee -a $INFO
  echo "Correction results:"         | tee -a $INFO
  echo "- homopol. errors:    $ERR"  | tee -a $INFO
  echo "- corrected homopol.: $CORR" | tee -a $INFO
  echo "- true positives:     $TP"   | tee -a $INFO
  echo "- false positives:    $FP"   | tee -a $INFO
  echo "- false negatives:    $FN"   | tee -a $INFO
  SN=`echo "$TP * 100 / ($TP + $FN)" | bc -l`
  PR=`echo "$TP * 100 / ($TP + $FP)" | bc -l`
  printf -- "- sensitivity:        %.2f %%\n" $SN | tee -a $INFO
  printf -- "- precision:          %.2f %%\n" $PR | tee -a $INFO
  echo | tee -a $INFO
  echo "==== done "
  echo
  echo "evaluation statistics:                $INFO"
  echo
}

function f_eval {
  f_stats
  if [ "$SHOWONLY" != "TRUE" ]; then
    stats2eds
    f_eval_stats
  fi
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

if [ "$ACTION" != "help" -a "$ACTION" != "eval-help" ]; then
  if [ $# -lt 4 -o $# -gt 5 ]; then f_die; fi
fi
case "$ACTION" in
  'prepare')            f_prepare           ;;
  'correct')            f_correct           ;;
  'clean')              f_clean             ;;
  'run')                f_run               ;;
  'stats')              f_stats             ;;
  'dists')              f_dists             ;;
  'help')               f_help              ;;
  'eval-help')          f_eval_help         ;;
  'eval-make-true')     f_eval_make_true    ;;
  'eval-prepare')       f_eval_prepare      ;;
  'eval')               f_eval              ;;
  # parts of other actions:
  'index')              f_index             ;;
  'map')                f_map               ;;
  'sortmap')            f_sortmap           ;;
  'encode')             f_encode            ;;
  'sortgff')            f_sortgff           ;;
  'eval-stats')         f_eval_stats        ;;
  'eval-mark')          f_eval_mark         ;;
  'eval-proc-marked')   f_eval_proc_marked  ;;
  *)                    f_die               ;;
esac
