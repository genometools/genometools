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
MIN_HLEN=${MIN_HLEN:=3}
BWA_T=${BWA_T:=4}
HOP_PARAMS=${HOP_PARAMS:=""}

if [ $# -lt 2 -o $# -gt 3 ]; then
  echo "Reference based homopolymer length correction in sequencing reads."
  echo
  echo "  Usage: $0 <reads.fastq> <reference.fasta> [<CDS.gff>]"
  echo
  echo " - reference.fasta must be a single fasta sequence which"
  echo "   defines the reference against which the correction is done"
  echo " - if a CDS.gff is provided, only position which map to a CDS"
  echo "   can be corrected (conservative mode, reccomended)"
  echo
  echo "Further parameters can be optionally passed using env variables:"
  echo " - BWA          path to bwa binary           (default: $BWA)"
  echo " - GT           path to GenomeTools binary   (default: $GT)"
  echo " - SAMTOOLS     path to samtools binary      (default: $SAMTOOLS)"
  echo " - MIN_HLEN     minimal homopolymer length   (default: $MIN_HLEN)"
  echo " - BWA_T        number of bwa threads        (default: $BWA_T)"
  echo " - HOP_PARAMS   additional hopcorrect params (default: $HOP_PARAMS)"
  exit $E_BADARGS
fi

# command line parameters:
READS=$1
GENOME=$2
ANNOTATION=$3

echo "==== Index the genome..."
echo
$BWA index $GENOME
echo
echo "==== Map reads to the genome..."
echo
$BWA aln -t ${BWA_T} $GENOME $READS > $READS.bwa
$BWA samse $GENOME $READS.bwa $READS > $READS.sam
echo
echo "==== Sort the mapping results..."
echo
$SAMTOOLS view -Shu $READS.sam > $READS.bam
$SAMTOOLS sort $READS.bam sorted.$READS
if [ "$ANNOTATION" != "" ]; then
  echo
  echo "==== Sort the annotation..."
  echo
  $GT gff3 -force -o sorted.$ANNOTATION -sort -retainids $ANNOTATION
  AOPT="-a sorted.$ANNOTATION"
else
  AOPT=""
fi
echo
echo "==== Encode genome in GtEncseq format..."
echo
$GT encseq encode -v $GENOME
echo
echo "==== Correct homopolymers..."
echo
$GT dev hopcorrect -v -r $GENOME -m sorted.$READS.bam $AOPT \
    -o hopcorrect_$READS $HOP_PARAMS
echo
echo "==== done "
echo
echo "the corrected reads are: hopcorrect_$READS"
# cleanup
rm -f $GENOME.amb $GENOME.ann $GENOME.bwt $GENOME.pac $GENOME.sa
rm -f $READS.bwa $READS.sam
rm -f $READS.bam sorted.$READS.bam
rm -f $GENOME.esq $GENOME.des $GENOME.md5 $GENOME.sds
rm -f sorted.$ANNOTATION
