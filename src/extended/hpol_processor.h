/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef HPOL_PROCESSOR_H
#define HPOL_PROCESSOR_H

#include "core/encseq.h"
#include "extended/seqpos_classifier.h"
#include "extended/aligned_segments_pile.h"
#include "core/seq_iterator_api.h"
#include "core/logger.h"

/* The <GtHpolProcessor> class implements an homopolymer processor.
   It scans a sequence and searches stretches of at least <hmin> identical
   symbols and processes them efficiently. Its main purpose is e.g. the
   correction of homopolymer errors in sequencing reads. */
typedef struct GtHpolProcessor GtHpolProcessor;

/* Creates a new <GtHpolProcessor> for <encseq>, recognising homopolymers of
   length at least <hmin>. */
GtHpolProcessor* gt_hpol_processor_new(GtEncseq *encseq, unsigned long hmin);

/* Start a scanning run over the sequence given at the instantiation of <hpp>.
   If <logger> is not NULL, then verbose information is directed there.
   Returns 0 on success, on error a negative number is returned and <err> is
   set accordingly. */
int              gt_hpol_processor_run(GtHpolProcessor *hpp, GtLogger *logger,
                                       GtError *err);

/* Enables the correction of homopolymers errors in the segments (usually
   sequencing reads) provided by <asp>.
   Some thresholds are defined for the correction to take place:
   the <read_hmin> parameter is the minimal length of the homopolymer in the
   reads; the <altmax> parameter (maximum alternative consensus) is a value
   between 0.0 and 1.0 and is the maximal portion of the reads with agree
   on a different length than the reference; the <refmin> parameter
   (minimal reference support) is a value between 0.0 and 1.0 and is the
   minimal portion of the reads which must agree with the reference;
   <mapqmin> is the minimal mapping quality; <covmin> is the minimal
   number of reads covering the homopolymer; <allow_partial> defines the
   behaviour when an insertion must take place and there are not enough gaps
   are available in the aligned read (true: insert at most number_of_gaps
   symbols; false: do not correct); <allow_multiple> defines the behaviour
   when a read requires multiple corrections (true: correct each time,
   false: correct only first homopolymer error found); <clenmax> defines
   the maximal correction length. */
void             gt_hpol_processor_enable_segments_hlen_adjustment(
                                                     GtHpolProcessor *hpp,
                                                     GtAlignedSegmentsPile *asp,
                                                     unsigned long read_hmin,
                                                     double altmax,
                                                     double refmin,
                                                     unsigned long mapqmin,
                                                     unsigned long covmin,
                                                     bool allow_partial,
                                                     bool allow_multiple,
                                                     unsigned long clenmax);

/* Restrict processing to homopolymers whose end position is classified
   by <spc> as inside a feature (e.g. inside a CDS). */
void             gt_hpol_processor_restrict_to_feature_type(
                                                       GtHpolProcessor *hpp,
                                                       GtSeqposClassifier *spc);

/* Make <hpp> output the segments sequence in FastQ format to <outfile> during
   run. This mainly is useful for the cases in which each read has exactly one
   alignment (e.g. output of bwa samse/sampe but not bwasw) - otherwise
   a single read is output several times. */
void             gt_hpol_processor_enable_direct_segments_output(
                                                           GtHpolProcessor *hpp,
                                                           GtFile *outfile);

/* Make <hpp> output the segments sorted by the order in which sequences
   are returned by the <reads_iters>. This is a more general output mode
   than the direct segment output, but it is slower and requires more memory
   as the information is stored during the run and output at the end.
   The <read_iters> are <nfiles> GtSeqIterator objects (each corresponding
   to an input file). The reads from the i-th GtSeqIterator will
   be output to the i-th element of <outfiles>. This method assumes that the
   sequence ID (first word of the description line) of each sequence is
   unique.
   XXX: use a single input iterator and gt_seq_iterator_fastq_get_file_index()*/
void             gt_hpol_processor_enable_sorted_segments_output(
                                                    GtHpolProcessor *hpp,
                                                    unsigned long nfiles,
                                                    GtSeqIterator **reads_iters,
                                                    GtFile **outfiles);

/* Output some statistics about each correction position to <outfile>.
   Data are output as a tab-separated table, one row per correction
   position. If <output_multihit_stats> is set, stats are output using
   information from multiple reference matches, not just the one being actually
   corrected.  */
void             gt_hpol_processor_enable_statistics_output(
                                                     GtHpolProcessor *hpp,
                                                     bool output_multihit_stats,
                                                     GtFile *outfile);

/* Enable debug mode which is useful to compare refregion of the segments
   with the reference sequence. */
void             gt_hpol_processor_enable_aligned_segments_refregionscheck(
                                                    GtHpolProcessor *hpp,
                                                    GtAlignedSegmentsPile *asp);

/* Deletes <hpp> and frees all associated memory. */
void             gt_hpol_processor_delete(GtHpolProcessor *hpp);

#endif
