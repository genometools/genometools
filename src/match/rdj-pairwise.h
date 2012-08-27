/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_PAIRWISE_H
#define RDJ_PAIRWISE_H

#include "match/rdj-spmproc.h"           /* GtSpmproc and GtSpmprocA */
#include "core/intbits.h"                /* GtBitsequence */
#include "core/encseq.h"                 /* GtEncseq */

typedef enum
{
  GT_OVLFIND_SPM,        /* containments are not determined;
                            suf-pref matches are determined,
                            including the one that are also containments */
  GT_OVLFIND_CNT,        /* containments are determined;
                            suf-pref matches are not determined */
  GT_OVLFIND_ALL,        /* containments are determined;
                            suf-pref matches are detemined,
                            including the one that are also containments */
  GT_OVLFIND_PROPER_SPM  /* containments are determined;
                            suf-pref matches are determined
                            only if there is no containment */
}
GtOvlfindMode;

/*
(1) general arguments

m:          mode, see GtOvlfindMode definition
encseq:     pointer to a GtEncseq
revcompl:   if true, the second half of the fasta file must contain
            the reverse complements of the sequences in the first half
show_progressbar: if true, a progressbar is displayed

use_kmp:    [gt_rdj_pairwise_exact]
            if true, the KMP-based algorithm is used,
            if false, the brute force algorithm is used

max_error:  [gt_rdj_pairwise_approx]
            the maximal error rate allowed for both
            suffix-prefix matches and containments

(2) suffix-prefix matches related

min_length: only suffix-prefix matches over this length are considered
find_nonmaximal: if true, suffix-prefix matches shorter than the maximal
                 one are also considered
proc:       matches are processed using this function;
            it must be set to NULL in gt_rdj_CNT mode
procdata:   generic data pointer passed to proc;
            it must be set to NULL in gt_rdj_CNT mode

(3) containment related

cntfilter:   if true, contained reads are eliminated from the
             search for further suffix-prefix matches once identified;
             it must be set to false in gt_rdj_SPM mode
cntreads_in: if not NULL, this will be used as the contained reads list,
             i.e. any bit already set will prevent matches involving the
             corresponding read; it must be NULL in gt_rdj_SPM mode;
             the list will never be deleted
cntreads_out: if not NULL, the address of the contained reads list
              will be written here; it must be NULL in gt_rdj_SPM mode;
              if this and cntreads_in are both NULL, the list is deleted
nofreads:    if not NULL, the number of direct reads is written here
             (which is also the number of elements of the cntreads list)

*/

void gt_rdj_pairwise_exact(GtOvlfindMode m, GtEncseq *encseq,
    bool revcompl, bool show_progressbar, bool use_kmp,
    unsigned long min_length, bool find_nonmaximal, GtSpmproc proc,
    void *procdata, bool cntfilter, GtBitsequence *cntreads_in,
    GtBitsequence **cntreads_out, unsigned long *nofreads);

void gt_rdj_pairwise_approx(GtOvlfindMode m,  GtEncseq *encseq, bool revcompl,
    bool show_progressbar, double max_error, unsigned long min_length,
    bool find_nonmaximal, GtSpmprocA proc, void* procdata, bool cntfilter,
    GtBitsequence *cntreads_in, GtBitsequence **cntreads_out,
    unsigned long *nofreads);

#endif
