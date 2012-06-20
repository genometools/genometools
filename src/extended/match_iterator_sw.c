/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include "core/encseq.h"
#include "core/ma.h"
#include "core/seq.h"
#include "core/score_function.h"
#include "core/unused_api.h"
#include "extended/alignment.h"
#include "extended/match.h"
#include "extended/match_sw.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_sw.h"
#include "extended/match_iterator_rep.h"
#include "extended/swalign.h"

const GtMatchIteratorClass* gt_match_iterator_sw_class(void);

#define gt_match_iterator_sw_cast(M)\
        gt_match_iterator_cast(gt_match_iterator_sw_class(), M)

typedef struct {
  GtScoreFunction *sf;
  GtEncseq *es1, *es2;
  unsigned long seqno_es1,
                seqno_es2,
                min_len,
                max_edist;
  bool firstali;
} GtMatchIteratorSWMembers;

struct GtMatchIteratorSW {
  const GtMatchIterator parent_instance;
  GtMatchIteratorSWMembers *pvt;
};

static GtMatchIteratorStatus gt_match_iterator_sw_next(GtMatchIterator *mi,
                                                      GT_UNUSED GtMatch **match,
                                                      GT_UNUSED GtError *err)
{
  GtMatchIteratorSW *mis;
  GtSeq *seq_a, *seq_b;
  char *a, *b;
  const char *adesc, *bdesc;
  GtAlignment *ali = NULL;
  unsigned long seqlen_a, seqlen_b, seqpos;
  GtRange arng, brng;
  gt_assert(mi && match);

  mis = gt_match_iterator_sw_cast(mi);
  while (true) {
    if (!mis->pvt->firstali)
      mis->pvt->seqno_es2++;
    if (mis->pvt->seqno_es2 == gt_encseq_num_of_sequences(mis->pvt->es2)) {
      mis->pvt->seqno_es1++;
      if (mis->pvt->seqno_es1 == gt_encseq_num_of_sequences(mis->pvt->es1))
        return GT_MATCHER_STATUS_END;
      mis->pvt->seqno_es2 = 0;
    }
    seqlen_a = gt_encseq_seqlength(mis->pvt->es1, mis->pvt->seqno_es1);
    seqlen_b = gt_encseq_seqlength(mis->pvt->es2, mis->pvt->seqno_es2);
    /* XXX: reuse buffers for performance improvement */
    a = gt_malloc(seqlen_a * sizeof (char));
    seqpos = gt_encseq_seqstartpos(mis->pvt->es1, mis->pvt->seqno_es1);
    gt_encseq_extract_decoded(mis->pvt->es1, a, seqpos, seqpos + seqlen_a - 1);
    b = gt_malloc(seqlen_b * sizeof (char));
    seqpos = gt_encseq_seqstartpos(mis->pvt->es2, mis->pvt->seqno_es2);
    gt_encseq_extract_decoded(mis->pvt->es1, b, seqpos, seqpos + seqlen_b - 1);
    seq_a = gt_seq_new(a, seqlen_a, gt_encseq_alphabet(mis->pvt->es1));
    seq_b = gt_seq_new(b, seqlen_b, gt_encseq_alphabet(mis->pvt->es2));
    ali = gt_swalign(seq_a, seq_b, mis->pvt->sf);
    mis->pvt->firstali = false;
    if (ali && gt_alignment_get_length(ali) >= mis->pvt->min_len
          && gt_alignment_eval(ali) <= mis->pvt->max_edist) {
      break;
    }
    gt_alignment_delete(ali);
    gt_seq_delete(seq_a);
    gt_seq_delete(seq_b);
    gt_free(a);
    gt_free(b);
  }
  arng = gt_alignment_get_urange(ali);
  brng = gt_alignment_get_vrange(ali);
  adesc = gt_encseq_description(mis->pvt->es1, &seqlen_a, mis->pvt->seqno_es1);
  bdesc = gt_encseq_description(mis->pvt->es2, &seqlen_b, mis->pvt->seqno_es2);
  *match = gt_match_sw_new("", "",
                           mis->pvt->seqno_es1,
                           mis->pvt->seqno_es2,
                           gt_alignment_get_length(ali),
                           gt_alignment_eval(ali),
                           arng.start, brng.start,
                           arng.end, brng.end,
                           GT_MATCH_DIRECT);
  gt_match_set_seqid1_nt(*match, adesc, seqlen_a);
  gt_match_set_seqid2_nt(*match, bdesc, seqlen_b);
  gt_alignment_delete(ali);
  gt_seq_delete(seq_a);
  gt_seq_delete(seq_b);
  gt_free(a);
  gt_free(b);
  return GT_MATCHER_STATUS_OK;
}

GtMatchIterator* gt_match_iterator_sw_new(GtEncseq *es1, GtEncseq *es2,
                                          GtScoreFunction *sf,
                                          unsigned long min_len,
                                          unsigned long max_edist,
                                          GT_UNUSED GtError *err)
{
  GtMatchIterator *mi;
  GtMatchIteratorSW *mis;
  gt_assert(es1 && es2 && sf);
  gt_error_check(err);

  mi = gt_match_iterator_create(gt_match_iterator_sw_class());
  mis = (GtMatchIteratorSW*) mi;
  mis->pvt = gt_calloc(1, sizeof (GtMatchIteratorSWMembers));
  mis->pvt->es1 = gt_encseq_ref(es1);
  mis->pvt->es2 = gt_encseq_ref(es2);
  mis->pvt->sf = gt_score_function_ref(sf);
  mis->pvt->min_len = min_len;
  mis->pvt->max_edist = max_edist;
  mis->pvt->firstali = true;
  return mi;
}

static void gt_match_iterator_sw_free(GtMatchIterator *mi)
{
  GtMatchIteratorSW *mis;
  if (!mi) return;
  mis = gt_match_iterator_sw_cast(mi);
  gt_encseq_delete(mis->pvt->es1);
  gt_encseq_delete(mis->pvt->es2);
  gt_score_function_delete(mis->pvt->sf);
  gt_free(mis->pvt);
}

const GtMatchIteratorClass* gt_match_iterator_sw_class(void)
{
  static const GtMatchIteratorClass *mic;

  if (!mic) {
    mic = gt_match_iterator_class_new(sizeof (GtMatchIteratorSW),
                                      gt_match_iterator_sw_next,
                                      gt_match_iterator_sw_free);
  }
  return mic;
}
