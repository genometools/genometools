/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef PDOM_H
#define PDOM_H

#ifdef HAVE_HMMER

#include "core/array.h"
#include "core/hashmap.h"
#include "core/strand.h"
#include "core/str_array.h"
#include "ltr/ltrelement.h"
#include "structs.h"

typedef struct GtPdomOptions {
  double evalue_cutoff;
  GtStrArray *hmm_files;
  GtArray *plan7_ts;
  struct threshold_s thresh;
  unsigned int nof_threads,
               chain_max_gap_length;
} GtPdomOptions;

typedef struct GtPdomHit {
  struct tophit_s *hits_fwd, *hits_rev;
  GtStrand strand;
  GtArray *best_chain;
} GtPdomHit;

typedef struct GtPdomResults {
  GtHashmap *domains;
  double combined_e_value_fwd,
         combined_e_value_rev;
  bool empty;
} GtPdomResults;

int  gt_pdom_load_hmm_files(GtPdomOptions *opts, GtError *err);
void gt_pdom_convert_frame_position(GtRange *rng, int frame);
void gt_pdom_find(const char *seq, const char *rev_seq, GtLTRElement *element,
                  GtPdomResults *results, GtPdomOptions *opts);
void gt_pdom_clear_hmms(GtArray*);
void gt_pdom_clear_domain_hit(void*);

#endif
#endif
