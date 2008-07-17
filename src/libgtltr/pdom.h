/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#include "libgtcore/array.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/strand.h"
#include "libgtcore/strarray.h"
#include "libgtltr/ltrelement.h"
#include "structs.h"

typedef struct PdomOptions {
  double evalue_cutoff;
  StrArray *hmm_files;
  Array *plan7_ts;
  struct threshold_s thresh;
  unsigned int nof_threads,
               chain_max_gap_length;
} PdomOptions;

typedef struct PdomHit {
  struct tophit_s *hits_fwd, *hits_rev;
  Strand strand;
  Array *best_chain;
} PdomHit;

typedef struct PdomResults {
  Hashtable *domains;
  double combined_e_value_fwd,
         combined_e_value_rev;
  bool empty;
} PdomResults;

int  pdom_load_hmm_files(PdomOptions *opts, Error *err);
void pdom_convert_frame_position(Range *rng, int frame);
void pdom_find(const char *seq, const char *rev_seq, LTRElement *element,
               PdomResults *results, PdomOptions *opts);
void pdom_clear_hmms(Array*);
void pdom_clear_domain_hit(void*);

#endif
#endif
