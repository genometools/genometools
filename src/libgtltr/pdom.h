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
} PdomOptions;

typedef struct PdomHit {
  struct tophit_s *hits_fwd, *hits_rev;
  struct hit_s *best_hit;
} PdomHit;

typedef struct PdomResults {
  Hashtable *domains;                  /* maps models to PdomHits */
  double combined_e_value;
  bool empty;
} PdomResults;

int  load_hmm_files(StrArray *files, Array *models, Error *err);
int  pdom_domain_report_hits(void *key, void *value, UNUSED void *data,
                             UNUSED Error *err);
void pdom_convert_frame_position(Range *rng, int frame);
void pdom_find(const char *seq, const char *rev_seq, LTRElement *element,
               PdomResults *results, PdomOptions *opts);
void pdom_clear_hmms(Array*);
void pdom_clear_domain_hit(void*);

#endif
