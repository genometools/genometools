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

#include "libgtcore/dlist.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/strand.h"
#include "libgtcore/strarray.h"
#include "structs.h"

typedef struct PdomOptions {
  double evalue_cutoff;
  StrArray *hmm_files;
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

#endif
