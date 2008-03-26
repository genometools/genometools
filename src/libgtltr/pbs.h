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

#ifndef PBS_H
#define PBS_H

#include "libgtcore/undef.h"
#include "libgtcore/bioseq.h"
#include "libgtcore/dlist.h"
#include "libgtcore/seq.h"
#include "libgtcore/strand.h"
#include "libgtcore/scorefunction.h"
#include "libgtext/alignment.h"
#include "libgtltr/ltrelement.h"

typedef struct PBSOptions {
  unsigned int radius,
               max_edist;
  Range alilen,
        offsetlen,
        trnaoffsetlen;
  int ali_score_match,
      ali_score_mismatch,
      ali_score_insertion,
      ali_score_deletion;
  Bioseq *trna_lib;
} PBSOptions;

/* This struct holds information about a primer binding site (PBS). */
typedef struct PBS_Hit {
  unsigned long start,
                end,
                edist,
                offset,
                tstart,
                alilen;
  Strand strand;
  double score;
  const char *trna;
} PBS_Hit;

typedef struct PBSResults {
  Dlist *hits_fwd, *hits_rev;
  PBS_Hit *best_hit;
} PBSResults;

/* Aligns tRNA from a library to the LTR retrotransposon candidate and
   returns highest-scoring hit (newly created). */
void  pbs_find(const char *seq,
               const char *rev_seq,
               LTRElement *element,
               PBSResults *results,
               PBSOptions *lo,
               Error *err);

void  pbs_clear_results(PBSResults *results);

#endif
