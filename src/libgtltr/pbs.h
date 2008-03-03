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
#include "libgtcore/seq.h"
#include "libgtcore/strand.h"
#include "libgtcore/scorefunction.h"
#include "libgtext/alignment.h"
#include "libgtltr/repeattypes.h"
#include "libgtltr/ltrharvest-opt.h"

/* This struct holds information about a primer binding site (PBS). */
typedef struct PBS_Hit PBS_Hit;

struct PBS_Hit {
  unsigned long start,
                end,
                edist,
                offset,
                tstart,
                strand;
  double score;
  const char *trna;
  Alignment *ali;
};

/* Aligns tRNA from a library to the LTR retrotransposon candidate and
   returns highest-scoring hit (newly created). */
PBS_Hit* pbs_find(const char *seq,
                  LTRboundaries *line,
                  unsigned long seqlen,
                  unsigned long ltrlen,
                  Bioseq *trna_lib,
                  LTRharvestoptions *lo,
                  Error *err);

#endif
