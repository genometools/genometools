/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef CSA_VARIABLE_STRANDS_H
#define CSA_VARIABLE_STRANDS_H

#include "libgtext/consensus_sa.h"
#include "libgtext/csa_gene.h"

/*
  This module wraps the method to construct consensus spliced alignments
  from the ``consensus_sa'' module to handle splice forms on variable strands
  conveniently.

  That is, after the call to consensus_sa() the (collected) splice forms are
  postprocessed into CSAGenes representing genes on variable strands.
*/

/* Returns an array of CSAGenes. */
Array* csa_variable_strands(const void *set_of_sas, unsigned long number_of_sas,
                            size_t size_of_sa, GetGenomicRangeFunc,
                            GetStrandFunc, GetExonsFunc);

#endif
