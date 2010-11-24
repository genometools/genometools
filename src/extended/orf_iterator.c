/*
  Copyright (c) 2010      Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/ma.h"
#include "core/codon_iterator_api.h"
#include "core/translator.h"
#include "core/trans_table.h"
#include "core/strand_api.h"
#include "core/undef_api.h"
#include "extended/orf_iterator_api.h"

struct GtORFIterator {
  GtCodonIterator *ci;
  GtTranslator *translator;
  unsigned long orf_start[3];
  bool found_start[3];
};

GtORFIterator* gt_orf_iterator_new(GtCodonIterator *ci,
                                   GtTranslator *translator)
{
  gt_assert(ci && translator);
  GtORFIterator *orfi = gt_malloc(sizeof (GtORFIterator));
  int i;

  orfi->ci = ci;
  orfi->translator = translator;

  for (i = 0; i < 3; i++) {
    orfi->orf_start[i] = GT_UNDEF_ULONG;
    orfi->found_start[i] = false;
  }
  return orfi;
}

GtORFIteratorStatus gt_orf_iterator_next(
                        GtORFIterator *orfi, GtRange *orf_rng,
                        unsigned int *orf_frame, GtError *err)
{
  gt_assert(orfi);
  unsigned int frame;
  char translated;
  GtTranslatorStatus state;

  while ((state = gt_translator_next(orfi->translator,
                            &translated, &frame, err)) != GT_TRANSLATOR_END) {
    if (state == GT_TRANSLATOR_ERROR)
      return GT_ORF_ITERATOR_ERROR;
    if (translated == GT_START_AMINO && !orfi->found_start[frame]) {
      orfi->orf_start[frame] = gt_codon_iterator_current_position(orfi->ci) - 1;
      orfi->found_start[frame] = true;
    }
    else if ((translated == GT_STOP_AMINO) &&
             (orfi->found_start[frame] == true)) {
      orf_rng->start = orfi->orf_start[frame];
      /* ORF ends before GT_STOP_AMINO */
      orf_rng->end = gt_codon_iterator_current_position(orfi->ci) - 2;
      *orf_frame = frame;
      orfi->found_start[frame] = false;
      return GT_ORF_ITERATOR_OK;
    }
  }
  return GT_ORF_ITERATOR_END;
}

void gt_orf_iterator_delete(GtORFIterator *ofi)
{
  if (!ofi) return;
  gt_free(ofi);
}
