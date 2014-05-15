/*
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c) 2011-2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2014      Genome Research Ltd.

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
  GtArray *orf_start[3],
          *starts_to_output;
  bool found_start[3],
       find_all;
};

GtORFIterator* gt_orf_iterator_new(GtCodonIterator *ci,
                                   GtTranslator *translator)
{
  GtORFIterator *orfi;
  gt_assert(ci && translator);
  orfi = gt_malloc(sizeof (GtORFIterator));
  int i;

  orfi->ci = ci;
  orfi->translator = translator;
  orfi->find_all = false;
  orfi->starts_to_output = NULL;
  for (i = 0; i < 3; i++) {
    orfi->orf_start[i] = gt_array_new(sizeof (GtUword));
    orfi->found_start[i] = false;
  }

  return orfi;
}

void gt_orf_iterator_find_all(GtORFIterator *orfi)
{
  gt_assert(orfi);
  orfi->find_all = true;
}

void gt_orf_iterator_find_longest(GtORFIterator *orfi)
{
  gt_assert(orfi);
  orfi->find_all = false;
}

GtORFIteratorStatus gt_orf_iterator_next(GtORFIterator *orfi, GtRange *orf_rng,
                                         unsigned int *orf_frame, GtError *err)
{
  unsigned int frame;
  char translated;
  GtTranslatorStatus state;
  gt_assert(orfi);

  if (orfi->starts_to_output != NULL) {
    if (gt_array_size(orfi->starts_to_output) > 0) {
      GtUword start = *(GtUword*) gt_array_pop(orfi->starts_to_output);
      orf_rng->start = start;
      /* ORF ends before GT_STOP_AMINO */
      orf_rng->end = gt_codon_iterator_current_position(orfi->ci) - 2;
      return GT_ORF_ITERATOR_OK;
    } else {
      orfi->starts_to_output = NULL;
    }
  }

  while ((state = gt_translator_next(orfi->translator,
                            &translated, &frame, err)) != GT_TRANSLATOR_END) {
    if (state == GT_TRANSLATOR_ERROR)
      return GT_ORF_ITERATOR_ERROR;
    if (translated == GT_START_AMINO
          && (orfi->find_all || !orfi->found_start[frame])) {
      GtUword startv = gt_codon_iterator_current_position(orfi->ci) - 1;
      gt_array_add(orfi->orf_start[frame], startv);
      orfi->found_start[frame] = true;
    }
    else if ((translated == GT_STOP_AMINO) &&
             (orfi->found_start[frame] == true)) {
      GtUword start = *(GtUword*) gt_array_pop(orfi->orf_start[frame]);
      orf_rng->start = start;
      /* ORF ends before GT_STOP_AMINO */
      orf_rng->end = gt_codon_iterator_current_position(orfi->ci) - 2;
      *orf_frame = frame;
      orfi->found_start[frame] = false;
      if (gt_array_size(orfi->orf_start[frame]) > 0) {
        orfi->starts_to_output = orfi->orf_start[frame];
      }
      return GT_ORF_ITERATOR_OK;
    }
  }
  return GT_ORF_ITERATOR_END;
}

void gt_orf_iterator_delete(GtORFIterator *ofi)
{
  int i;
  if (!ofi) return;
  for (i = 0; i < 3; i++) {
    gt_array_delete(ofi->orf_start[i]);
  }
  gt_free(ofi);
}
