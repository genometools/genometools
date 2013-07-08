/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include "core/xansi_api.h"
#include "extended/multieoplist.h"
#include "core/ma.h"
#include "core/chardef.h"

struct GtMultieoplist {
  GtArray *meoplist;
};

GtMultieoplist *gt_multieoplist_new(void)
{
  GtMultieoplist *eops;
  eops = gt_malloc(sizeof (GtMultieoplist));
  eops->meoplist = gt_array_new(sizeof (GtMultieop));
  return eops;
}

void gt_multieoplist_delete(GtMultieoplist *eops)
{
  gt_array_delete(eops->meoplist);
  gt_free(eops);
}

static void gt_multieoplist_add_eop(GtMultieoplist *eops, AlignmentEoptype type)
{
  GtMultieop meop, *meop_ptr;
  if (!gt_array_size(eops->meoplist)) {
    meop.type = type;
    meop.steps = 1;
    gt_array_add(eops->meoplist, meop);
  }
  else {
    meop_ptr = gt_array_get_last(eops->meoplist);
    if (meop_ptr->type == type)
      meop_ptr->steps++; /* XXX: check for overflow */
    else {
      meop.type = type;
      meop.steps = 1;
      gt_array_add(eops->meoplist, meop);
    }
  }
}

void gt_multieoplist_add_replacement(GtMultieoplist *eops)
{
  gt_multieoplist_add_eop(eops, Replacement);
}

void gt_multieoplist_add_insertion(GtMultieoplist *eops)
{
  gt_multieoplist_add_eop(eops, Insertion);
}

void gt_multieoplist_add_deletion(GtMultieoplist *eops)
{
  gt_multieoplist_add_eop(eops, Deletion);
}

void gt_multieoplist_add_mismatch(GtMultieoplist *eops)
{
  gt_multieoplist_add_eop(eops, Mismatch);
}

void gt_multieoplist_add_match(GtMultieoplist *eops)
{
  gt_multieoplist_add_eop(eops, Match);
}

void gt_multieoplist_reset(GtMultieoplist *eops)
{
  gt_array_reset(eops->meoplist);
}

void gt_multieoplist_remove_last(GtMultieoplist *eops)
{
  GtMultieop *meop_ptr;
  gt_assert(eops && gt_array_size(eops->meoplist));
  meop_ptr = gt_array_get_last(eops->meoplist);
  gt_assert(meop_ptr->steps);
  if (meop_ptr->steps == 1)
    (void) gt_array_pop(eops->meoplist);
  else
    meop_ptr->steps--;
}

unsigned long gt_multieoplist_get_repdel_length(GtMultieoplist *eops)
{
  unsigned long len = 0, i;
  GtMultieop meop;
  for (i = gt_array_size(eops->meoplist); i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(eops->meoplist, i-1);
    if (meop.type == Replacement || meop.type == Deletion)
      len += meop.steps;
  }
  return len;
}

unsigned long gt_multieoplist_get_repins_length(GtMultieoplist *eops)
{
  unsigned long len = 0, i;
  GtMultieop meop;
  for (i = gt_array_size(eops->meoplist); i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(eops->meoplist, i-1);
    if (meop.type == Replacement || meop.type == Insertion)
      len += meop.steps;
  }
  return len;
}

unsigned long gt_multieoplist_get_length(GtMultieoplist *eops)
{
  return(gt_array_size(eops->meoplist));
}

GtMultieop *gt_multieoplist_get_entry(GtMultieoplist *eops, unsigned long index)
{
  gt_assert(index < gt_array_size(eops->meoplist));
  return((GtMultieop *) gt_array_get(eops->meoplist, index));
}

void gt_multieoplist_show(GtMultieoplist *eops, FILE *fp)
{
  unsigned long i, size;
  GtMultieop meop;

  gt_assert(eops != NULL);

  size = gt_array_size(eops->meoplist);

  gt_xfputc('[', fp);
  for (i = size; i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(eops->meoplist, i-1);
    switch (meop.type) {
      case Mismatch:
      case Match:
      case Replacement:
        gt_xfputc('R', fp);
        break;
      case Insertion:
        gt_xfputc('I', fp);
        break;
      case Deletion:
        gt_xfputc('D', fp);
        break;
    }
    fprintf(fp, " %lu", meop.steps);
    if (i != 1UL) {
      gt_xfputc(',', fp);
    }
  }
  gt_xfputs("]\n", fp);
}
