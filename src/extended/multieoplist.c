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
  GtUword refcount;
};

GtMultieoplist *gt_multieoplist_new(void)
{
  GtMultieoplist *multieops;
  multieops = gt_malloc(sizeof (GtMultieoplist));
  multieops->meoplist = gt_array_new(sizeof (GtMultieop));
  multieops->refcount = 0;
  return multieops;
}

GtMultieoplist *gt_multieoplist_ref(GtMultieoplist *multieops)
{
  multieops->refcount++;
  return multieops;
}

void gt_multieoplist_delete(GtMultieoplist *multieops)
{
  if (multieops != NULL) {
    if (multieops->refcount != 0) {
      multieops->refcount--;
    }
    else {
      gt_array_delete(multieops->meoplist);
      gt_free(multieops);
    }
  }
}

static void gt_multieoplist_add_eop(GtMultieoplist *multieops,
                                    AlignmentEoptype type)
{
  GtMultieop meop, *meop_ptr;
  if (gt_array_size(multieops->meoplist) != 0) {
    meop_ptr = gt_array_get_last(multieops->meoplist);
    if (meop_ptr->type == type) {
      meop_ptr->steps++; /* XXX: check for overflow */
      return;
    }
  }
  meop.type = type;
  meop.steps = 1UL;
  gt_array_add(multieops->meoplist, meop);
}

void gt_multieoplist_add_replacement(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eop(multieops, Replacement);
}

void gt_multieoplist_add_insertion(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eop(multieops, Insertion);
}

void gt_multieoplist_add_deletion(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eop(multieops, Deletion);
}

void gt_multieoplist_add_mismatch(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eop(multieops, Mismatch);
}

void gt_multieoplist_add_match(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eop(multieops, Match);
}

void gt_multieoplist_reset(GtMultieoplist *multieops)
{
  gt_array_reset(multieops->meoplist);
}

void gt_multieoplist_remove_last(GtMultieoplist *multieops)
{
  GtMultieop *meop_ptr;
  gt_assert(multieops && gt_array_size(multieops->meoplist));
  meop_ptr = gt_array_get_last(multieops->meoplist);
  gt_assert(meop_ptr->steps);
  if (meop_ptr->steps == 1UL)
    (void) gt_array_pop(multieops->meoplist);
  else
    meop_ptr->steps--;
}

GtUword gt_multieoplist_get_repdel_length(GtMultieoplist *multieops)
{
  GtUword len = 0, i;
  GtMultieop meop;
  for (i = gt_array_size(multieops->meoplist); i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(multieops->meoplist, i-1);
    switch (meop.type) {
      case Match:
      case Mismatch:
      case Replacement:
      case Deletion:
        len += meop.steps;
        break;
      default:
        /* nothing */;
    }
  }
  return len;
}

GtUword gt_multieoplist_get_repins_length(GtMultieoplist *multieops)
{
  GtUword len = 0, i;
  GtMultieop meop;
  for (i = gt_array_size(multieops->meoplist); i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(multieops->meoplist, i-1);
    switch (meop.type) {
      case Match:
      case Mismatch:
      case Replacement:
      case Insertion:
        len += meop.steps;
        break;
      default:
        /* nothing */;
    }
  }
  return len;
}

GtUword gt_multieoplist_get_length(GtMultieoplist *multieops)
{
  return(gt_array_size(multieops->meoplist));
}

GtMultieop *gt_multieoplist_get_entry(GtMultieoplist *multieops,
                                      GtUword index)
{
  gt_assert(index < gt_array_size(multieops->meoplist));
  return((GtMultieop *) gt_array_get(multieops->meoplist, index));
}

void gt_multieoplist_show(GtMultieoplist *multieops, FILE *fp)
{
  GtUword i, size;
  GtMultieop meop;

  gt_assert(multieops != NULL);

  size = gt_array_size(multieops->meoplist);

  gt_xfputc('[', fp);
  for (i = size; i > 0; i--) {
    meop = *(GtMultieop*) gt_array_get(multieops->meoplist, i-1);
    switch (meop.type) {
      case Match:
        gt_xfputc('M', fp);
        break;
      case Mismatch:
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
    fprintf(fp, " "GT_LU"", meop.steps);
    if (i != 1UL)
      gt_xfputc(',', fp);
  }
  gt_xfputs("]\n", fp);
}
