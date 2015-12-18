/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/xansi_api.h"
#include "extended/multieoplist.h"
#include "core/log_api.h"

typedef uint8_t Eop;

GT_DECLAREARRAYSTRUCT(Eop);

struct GtMultieoplist {
  GtUword refcount;
  GtArrayEop meoplist;
};

GtMultieoplist *gt_multieoplist_new(void)
{
  GtMultieoplist *multieops;
  multieops = gt_malloc(sizeof (GtMultieoplist));
  multieops->refcount = 0;
  GT_INITARRAY(&multieops->meoplist, Eop);
  return multieops;
}

GtMultieoplist *gt_multieoplist_new_with_size(GtUword size)
{
  GtMultieoplist *multieops;
  multieops = gt_malloc(sizeof (GtMultieoplist));
  multieops->refcount = 0;
  GT_INITARRAY(&multieops->meoplist, Eop);
  multieops->meoplist.spaceEop =
    gt_calloc((size_t) size, sizeof (*(multieops->meoplist.spaceEop)));
  multieops->meoplist.allocatedEop = size;
  return multieops;
}

GtMultieoplist *gt_multieoplist_clone(GtMultieoplist *copy,
                                      GtMultieoplist *source)
{
  GtUword i;
  gt_assert(source != NULL);
  if (copy == NULL) {
    copy = gt_multieoplist_new();
  }
  if (copy->meoplist.allocatedEop < source->meoplist.nextfreeEop) {
    copy->meoplist.spaceEop = gt_realloc(copy->meoplist.spaceEop,
                                         (size_t) source->meoplist.nextfreeEop);
    copy->meoplist.allocatedEop = source->meoplist.nextfreeEop;
  }
  copy->refcount = 0;
  copy->meoplist.nextfreeEop = source->meoplist.nextfreeEop;
  for (i = 0; i < copy->meoplist.nextfreeEop; i++) {
    copy->meoplist.spaceEop[i] = source->meoplist.spaceEop[i];
  }
  return copy;
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
      GT_FREEARRAY(&multieops->meoplist, Eop);
      gt_free(multieops);
    }
  }
}

#define GT_MEOPS_MATCH ((Eop) 0)
#define GT_MEOPS_MIS ((Eop) 1)
#define GT_MEOPS_DEL ((Eop) 2)
#define GT_MEOPS_INS ((Eop) 3)
#define GT_MEOPS_STEPS_BITS ((sizeof (Eop) * CHAR_BIT) - 2)
#define GT_MEOPS_STEPS_MASK ((((Eop) 1) << GT_MEOPS_STEPS_BITS) - 1)

static void gt_multieoplist_add_eops(GtMultieoplist *multieops,
                                     AlignmentEoptype type,
                                     GtUword steps)
{
  Eop *space, type_bits = 0;

  gt_assert(multieops != NULL);
  space = multieops->meoplist.spaceEop;
  switch (type) {
    case Deletion:
      type_bits = GT_MEOPS_DEL;
      break;
    case Insertion:
      type_bits = GT_MEOPS_INS;
      break;
    case Match:
    case Replacement:
      type_bits = GT_MEOPS_MATCH;
      break;
    case Mismatch:
      type_bits = GT_MEOPS_MIS;
      break;
    default:
      gt_assert(false);
  }
  if (multieops->meoplist.nextfreeEop != 0) {
    GtUword current = multieops->meoplist.nextfreeEop - 1;
    if (space[current] >> GT_MEOPS_STEPS_BITS == type_bits) {
      while (steps != 0 &&
             (space[current] & GT_MEOPS_STEPS_MASK) < GT_MEOPS_STEPS_MASK) {
        space[current]++;
        steps--;
      }
    }
  }
  if (steps != 0) {
    while (steps != 0) {
      Eop tmp = type_bits << GT_MEOPS_STEPS_BITS;
      if (steps < (GtUword) GT_MEOPS_STEPS_MASK) {
        tmp += steps;
        steps = 0;
      }
      else {
        tmp |= GT_MEOPS_STEPS_MASK;
        steps -= GT_MEOPS_STEPS_MASK;
      }
      GT_STOREINARRAY(&multieops->meoplist, Eop, 64<<2, tmp);
    }
  }
}

void gt_multieoplist_add_replacement(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eops(multieops, Replacement, (GtUword) 1);
}

void gt_multieoplist_add_replacement_multi(GtMultieoplist *multieops,
                                           GtUword num)
{
  gt_assert(num > 0);
  gt_multieoplist_add_eops(multieops, Replacement, num);
}

void gt_multieoplist_add_insertion(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eops(multieops, Insertion, (GtUword) 1);
}

void gt_multieoplist_add_deletion(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eops(multieops, Deletion, (GtUword) 1);
}

void gt_multieoplist_add_mismatch(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eops(multieops, Mismatch, (GtUword) 1);
}

void gt_multieoplist_add_match(GtMultieoplist *multieops)
{
  gt_multieoplist_add_eops(multieops, Match, (GtUword) 1);
}

void gt_multieoplist_reset(GtMultieoplist *multieops)
{
  multieops->meoplist.nextfreeEop = 0;
}

void gt_multieoplist_remove_last(GtMultieoplist *multieops)
{
  Eop *space;
  GtUword current;
  gt_assert(multieops != NULL && multieops->meoplist.nextfreeEop > 0);
  space = multieops->meoplist.spaceEop;
  current = multieops->meoplist.nextfreeEop - 1;
  if ((space[current] & GT_MEOPS_STEPS_MASK) == (Eop) 1) {
    multieops->meoplist.nextfreeEop--;
  }
  else {
    space[current]--;
  }
}

GtUword gt_multieoplist_get_repdel_length(GtMultieoplist *multieops)
{
  GtUword len = 0, i;
  Eop *space;
  gt_assert(multieops);
  space = multieops->meoplist.spaceEop;
  for (i = 0; i < multieops->meoplist.nextfreeEop ; i++) {
    switch (space[i] >> GT_MEOPS_STEPS_BITS) {
      case GT_MEOPS_INS:
        break;
      default:
        len += space[i] & GT_MEOPS_STEPS_MASK;
    }
  }
  return len;
}

GtUword gt_multieoplist_get_repins_length(GtMultieoplist *multieops)
{
  GtUword len = 0, i;
  Eop *space;
  gt_assert(multieops);
  space = multieops->meoplist.spaceEop;
  for (i = 0; i < multieops->meoplist.nextfreeEop ; i++) {
    switch (space[i] >> GT_MEOPS_STEPS_BITS) {
      case GT_MEOPS_DEL:
        break;
      default:
        len += space[i] & GT_MEOPS_STEPS_MASK;
    }
  }
  return len;
}

GtUword gt_multieoplist_get_length(GtMultieoplist *multieops)
{
  GtUword len = 0, i;
  Eop *space;
  gt_assert(multieops);
  space = multieops->meoplist.spaceEop;
  for (i = 0; i < multieops->meoplist.nextfreeEop ; i++) {
    len += space[i] & GT_MEOPS_STEPS_MASK;
  }
  return len;
}

GtUword gt_multieoplist_get_num_entries(GtMultieoplist *multieops)
{
  return(multieops->meoplist.nextfreeEop);
}

GtMultieop gt_multieoplist_get_entry(const GtMultieoplist *multieops,
                                     GtUword index)
{
  GtMultieop eop = { (AlignmentEoptype) 0, (GtUword) 0};
  Eop *space;
  gt_assert(multieops);
  space = multieops->meoplist.spaceEop;
  gt_assert(multieops->meoplist.nextfreeEop != 0);
  gt_assert(multieops->meoplist.nextfreeEop > index);
  eop.steps = (GtUword) space[index] & (GtUword) GT_MEOPS_STEPS_MASK;
  switch (space[index] >> GT_MEOPS_STEPS_BITS) {
    case GT_MEOPS_MATCH:
      eop.type = Match;
      break;
    case GT_MEOPS_MIS:
      eop.type = Mismatch;
      break;
    case GT_MEOPS_DEL:
      eop.type = Deletion;
      break;
    case GT_MEOPS_INS:
      eop.type = Insertion;
  }
  return eop;
}

void gt_multieoplist_show(GtMultieoplist *multieops, FILE *fp)
{
  GtUword stepssum,
          num = multieops->meoplist.nextfreeEop;
  Eop *space, *last, *next;

  if (num == 0)
  {
    gt_xfputs("[]\n", fp);
    return;
  }
  gt_assert(multieops != NULL);
  space = multieops->meoplist.spaceEop;

  gt_xfputc('[', fp);
  last = space + num - 1;
  next = space + num - 2;
  stepssum = (GtUword) (*last & GT_MEOPS_STEPS_MASK);
  while (next >= space) {
    if (*last >> GT_MEOPS_STEPS_BITS ==
        *next >> GT_MEOPS_STEPS_BITS) {
      stepssum += *next & GT_MEOPS_STEPS_MASK;
      last--;
      next--;
    }
    else {
      switch (*last >> GT_MEOPS_STEPS_BITS) {
        case GT_MEOPS_MATCH:
          gt_xfputc('M', fp);
          break;
        case GT_MEOPS_MIS:
          gt_xfputc('R', fp);
          break;
        case GT_MEOPS_INS:
          gt_xfputc('I', fp);
          break;
        case GT_MEOPS_DEL:
          gt_xfputc('D', fp);
          break;
      }
      fprintf(fp, " %u, ", (unsigned int) stepssum);
      last--;
      stepssum = (GtUword) (*last & GT_MEOPS_STEPS_MASK);
      next = last - 1;
    }
    gt_assert(next + 1 == last);
  }
  gt_assert(num == 0 || last == space);
  switch (*last >> GT_MEOPS_STEPS_BITS) {
    case GT_MEOPS_MATCH:
      gt_xfputc('M', fp);
      break;
    case GT_MEOPS_MIS:
      gt_xfputc('R', fp);
      break;
    case GT_MEOPS_INS:
      gt_xfputc('I', fp);
      break;
    case GT_MEOPS_DEL:
      gt_xfputc('D', fp);
      break;
  }
  fprintf(fp, " %u", (unsigned int) stepssum);
  gt_xfputs("]\n", fp);
}

void gt_multieoplist_combine(GtMultieoplist *multieops,
                             const GtMultieoplist *multieops_to_add,
                             bool forward) {
  GtUword idx,
          nelems = 0;
  GtMultieop mop;
  gt_assert(multieops_to_add != NULL && multieops != NULL);

  nelems = multieops_to_add->meoplist.nextfreeEop;
  for (idx = 0; idx < nelems; ++idx) {
    mop = gt_multieoplist_get_entry(multieops_to_add,
                                    forward ?
                                    idx :
                                    nelems - idx - 1);
    gt_multieoplist_add_eops(multieops, mop.type, mop.steps);
  }
}

static GtMultieoplist *gt_multieoplist_io_fp(GtMultieoplist *multieops,
                                             FILE *fp, GtIOFunc io_func,
                                             GtError *err)
{
  int had_err = 0;
  had_err = io_func(&multieops->meoplist.nextfreeEop,
                    sizeof (multieops->meoplist.nextfreeEop),
                    (size_t) 1, fp, err);
  if (!had_err) {
    gt_assert(multieops->meoplist.nextfreeEop != 0);
    if (multieops->meoplist.spaceEop == NULL) {
      multieops->meoplist.allocatedEop = multieops->meoplist.nextfreeEop;
      multieops->meoplist.spaceEop =
        gt_malloc((size_t) multieops->meoplist.nextfreeEop *
                  sizeof (*(multieops->meoplist.spaceEop)));
    }
  }
  had_err = io_func(multieops->meoplist.spaceEop,
                    sizeof (*(multieops->meoplist.spaceEop)),
                    (size_t) multieops->meoplist.nextfreeEop,
                    fp, err);
  if (had_err) {
    gt_multieoplist_delete(multieops);
    multieops = NULL;
  }
  return(multieops);
}

GtMultieoplist *gt_multieoplist_io(GtMultieoplist *multieops, FILE *fp,
                               GtError *err)
{
  if (multieops == NULL) {
    multieops = gt_calloc((size_t) 1, sizeof (GtMultieoplist));
    GT_INITARRAY(&multieops->meoplist, Eop);
    multieops = gt_multieoplist_io_fp(multieops, fp, gt_io_error_fread, err);
  }
  else
    multieops = gt_multieoplist_io_fp(multieops, fp, gt_io_error_fwrite, err);
  return multieops;
}

int gt_multieoplist_unit_test_2(GtError *err)
{
  /* char *u = "attttgatatcgtctctctatgctgtcac", */
       /* *v = "atttngatatcgtctctctatgctgtcac"; */
  GtMultieoplist *meop = gt_multieoplist_new();
  unsigned int i;
  int had_err = 0;

  for (i = 0; i < 4; ++i) {
    gt_multieoplist_add_match(meop);
  }
  gt_multieoplist_add_mismatch(meop);
  for (i = 0; i < 24; ++i) {
    gt_multieoplist_add_match(meop);
  }
  gt_ensure(gt_multieoplist_get_length(meop) == 29);
  gt_log_log("length: " GT_WU, gt_multieoplist_get_length(meop));

  gt_multieoplist_delete(meop);
  return had_err;
}

int gt_multieoplist_unit_test(GtError *err)
{
  int had_err = 0;
  GtUword idx, num;
  GtMultieoplist *list = gt_multieoplist_new(),
                 *list2 = NULL;

  gt_error_check(err);

  gt_multieoplist_add_deletion(list);
  gt_ensure(gt_multieoplist_get_num_entries(list) == (GtUword) 1);
  gt_ensure(list->meoplist.spaceEop[0] >> GT_MEOPS_STEPS_BITS == GT_MEOPS_DEL);
  gt_ensure((int) (list->meoplist.spaceEop[0] & GT_MEOPS_STEPS_MASK) == 1);
  if (!had_err) {
    GtMultieop meop = gt_multieoplist_get_entry(list, 0);
    gt_ensure(meop.type = Deletion);
    gt_ensure(meop.steps == (GtUword) 1UL);
  }

  for (idx = 0; idx < (GtUword) GT_MEOPS_STEPS_MASK; idx++) {
    gt_multieoplist_add_deletion(list);
  }
  gt_ensure(gt_multieoplist_get_num_entries(list) == (GtUword) 2);
  gt_ensure(list->meoplist.spaceEop[1] >> GT_MEOPS_STEPS_BITS == GT_MEOPS_DEL);
  gt_ensure((GtUword) (list->meoplist.spaceEop[1] & GT_MEOPS_STEPS_MASK) ==
            (GtUword) 1);
  gt_ensure((GtUword) (list->meoplist.spaceEop[0] & GT_MEOPS_STEPS_MASK) ==
            (GtUword) GT_MEOPS_STEPS_MASK);
  gt_ensure(gt_multieoplist_get_repdel_length(list) ==
            (GtUword) GT_MEOPS_STEPS_MASK + 1);

  if (!had_err)
    gt_multieoplist_add_eops(list, Match, (GtUword) (GT_MEOPS_STEPS_MASK >> 1));
  gt_ensure(gt_multieoplist_get_num_entries(list) == (GtUword) 3);
  gt_ensure(list->meoplist.spaceEop[2] >> GT_MEOPS_STEPS_BITS ==
            GT_MEOPS_MATCH);
  gt_ensure((GtUword) (list->meoplist.spaceEop[2] & GT_MEOPS_STEPS_MASK) ==
            (GtUword) GT_MEOPS_STEPS_MASK >> 1);

  if (!had_err)
    gt_multieoplist_add_eops(list, Match, (GtUword) (GT_MEOPS_STEPS_MASK << 1));
  gt_ensure(gt_multieoplist_get_num_entries(list) == (GtUword) 5);
  gt_ensure(list->meoplist.spaceEop[4] >> GT_MEOPS_STEPS_BITS ==
            GT_MEOPS_MATCH);
  gt_ensure((list->meoplist.spaceEop[4] & GT_MEOPS_STEPS_MASK) ==
            GT_MEOPS_STEPS_MASK >> 1);

  if (!had_err) {
    list2 = gt_multieoplist_clone(list2, list);
  }

  gt_ensure(gt_multieoplist_get_length(list) ==
            gt_multieoplist_get_length(list2));
  gt_ensure(gt_multieoplist_get_repdel_length(list) ==
            gt_multieoplist_get_repdel_length(list2));
  if (!had_err) {
    num = gt_multieoplist_get_num_entries(list);

    for (idx = 0; idx < num; idx++) {
      GtMultieop eop1 = gt_multieoplist_get_entry(list, idx),
                 eop2 = gt_multieoplist_get_entry(list2, idx);
      gt_ensure(eop1.type == eop2.type);
      gt_ensure(eop1.steps == eop2.steps);
    }
  }

  if (!had_err) {
    gt_multieoplist_combine(list, list2, true);
  }
  gt_multieoplist_delete(list);
  gt_multieoplist_delete(list2);

  had_err = gt_multieoplist_unit_test_2(err);
  return had_err;
}
