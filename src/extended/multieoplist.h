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

#ifndef MULTIEOPLIST_H
#define MULTIEOPLIST_H

#include "core/array.h"
#include "core/types_api.h"
#include "core/encseq_api.h"

typedef struct GtMultieoplist GtMultieoplist;

typedef enum {
  Replacement,
  Deletion,
  Insertion,
  Match,
  Mismatch
} AlignmentEoptype;

typedef struct {
  AlignmentEoptype type;
  unsigned long steps;
} GtMultieop;
/* XXX: possible improvement to save memory: combine both parts into a single
        variable (use bit shifting operations) */

GtMultieoplist *gt_multieoplist_new(void);
void gt_multieoplist_delete(GtMultieoplist *multieops);

void gt_multieoplist_add_replacement(GtMultieoplist *eops);
void gt_multieoplist_add_insertion(GtMultieoplist *eops);
void gt_multieoplist_add_deletion(GtMultieoplist *eops);
void gt_multieoplist_add_mismatch(GtMultieoplist *eops);
void gt_multieoplist_add_match(GtMultieoplist *eops);
void gt_multieoplist_reset(GtMultieoplist *eops);
void gt_multieoplist_remove_last(GtMultieoplist *eops);
unsigned long gt_multieoplist_get_length(GtMultieoplist *eops);
GtMultieop *gt_multieoplist_get_entry(GtMultieoplist *eops,
    unsigned long index);

unsigned long gt_multieoplist_get_repdel_length(GtMultieoplist *eops);
unsigned long gt_multieoplist_get_repins_length(GtMultieoplist *eops);

void gt_multieoplist_show(GtMultieoplist *eops, FILE *fp);

#endif
