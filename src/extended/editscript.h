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

#ifndef EDITSCRIPT_H
#define EDITSCRIPT_H

#include "extended/multieoplist.h"
#include "core/encseq_api.h"

typedef struct GtEditscript GtEditscript;

GtEditscript *gt_editscript_new(void);
void gt_editscript_delete(GtEditscript *editscript);

void gt_editscript_add_deletion(GtEditscript *editscript);
void gt_editscript_add_insertion(GtEditscript *editscript, char c);
void gt_editscript_add_match(GtEditscript *editscript);
void gt_editscript_add_mismatch(GtEditscript *editscript, char c);
void gt_editscript_reset(GtEditscript *editscript);
void gt_editscript_remove_last(GtEditscript *editscript);

unsigned long gt_editscript_get_orig_start(GtEditscript *editscript);
void gt_editscript_set_orig_start(GtEditscript *editscript,
    unsigned long value);
unsigned long gt_editscript_get_orig_end(GtEditscript *editscript);
void gt_editscript_set_orig_end(GtEditscript *editscript, unsigned long value);
unsigned long gt_editscript_get_align_start(GtEditscript *editscript);
void gt_editscript_set_align_start(GtEditscript *editscript,
    unsigned long value);
unsigned long gt_editscript_get_align_end(GtEditscript *editscript);
void gt_editscript_set_align_end(GtEditscript *editscript, unsigned long value);
unsigned long gt_editscript_get_seqlen(GtEditscript *editscript);

void gt_editscript_show_sequence(const GtEncseq *encseq,
    const GtEditscript *editscript, FILE *fp);

#endif
