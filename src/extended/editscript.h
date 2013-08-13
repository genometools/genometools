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

#ifndef EDITSCRIPT_H
#define EDITSCRIPT_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "extended/multieoplist.h"

/* Class <GtEditscript> stores everything from an alignment of two sequences
   that is needed to reproduce sequence 2 (usualy called v) from sequence one
   (u). The purpose is reduction in space for sequence v. */
typedef struct GtEditscript GtEditscript;

/* Function Pointer to read or write to/from file. Only wrappers around fread
   und fwrite to have similar signatures. */
typedef void (*EditscriptIOFunc)(void *ptr,
                                 size_t size,
                                 size_t nmemb,
                                 FILE *stream);

/* Returns new empty <GtEditscript> object. */
GtEditscript*   gt_editscript_new(void);

/* Add one deletion to <editscript>. */
void            gt_editscript_add_deletion(GtEditscript *editscript);
/* Add one insertion to <editscript>. <c> is the encoded character in v that
   was inserted. */
void            gt_editscript_add_insertion(GtEditscript *editscript,
                                            GtUchar c);
/* Add one match to <editscript>. */
void            gt_editscript_add_match(GtEditscript *editscript);
/* Add one mismatch to <editscript>. <c> is the encoded character in v that
   replaces one character in u. */
void            gt_editscript_add_mismatch(GtEditscript *editscript, GtUchar c);
/* Delete content of <editscript> */
void            gt_editscript_reset(GtEditscript *editscript);
/* Remove one edit operation from <editscript>, this will be the last one added
   with one of the gt_editscript_add_*() functions. */
void            gt_editscript_remove_last(GtEditscript *editscript);
/* Returns new <GtEditscript> object and fills it with <encseq> as "second" or v
   sequence and <multieops> representing the edit operations when aligned to a
   sequence u. <start> and <dir> define where sequence v starts within
   <encseq> */
GtEditscript*   gt_editscript_new_with_sequences(const GtEncseq *encseq,
                                                 GtMultieoplist *multieops,
                                                 GtUword start,
                                                 GtReadmode dir);
/* Returns the multi edit operation list from <editscript> */
GtMultieoplist* gt_editscript_get_multieops(GtEditscript *editscript);
/* Returns 0 if ulen = Match + Mismatch + Deletion and vlen = Match + Mismatch
   + Insertion, 1 otherwise. */
int             gt_editscript_is_valid(const GtEditscript *editscript,
                                       GtUword ulen,
                                       GtUword vlen);

/* Transform sequence u represeted by <esr> using the edit operations in
   <editscript> and stores the result in <buffer> which has to be large enough
   to hold v. Returns the number of <GtUchar>s written to buffer. */
GtUword         gt_editscript_get_sequence(const GtEditscript *editscript,
                                           GtEncseqReader *esr,
                                           GtUchar *buffer);

/* Write or read to <fp> depending on what <io_func> is given. */
GtEditscript*   gt_editscript_io(GtEditscript *es, FILE *fp,
                                 EditscriptIOFunc io_func);

void            gt_editscript_delete(GtEditscript *editscript);

int             gt_editscript_unit_test(GtError *err);
#endif
