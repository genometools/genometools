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

#include "extended/editscript.h"
#include "core/ma.h"
#include "core/xansi_api.h"

struct GtEditscript
{
  GtUword encseq_orig_start;
  GtUword encseq_orig_end;
  GtUword encseq_align_start;
  GtUword encseq_align_end;
  GtUword encseq_compressed_start;
  GtUword seqlen;
  GtMultieoplist *eops;
  char *charlist;
  GtUword charidx;
  GtUword maxchar;
};

GtEditscript *gt_editscript_new(void)
{
  GtEditscript *editscript;
  editscript = gt_malloc(sizeof (GtEditscript));
  editscript->charlist = gt_malloc(sizeof (GtUchar) * 100);
  editscript->eops = gt_multieoplist_new();
  editscript->seqlen = 0;
  editscript->charidx = 0;
  editscript->maxchar = 100UL;
  return editscript;
}

void gt_editscript_delete(GtEditscript *editscript)
{
  gt_multieoplist_delete(editscript->eops);
  gt_free(editscript->charlist);
  gt_free(editscript);
}

static void gt_editscript_add_char_to_charlist(GtEditscript *editscript, char c)
{
  if (editscript->charidx == editscript->maxchar) {
    char *charlist = gt_realloc(editscript->charlist,
                                sizeof (char) * (editscript->maxchar + 100));
    editscript->charlist = charlist;
    editscript->maxchar += 100;
  }
  editscript->charlist[editscript->charidx] = c;
  editscript->charidx++;
}

void gt_editscript_add_mismatch(GtEditscript *editscript, char c)
{
  gt_assert(editscript && editscript->eops);
  gt_multieoplist_add_mismatch(editscript->eops);
  gt_editscript_add_char_to_charlist(editscript, c);
  editscript->seqlen++;
}

void gt_editscript_add_match(GtEditscript *editscript)
{
  gt_assert(editscript && editscript->eops);
  gt_multieoplist_add_match(editscript->eops);
  editscript->seqlen++;
}

void gt_editscript_add_deletion(GtEditscript *editscript)
{
  gt_assert(editscript && editscript->eops);
  gt_multieoplist_add_deletion(editscript->eops);
}

void gt_editscript_add_insertion(GtEditscript *editscript, char c)
{
  gt_assert(editscript && editscript->eops);
  gt_multieoplist_add_insertion(editscript->eops);
  gt_editscript_add_char_to_charlist(editscript, c);
  editscript->seqlen++;
}

void gt_editscript_reset(GtEditscript *editscript)
{
  gt_multieoplist_reset(editscript->eops);
  editscript->seqlen = 0;
}

void gt_editscript_remove_last(GtEditscript *editscript)
{
  gt_assert(editscript && editscript->eops);
  gt_multieoplist_remove_last(editscript->eops);
  editscript->seqlen--;
}

void gt_editscript_show(const GtEncseq *encseq, const GtEditscript *editscript,
                        FILE* fp)
{
  GtUword i, j, charidx = editscript->seqlen - 1, meoplen;
  GtEncseqReader *esr;
  GtMultieop *meop;
  char c;

  gt_assert(encseq != NULL && editscript != NULL);

  meoplen = gt_multieoplist_get_length(editscript->eops);
  esr =
    gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_FORWARD,
                                          editscript->encseq_compressed_start);

  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(editscript->eops, i - 1);
    if (meop->type == Deletion) {
      for (j = 0; j < meop->steps; j++) {
        (void) gt_encseq_reader_next_decoded_char(esr);
      }
    }
    else if (meop->type == Insertion) {
      for (j = 0; j < meop->steps; j++) {
        gt_xfputc(editscript->charlist[charidx], fp);
        charidx--;
      }
    }
    else if (meop->type == Match) {
      for (j = 0; j < meop->steps; j++) {
        c = gt_encseq_reader_next_decoded_char(esr);
        gt_xfputc(c, fp);
      }
    }
    else if (meop->type == Mismatch) {
      for (j = 0; j < meop->steps; j++) {
        (void) gt_encseq_reader_next_decoded_char(esr);
        gt_xfputc(editscript->charlist[charidx], fp);
        charidx--;
      }
    }
  }
  gt_encseq_reader_delete(esr);
}

GtUword gt_editscript_get_orig_start(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->encseq_orig_start);
}

void gt_editscript_set_orig_start(GtEditscript *editscript, GtUword value)
{
  gt_assert(editscript);
  editscript->encseq_orig_start = value;
}

GtUword gt_editscript_get_orig_end(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->encseq_orig_end);
}

void gt_editscript_set_orig_end(GtEditscript *editscript, GtUword value)
{
  gt_assert(editscript);
  editscript->encseq_orig_end = value;
}

GtUword gt_editscript_get_align_start(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->encseq_align_start);
}

void gt_editscript_set_align_start(GtEditscript *editscript,
    GtUword value)
{
  gt_assert(editscript);
  editscript->encseq_align_start = value;
}

GtUword gt_editscript_get_align_end(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->encseq_align_end);
}

void gt_editscript_set_align_end(GtEditscript *editscript, GtUword value)
{
  gt_assert(editscript);
  editscript->encseq_align_end = value;
}

GtUword gt_editscript_get_seqlen(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->seqlen);
}
