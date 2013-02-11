/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/fa.h"
#include "core/ma.h"
#include "match/rdj-twobitenc-editor.h"

struct GtTwobitencEditor
{
  GtTwobitencoding *twobitencoding;
  unsigned long *charcount;
  unsigned char *mapptr;
};

static int gt_twobitenc_editor_check(const GtEncseq *encseq, GtError *err)
{
  int had_err = 0;
  if (!gt_encseq_has_twobitencoding(encseq))
  {
    gt_error_set(err, "encseq correction is only "
        "possible if the sequence has a twobitencoding");
    had_err = -1;
  }
  else if (gt_encseq_accesstype_get(encseq) != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_error_set(err, "twobitencoding correction is currently only "
        "implemented if the sequence access type is EQUALLENGTH");
    had_err = -1;
  }
  return had_err;
}

GtTwobitencEditor *gt_twobitenc_editor_new(const GtEncseq *encseq,
    const char* indexname, GtError *err)
{
  GtTwobitencEditor *twobitenc_editor;
  size_t t_offset;
  size_t c_offset;
  GtStr *encseqfilename;
  int had_err = 0;

  twobitenc_editor = gt_malloc(sizeof (GtTwobitencEditor));
  had_err = gt_twobitenc_editor_check(encseq, err);
  if (had_err == 0)
  {
    t_offset = gt_encseq_twobitencoding_mapoffset(encseq);
    c_offset = gt_encseq_chardistri_mapoffset(encseq);
    encseqfilename = gt_str_new_cstr(indexname);
    gt_str_append_cstr(encseqfilename, GT_ENCSEQFILESUFFIX);
    twobitenc_editor->mapptr = (unsigned char*)gt_fa_mmap_write(gt_str_get(
          encseqfilename), NULL, err);
    twobitenc_editor->twobitencoding = (GtTwobitencoding*)
      (twobitenc_editor->mapptr + t_offset);
    twobitenc_editor->charcount = (unsigned long*)
      (twobitenc_editor->mapptr + c_offset);
    gt_str_delete(encseqfilename);
  }
  return (had_err == 0) ? twobitenc_editor : NULL;
}

void gt_twobitenc_editor_edit(GtTwobitencEditor *twobitenc_editor,
    unsigned long pos, GtUchar newchar)
{
  size_t codenum, posincode;
  GtTwobitencoding oldcode, newcode;
  GtUchar oldchar;

  gt_assert(twobitenc_editor);
  codenum = (size_t)pos / GT_UNITSIN2BITENC;
  oldcode = twobitenc_editor->twobitencoding[codenum];
  posincode = (GT_UNITSIN2BITENC - 1 -
      ((size_t)pos % GT_UNITSIN2BITENC)) << 1;
  oldchar = (oldcode & ((GtTwobitencoding)3 << posincode)) >> posincode;
  newcode = (oldcode & (~((GtTwobitencoding)3 << posincode)));
  newcode |= ((GtTwobitencoding)newchar << posincode);
  twobitenc_editor->twobitencoding[codenum] = newcode;

  /* fix counts */
  twobitenc_editor->charcount[oldchar]--;
  twobitenc_editor->charcount[newchar]++;
}

void gt_twobitenc_editor_delete(GtTwobitencEditor *twobitenc_editor)
{
  gt_assert(twobitenc_editor);
  gt_fa_xmunmap(twobitenc_editor->mapptr);
  gt_free(twobitenc_editor);
}
