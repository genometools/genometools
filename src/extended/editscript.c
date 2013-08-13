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

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/editscript.h"
#include "core/ensure.h"

struct GtEditscript
{
  GtMultieoplist *multieops;
  /* XXX instead of GtArray use compactulongstore */
  GtArrayGtUchar mis_ins_chars;
};

GtEditscript *gt_editscript_new(void)
{
  GtEditscript *editscript;
  editscript = gt_calloc((size_t) 1, sizeof (GtEditscript));
  GT_INITARRAY(&editscript->mis_ins_chars, GtUchar);
  editscript->multieops = gt_multieoplist_new();
  return editscript;
}

void gt_editscript_delete(GtEditscript *editscript)
{
  if (editscript != NULL) {
    gt_multieoplist_delete(editscript->multieops);
    GT_FREEARRAY(&editscript->mis_ins_chars, GtUchar);
    gt_free(editscript);
  }
}

void gt_editscript_add_mismatch(GtEditscript *editscript, GtUchar c)
{
  gt_assert(editscript && editscript->multieops);
  gt_multieoplist_add_mismatch(editscript->multieops);
  GT_STOREINARRAY(&editscript->mis_ins_chars, GtUchar, 128, c);
}

void gt_editscript_add_match(GtEditscript *editscript)
{
  gt_assert(editscript && editscript->multieops);
  gt_multieoplist_add_match(editscript->multieops);
}

void gt_editscript_add_deletion(GtEditscript *editscript)
{
  gt_assert(editscript && editscript->multieops);
  gt_multieoplist_add_deletion(editscript->multieops);
}

void gt_editscript_add_insertion(GtEditscript *editscript, GtUchar c)
{
  gt_assert(editscript && editscript->multieops);
  gt_multieoplist_add_insertion(editscript->multieops);
  GT_STOREINARRAY(&editscript->mis_ins_chars, GtUchar, 128, c);
}

void gt_editscript_reset(GtEditscript *editscript)
{
  editscript->mis_ins_chars.nextfreeGtUchar = 0;
  gt_multieoplist_reset(editscript->multieops);
}

void gt_editscript_remove_last(GtEditscript *editscript)
{
  GtUword lastidx;
  GtMultieop meop;
  gt_assert(editscript && editscript->multieops &&
            gt_multieoplist_get_length(editscript->multieops) > 0);
  lastidx = gt_multieoplist_get_length(editscript->multieops);
  meop = gt_multieoplist_get_entry(editscript->multieops, lastidx);
  if (meop.type == Insertion || meop.type == Mismatch) {
    editscript->mis_ins_chars.nextfreeGtUchar--;
  }
  gt_multieoplist_remove_last(editscript->multieops);
}

/* XXX add check for right number of chars for mismatch and insertion. */
int gt_editscript_is_valid(const GtEditscript *editscript,
                           GtUword ulen, GtUword vlen)
{
  int had_err = 0;
  GtUword len;
  /* check ulen */
  len = gt_multieoplist_get_repdel_length(editscript->multieops);
  if (len != ulen) {
    printf("ulen: " GT_WU ", repdel: " GT_WU "\n", ulen, len);
    had_err = 1;
  }
  /* check vlen */
  len = gt_multieoplist_get_repins_length(editscript->multieops);
  if (len != vlen) {
    printf("vlen: " GT_WU ", repins: " GT_WU "\n", vlen, len);
    had_err = 1;
  }
  return had_err;
}

GtEditscript *gt_editscript_new_with_sequences(const GtEncseq *encseq,
                                               GtMultieoplist *multieops,
                                               GtUword start,
                                               GtReadmode dir)
{
  GtUword vidx, idx, meopidx, meoplen;
  GtEditscript *editscript = gt_calloc((size_t) 1, sizeof (GtEditscript));
  GtMultieop meop;
  GtUchar cchar;

  GT_INITARRAY(&editscript->mis_ins_chars, GtUchar);
  /* len of v */
  vidx = gt_multieoplist_get_repins_length(multieops) - 1;
  meoplen = gt_multieoplist_get_length(multieops);

  for (idx = 0; idx < meoplen; idx++) {
    meop = gt_multieoplist_get_entry(multieops, idx - 1);
    switch (meop.type) {
      case Match:
        for (meopidx = 0; meopidx < meop.steps; meopidx++) {
          vidx--;
        }
        break;
      case Replacement:
      case Mismatch:
        for (meopidx = 0; meopidx < meop.steps; meopidx++) {
          cchar = gt_encseq_get_encoded_char(encseq, start + vidx, dir);
          GT_STOREINARRAY(&editscript->mis_ins_chars, GtUchar, 128, cchar);
          vidx--;
        }
        break;
      case Insertion:
        for (meopidx = 0; meopidx < meop.steps; meopidx++) {
          cchar = gt_encseq_get_encoded_char(encseq, start + vidx, dir);
          GT_STOREINARRAY(&editscript->mis_ins_chars, GtUchar, 128, cchar);
          vidx--;
        }
        break;
      case Deletion:
      default:
        break;
    }
  }
  editscript->multieops = multieops;
  return(editscript);
}

GtUword gt_editscript_get_sequence(const GtEditscript *editscript,
                                   GtEncseqReader *esr,
                                   GtUchar *buffer)
{
  GtUword i, j, charidx, meoplen,
      count = 0;
  GtMultieop meop;
  GtUchar *space;

  gt_assert(esr != NULL && editscript != NULL);

  meoplen = gt_multieoplist_get_length(editscript->multieops);
  space = editscript->mis_ins_chars.spaceGtUchar;
  charidx = editscript->mis_ins_chars.nextfreeGtUchar - 1;

  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(editscript->multieops, i - 1);
    if (meop.type == Deletion) {
      for (j = 0; j < meop.steps; j++) {
        (void) gt_encseq_reader_next_encoded_char(esr);
      }
    }
    else if (meop.type == Insertion) {
      for (j = 0; j < meop.steps; j++) {
        gt_assert(charidx < editscript->mis_ins_chars.nextfreeGtUchar);
        buffer[count] = space[charidx];
        count++;
        charidx--;
      }
    }
    else if (meop.type == Match) {
      for (j = 0; j < meop.steps; j++) {
        buffer[count] = gt_encseq_reader_next_encoded_char(esr);
        count++;
      }
    }
    else if (meop.type == Mismatch) {
      for (j = 0; j < meop.steps; j++) {
        gt_assert(charidx < editscript->mis_ins_chars.nextfreeGtUchar);
        (void) gt_encseq_reader_next_encoded_char(esr);
        buffer[count] = space[charidx];
        count++;
        charidx--;
      }
    }
  }
  return(count);
}

GtMultieoplist *gt_editscript_get_multieops(GtEditscript *editscript)
{
  gt_assert(editscript);
  return(editscript->multieops);
}

GtEditscript *gt_editscript_io(GtEditscript *es, FILE *fp,
                               EditscriptIOFunc io_func)
{
  if (es == NULL) {
    es = gt_calloc((size_t) 1, sizeof (GtEditscript));
  }
  io_func(&es->mis_ins_chars.nextfreeGtUchar,
          sizeof (es->mis_ins_chars.nextfreeGtUchar), (size_t) 1, fp);
  if (es->mis_ins_chars.spaceGtUchar == NULL) {
    es->mis_ins_chars.spaceGtUchar =
      gt_malloc((size_t) es->mis_ins_chars.nextfreeGtUchar *
                sizeof (*(es->mis_ins_chars.spaceGtUchar)));
    es->mis_ins_chars.allocatedGtUchar = es->mis_ins_chars.nextfreeGtUchar;
  }
  io_func(es->mis_ins_chars.spaceGtUchar,
          sizeof (*(es->mis_ins_chars.spaceGtUchar)),
          (size_t) es->mis_ins_chars.nextfreeGtUchar,
          fp);
  es->multieops = gt_meoplist_io(es->multieops, io_func, fp);
  return(es);
}

int gt_editscript_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtUword length, i;
  GtMultieoplist *meops = gt_multieoplist_new();
  GtEncseq *v, *u;
  GtAlphabet *dna = gt_alphabet_new_dna();
  GtEncseqBuilder *esb = gt_encseq_builder_new(dna);
  GtUchar buffer[12];
  /*
    u: aaacccg-ggttt
       || | || ||| |
    v: aatc-cggggtat
    meops(reverse): MmMMMIMMDMmMM
   */
  const char *seq1 = "aaacccgggttt",
             *seq2 = "aatccggggtat";
  GtEditscript *es = gt_editscript_new();
  GtEncseqReader *esr;

  gt_error_check(err);

  gt_encseq_builder_add_cstr(esb, seq1, 12UL, "u");
  u = gt_encseq_builder_build(esb, err);
  gt_encseq_builder_reset(esb);
  gt_encseq_builder_add_cstr(esb, seq2, 12UL, "v");
  v = gt_encseq_builder_build(esb, err);
  gt_encseq_builder_delete(esb);

  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);
  gt_multieoplist_add_mismatch(meops);
  gt_editscript_add_mismatch(es,
                             gt_encseq_get_encoded_char(v, 10UL,
                                                        GT_READMODE_FORWARD));
  for (i = 0; i < 3UL; ++i) {
    gt_multieoplist_add_match(meops);
    gt_editscript_add_match(es);
  }
  gt_multieoplist_add_insertion(meops);
  gt_editscript_add_insertion(es,
                              gt_encseq_get_encoded_char(v, 6UL,
                                                         GT_READMODE_FORWARD));
  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);
  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);
  gt_multieoplist_add_deletion(meops);
  gt_editscript_add_deletion(es);
  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);
  gt_multieoplist_add_mismatch(meops);
  gt_editscript_add_mismatch(es,
                             gt_encseq_get_encoded_char(v, 2UL,
                                                        GT_READMODE_FORWARD));
  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);
  gt_multieoplist_add_match(meops);
  gt_editscript_add_match(es);

  for (i = 0; i < es->mis_ins_chars.nextfreeGtUchar; ++i) {
    printf("stored: %c\n",
           gt_alphabet_decode(dna, es->mis_ins_chars.spaceGtUchar[i]));
  }
  esr = gt_encseq_create_reader_with_readmode(u, GT_READMODE_FORWARD, 0);
  length = gt_editscript_get_sequence(es, esr, buffer);
  gt_ensure(length == 12UL);
  for (i = 0; !had_err && i < 12UL; ++i) {
    printf("in buf: %c\n", gt_alphabet_decode(dna, buffer[i]));
    printf("in encseq: %c\n", gt_encseq_get_decoded_char(v, (GtUword) i,
                                              GT_READMODE_FORWARD));
    gt_ensure(buffer[i] == gt_encseq_get_encoded_char(v, (GtUword) i,
                                                      GT_READMODE_FORWARD));
  }
  gt_editscript_delete(es);
  if (!had_err) {
    gt_encseq_reader_reinit_with_readmode(esr, u, GT_READMODE_FORWARD, 0);
    es = gt_editscript_new_with_sequences(v, meops, 0, GT_READMODE_FORWARD);
    for (i = 0; i < es->mis_ins_chars.nextfreeGtUchar; ++i) {
      printf("stored: %c\n",
             gt_alphabet_decode(dna, es->mis_ins_chars.spaceGtUchar[i]));
    }
    length = gt_editscript_get_sequence(es, esr, buffer);
    gt_ensure(length == 12UL);
    for (i = 0; !had_err && i < 12UL; ++i) {
      printf("in buf: %c\n", gt_alphabet_decode(dna, buffer[i]));
      printf("in encseq: %c\n", gt_encseq_get_decoded_char(v, (GtUword) i,
                                                GT_READMODE_FORWARD));
      gt_ensure(buffer[i] == gt_encseq_get_encoded_char(v, (GtUword) i,
                                                        GT_READMODE_FORWARD));
    }
    gt_editscript_delete(es);
  }
  gt_encseq_reader_delete(esr);
  gt_alphabet_delete(dna);
  gt_encseq_delete(u);
  gt_encseq_delete(v);
  return had_err;
}
