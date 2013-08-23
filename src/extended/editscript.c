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

#include <stdbool.h>

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/chardef.h"
#include "core/ensure.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/editscript.h"

typedef struct GtEditscriptPos {
  GtUword      cur_word;
  unsigned int bitsleft;
}GtEditscriptPos;

struct GtEditscript {
  GtEditscriptPos fillpos;
  GtBitsequence  *space;
  GtBitsequence   firstmask,
                  fullmask,
                  del,
                  misdel,
                  ins,
                  last;
  size_t          size;
  GtUword         op_count,
                  space_elems,
                  trailing_matches,
                  ulen,
                  vlen;
  unsigned int    entry_size;
};

static inline void gt_editscript_pos_reset(GtEditscriptPos *pos)
{
  pos->cur_word = 0;
  pos->bitsleft = (unsigned int) GT_INTWORDSIZE;
}

GtEditscript *gt_editscript_new(GtAlphabet *alphabet)
{
  unsigned int alphabet_size;
  GtEditscript *es;
  es = gt_calloc((size_t) 1, sizeof (*es));
  alphabet_size = gt_alphabet_size(alphabet);
  es->entry_size =
    gt_determinebitspervalue((GtUword) alphabet_size + 3);
  es->size = (size_t) ((GtUword) 1 +
                     ((GtUword) 10 / (GT_INTWORDSIZE / es->entry_size)));
  gt_editscript_pos_reset(&es->fillpos);
  es->space = gt_calloc(es->size, sizeof (*(es->space)));
  es->fullmask = (GtBitsequence) ((GtUword) 1 << es->entry_size) - (GtUword) 1;
  es->firstmask = (GtBitsequence) ((GtUword) 1 <<
                                   (es->entry_size - (GtUword) 1));
  es->del = (GtBitsequence) alphabet_size;
  es->misdel = (GtBitsequence) alphabet_size + 1;
  es->ins = (GtBitsequence) alphabet_size + 2;
  es->op_count = 0;
  es->ulen = 0;
  es->space_elems = 0;
  return es;
}

void gt_editscript_delete(GtEditscript *editscript)
{
  if (editscript != NULL) {
    gt_free(editscript->space);
    gt_free(editscript);
  }
}

static inline void gt_editscript_space_add_next(GtEditscript *es,
                                                GtBitsequence value)
{
  GtUword cur_word = es->fillpos.cur_word;
  GtBitsequence remaining = (GtBitsequence) es->fillpos.bitsleft,
                bits2store = (GtBitsequence) es->entry_size;
  gt_assert(value <= es->fullmask);
  es->space_elems++;
  if (remaining == 0) {
    if (cur_word == (GtUword) es->size - 1) {
      es->size = GT_MULT2(es->size);
      es->space = gt_realloc(es->space, es->size * sizeof (*(es->space)));
    }
    cur_word++;
    es->space[cur_word] = 0;
    remaining = (GtBitsequence) GT_INTWORDSIZE;
  }
  if (remaining >= bits2store) {
    es->space[cur_word] |= value << (remaining - bits2store);
    remaining -= bits2store;
  }
  else {
    es->space[cur_word] |= value >> (bits2store - remaining);
    if (cur_word == (GtUword) es->size - 1) {
      es->size = GT_MULT2(es->size);
      es->space = gt_realloc(es->space, es->size * sizeof (*(es->space)));
    }
    cur_word++;
    es->space[cur_word] = 0;
    bits2store -= remaining;
    remaining = (GtBitsequence) GT_INTWORDSIZE;
    es->space[cur_word] |= value << (remaining - bits2store);
    remaining -= bits2store;
  }
  es->fillpos.cur_word = cur_word;
  es->fillpos.bitsleft = (unsigned int) remaining;
}

static inline GtBitsequence gt_editscript_space_get_next(const GtEditscript *es,
                                                         GtEditscriptPos *pos)
{
  GtBitsequence ret = 0;
  unsigned int shift;

  if (pos->bitsleft == 0) {
    gt_assert(pos->cur_word < es->fillpos.cur_word);
    pos->cur_word++;
    pos->bitsleft = (unsigned int) GT_INTWORDSIZE;
  }
  shift = pos->bitsleft;

  ret = es->space[pos->cur_word] << (GT_INTWORDSIZE - shift);
  if (shift < es->entry_size) {
    gt_assert(pos->cur_word < es->fillpos.cur_word);
    pos->cur_word++;
    ret |= es->space[pos->cur_word] >> (es->entry_size - shift);
    pos->bitsleft = (unsigned int) GT_INTWORDSIZE -
                          (es->entry_size - shift);
  }
  else {
    pos->bitsleft -= es->entry_size;
  }
  ret >>= GT_INTWORDSIZE - es->entry_size;
  return ret;
}

static inline void gt_editscript_space_add_length(GtEditscript *es,
                                                  GtBitsequence value)
{
  /* scheme: [1111][x][x][x][x] first store consecutive ones for each element
     needed, (like elias gamma) then store the value seperated to elements */
  unsigned int num_elems = 0,
               shift;
  GtBitsequence tmp = value;

  /* needs just one element not starting with 1, store it */
  if (tmp <= es->firstmask - 1) {
    gt_editscript_space_add_next(es, tmp);
  }
  else {
    while (tmp != 0) {
      num_elems++;
      tmp >>= es->entry_size;
    }
    /* number of one bits corresponds to num of elements needed */
    tmp = (GtBitsequence) (1 << num_elems) - 1;
    while ((tmp & es->firstmask) != 0) {
      /* add ~0 elements */
      gt_editscript_space_add_next(es, es->fullmask);
      tmp >>= es->entry_size;
    }
    if (tmp != 0) {
      /* move bits to front */
      while ((tmp & es->firstmask) == 0) {
        tmp <<= 1;
      }
    }
    /* either add remaining bits or one zero element to seperated */
    gt_editscript_space_add_next(es, tmp);

    /* store actual value */
    shift = (num_elems * es->entry_size);
    while (num_elems != 0) {
      num_elems--;
      shift -= es->entry_size;
      gt_editscript_space_add_next(es,
                                   (value >> (GtBitsequence) shift) &
                                     es->fullmask);
    }
  }
}

static inline GtBitsequence
gt_editscript_space_get_length(const GtEditscript *es,
                               GtEditscriptPos *pos,
                               GtUword *elems_used)
{

  GtBitsequence ret = 0,
                elem = 0,
                num_elems = 0;

  gt_assert(elems_used != NULL);
  gt_assert(es != NULL);
  gt_assert(pos != NULL);
  num_elems = elem = gt_editscript_space_get_next(es, pos);
  (*elems_used)++;
  /* first bit not set, the value itself. */
  if ((elem & es->firstmask) == 0) {
    return elem;
  }
  while (elem == es->fullmask) {
    elem = gt_editscript_space_get_next(es, pos);
    (*elems_used)++;
    if (elem != 0) {
      num_elems <<= es->entry_size;
      num_elems |= elem;
    }
  }
  /* move bits to right, max entry_size steps. */
  while ((num_elems & 1) == 0) {
    num_elems >>= 1;
  }
  while (num_elems != 0) {
    num_elems >>= 1;
    ret <<= es->entry_size;
    elem = gt_editscript_space_get_next(es, pos);
    (*elems_used)++;
    ret |= elem;
  }
  return ret;
}

void gt_editscript_add_match(GtEditscript *editscript)
{
  gt_assert(editscript);
  editscript->last = 0;
  editscript->trailing_matches++;
  editscript->ulen++;
  editscript->vlen++;
}

void gt_editscript_add_mismatch(GtEditscript *editscript, GtUchar c)
{
  gt_assert(editscript != NULL);
  gt_assert(c < (GtUchar) editscript->del || c == (GtUchar) WILDCARD);
  if (c == (GtUchar) WILDCARD) {
    c = (GtUchar) editscript->del - 1;
  }
  if (editscript->last != editscript->misdel) {
    editscript->op_count += 1;
    editscript->last = editscript->misdel;
    gt_editscript_space_add_next(editscript, editscript->misdel);
    gt_assert(editscript->vlen >= editscript->trailing_matches);
    gt_assert(editscript->ulen >= editscript->trailing_matches);
    gt_editscript_space_add_length(
                                  editscript,
                                  (GtBitsequence) editscript->trailing_matches);
    editscript->trailing_matches = 0;
  }
  editscript->ulen++;
  editscript->vlen++;
  gt_editscript_space_add_next(editscript, (GtBitsequence) c);
}

void gt_editscript_add_deletion(GtEditscript *editscript)
{
  gt_assert(editscript);
  if (editscript->last != editscript->misdel) {
    editscript->op_count += 1;
    editscript->last = editscript->misdel;
    gt_editscript_space_add_next(editscript, editscript->misdel);
    gt_assert(editscript->vlen >= editscript->trailing_matches);
    gt_assert(editscript->ulen >= editscript->trailing_matches);
    gt_editscript_space_add_length(
                                  editscript,
                                  (GtBitsequence) editscript->trailing_matches);
    editscript->trailing_matches = 0;
  }
  editscript->ulen++;
  gt_editscript_space_add_next(editscript, editscript->del);
}

void gt_editscript_add_insertion(GtEditscript *editscript, GtUchar c)
{
  gt_assert(editscript);
  gt_assert(c < (GtUchar) editscript->del || c == (GtUchar) WILDCARD);
  if (c == (GtUchar) WILDCARD) {
  c = (GtUchar) editscript->del - 1;
  }
  if (editscript->last != editscript->ins) {
    editscript->op_count += 1;
    editscript->last = editscript->ins;
    gt_editscript_space_add_next(editscript, editscript->ins);
    gt_assert(editscript->vlen >= editscript->trailing_matches);
    gt_assert(editscript->ulen >= editscript->trailing_matches);
    gt_editscript_space_add_length(
                                  editscript,
                                  (GtBitsequence) editscript->trailing_matches);
    editscript->trailing_matches = 0;
  }
  editscript->vlen++;
  gt_editscript_space_add_next(editscript, (GtBitsequence) c);
}

void gt_editscript_reset(GtEditscript *editscript)
{
  GtUword i;
  gt_assert(editscript);
  for (i = 0; i <= editscript->fillpos.cur_word; ++i) {
    editscript->space[i] = 0;
  }
  editscript->space_elems = 0;
  editscript->op_count = 0;
  gt_editscript_pos_reset(&(editscript->fillpos));
}

GtEditscript *gt_editscript_new_with_sequences(const GtEncseq *encseq,
                                               GtMultieoplist *multieops,
                                               GtUword start,
                                               GtReadmode dir)
{
  GtUword vlen, idx, meopidx, meoplen;
  GtEditscript *es;
  GtMultieop meop;
  GtUchar cchar;

  gt_assert(encseq != NULL && multieops != NULL);

  es = gt_editscript_new(gt_encseq_alphabet(encseq));

  vlen = gt_multieoplist_get_repins_length(multieops);
  meoplen = gt_multieoplist_get_length(multieops);

  for (idx = 0; idx < meoplen; idx++) {
    meop = gt_multieoplist_get_entry(multieops, idx);
    for (meopidx = 0; meopidx < meop.steps; meopidx++) {
      switch (meop.type) {
        case Match:
          vlen--;
          gt_editscript_add_match(es);
          break;
        case Replacement:
        case Mismatch:
          cchar = gt_encseq_get_encoded_char(encseq, start + vlen - 1, dir);
          vlen--;
          gt_editscript_add_mismatch(es, cchar);
          break;
        case Insertion:
          cchar = gt_encseq_get_encoded_char(encseq, start + vlen - 1, dir);
          vlen--;
          gt_editscript_add_insertion(es, cchar);
          break;
        case Deletion:
          gt_editscript_add_deletion(es);
          break;
        default:
          break;
      }
    }
  }
  gt_assert(vlen == 0);
  gt_assert(es->ulen == gt_multieoplist_get_repdel_length(multieops));
  gt_assert(es->vlen == gt_multieoplist_get_repins_length(multieops));
  gt_assert(es->trailing_matches <= es->vlen);
  gt_assert(es->trailing_matches <= es->ulen);
  es->space = gt_realloc(es->space,
                         sizeof (*(es->space)) * (es->fillpos.cur_word + 1));
  es->size = (size_t) es->fillpos.cur_word + (size_t) 1;
  return(es);
}

void gt_editscript_get_stats(const GtEditscript *editscript,
                             GtUword *match,
                             GtUword *mismatch,
                             GtUword *insertion,
                             GtUword *deletion)
{
  GtUword elems_served = 0;
  bool misdel = true; /* shut up scan-build */
  GtBitsequence elem;
  GtEditscriptPos getpos;

  gt_assert(editscript != NULL);

  gt_editscript_pos_reset(&getpos);
  *match = *mismatch = *insertion = *deletion = 0;
  if (editscript->op_count != 0) {
    elem = gt_editscript_space_get_next(editscript, &getpos);
    gt_assert(elem > editscript->del);
    elems_served++;
    if (elem == editscript->misdel) {
      misdel = true;
    }
    else if (elem == editscript->ins) {
      misdel = false;
    }
    (*match) += gt_editscript_space_get_length(editscript,
                                             &getpos,
                                             &elems_served);
    while (elems_served < editscript->space_elems) {
      gt_assert(editscript->vlen > *match + *mismatch + *insertion);
      gt_assert(editscript->ulen > *match + *mismatch + *deletion);
      elem = gt_editscript_space_get_next(editscript, &getpos);
      elems_served++;
      if (elem == editscript->misdel) {
        misdel = true;
        (*match) += gt_editscript_space_get_length(editscript,
                                                   &getpos,
                                                   &elems_served);
      }
      else if (elem == editscript->ins) {
        misdel = false;
        (*match) += gt_editscript_space_get_length(editscript,
                                                   &getpos,
                                                   &elems_served);
      }
      else {
        if (misdel) {
          if (elem == editscript->del)
            (*deletion)++;
          else
            (*mismatch)++;
        }
        else
          (*insertion)++;
      }
    }
    gt_assert(elems_served == editscript->space_elems);
  }
  (*match) += editscript->trailing_matches;
}

GtUword gt_editscript_get_ref_len(const GtEditscript *editscript)
{
  gt_assert(editscript != NULL);
  return editscript->ulen;
}

GtUword gt_editscript_get_target_len(const GtEditscript *editscript)
{
  gt_assert(editscript != NULL);
  return editscript->vlen;
}

GtUword gt_editscript_get_sequence(const GtEditscript *editscript,
                                   const GtEncseq *encseq,
                                   GtUword start,
                                   GtReadmode dir,
                                   GtUchar *buffer)
{
  GtUword i, j,
          vlen = 0,
          matchcount,
          elems_served = 0,
          uidx;
  GtBitsequence elem;
  GtEditscriptPos pos;

  gt_editscript_pos_reset(&pos);

  gt_assert(encseq != NULL && editscript != NULL);

  vlen = editscript->vlen;
  uidx = start + editscript->ulen;

  elem = gt_editscript_space_get_next(editscript, &pos);
  elems_served++;
  for (i = 0; i < editscript->op_count; ++i) {
    matchcount =
      (GtUword) gt_editscript_space_get_length(editscript, &pos,
                                                     &elems_served);
    gt_log_log("elem: "GT_WU", matches: "GT_WU,
               (GtUword) elem, matchcount);
    gt_log_log("current word: (" GT_WU ":%u " GT_WU,
               pos.cur_word,
               pos.bitsleft,
               (GtUword) editscript->space[pos.cur_word]);
    for (j = 0; j < matchcount; ++j) {
      buffer[vlen - 1] = gt_encseq_get_encoded_char(encseq, uidx - 1, dir);
      vlen--;
      uidx--;
    }
    if (elem == editscript->ins) {
      while (elems_served < editscript->space_elems) {
        elem = gt_editscript_space_get_next(editscript, &pos);
        elems_served++;
        gt_assert(elem != editscript->del);
        if (elem > editscript->del) {
          break;
        }
        buffer[vlen - 1] = elem == editscript->del - 1 ?
                           (GtUchar) WILDCARD :
                           (GtUchar) elem;
        vlen--;
      }
    }
    else {
      gt_assert(elem == editscript->misdel);
      while (elems_served < editscript->space_elems) {
        elem = gt_editscript_space_get_next(editscript, &pos);
        elems_served++;
        if (elem > editscript->del) {
          break;
        }
        if (elem < editscript->del) {
          buffer[vlen - 1] = elem == editscript->del - 1 ?
                             (GtUchar) WILDCARD :
                             (GtUchar) elem;
          vlen--;
        }
        uidx--;
      }
    }
  }
  for (j = 0; j < editscript->trailing_matches; ++j) {
    gt_assert(vlen != 0);
    buffer[vlen - 1] = gt_encseq_get_encoded_char(encseq, uidx - 1, dir);
    vlen--;
    uidx--;
  }
  gt_assert(uidx == start);
  gt_assert(vlen == 0);
  return(editscript->vlen);
}

size_t gt_editscript_size(GtEditscript *editscript)
{
  gt_assert(editscript != NULL);
  return sizeof (*editscript) +
         sizeof (*(editscript->space)) * editscript->size;
}

#define EDITSCRIPT_IO_ONE(elemptr, fp) \
  io_func(elemptr, sizeof (*(elemptr)), (size_t) 1, fp)

GtEditscript *gt_editscript_io(GtEditscript *editscript, FILE *fp,
                               EditscriptIOFunc io_func)
{
  if (editscript == NULL) {
    editscript = gt_calloc((size_t) 1, sizeof (GtEditscript));
    gt_editscript_pos_reset(&editscript->fillpos);
  }

  EDITSCRIPT_IO_ONE(&editscript->entry_size, fp);
  editscript->fullmask = (GtBitsequence) (1UL << editscript->entry_size) - 1UL;
  editscript->firstmask =
    (GtBitsequence) (1UL << (editscript->entry_size - 1UL));
  EDITSCRIPT_IO_ONE(&editscript->op_count, fp);
  EDITSCRIPT_IO_ONE(&editscript->trailing_matches, fp);
  EDITSCRIPT_IO_ONE(&editscript->ulen, fp);
  EDITSCRIPT_IO_ONE(&editscript->vlen, fp);
  EDITSCRIPT_IO_ONE(&editscript->del, fp);
  editscript->misdel = editscript->del + 1;
  editscript->ins = editscript->misdel + 1;
  EDITSCRIPT_IO_ONE(&editscript->space_elems, fp);

  EDITSCRIPT_IO_ONE(&editscript->fillpos.cur_word, fp);
  if (editscript->space == NULL) {
    editscript->size = (size_t) editscript->fillpos.cur_word + 1;
    editscript->space =
      gt_malloc(editscript->size * sizeof (*(editscript->space)));
  }
  io_func(editscript->space,
          sizeof (*(editscript->space)),
          (size_t) (editscript->fillpos.cur_word + 1),
          fp);

  return(editscript);
}

#define EDITSCRIPT_TEST_SEQLEN 20UL
#define gt_editscript_t_add_m()     \
  gt_multieoplist_add_match(meops); \
  gt_editscript_add_match(es)
#define gt_editscript_t_add_mm(POS)                             \
  gt_multieoplist_add_mismatch(meops);                          \
  gt_editscript_add_mismatch(es,                                \
                             gt_encseq_get_encoded_char(v, POS, \
                                                        GT_READMODE_FORWARD))
#define gt_editscript_t_add_d()        \
  gt_multieoplist_add_deletion(meops); \
  gt_editscript_add_deletion(es)
#define gt_editscript_t_add_i(POS)                               \
  gt_multieoplist_add_insertion(meops);                          \
  gt_editscript_add_insertion(es,                                \
                              gt_encseq_get_encoded_char(v, POS, \
                                                         GT_READMODE_FORWARD))

int gt_editscript_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtUword length, i, max;
  GtMultieoplist *meops = gt_multieoplist_new();
  GtEncseq *v = NULL, *u = NULL;
  GtAlphabet *dna = gt_alphabet_new_dna();
  GtEncseqBuilder *esb = gt_encseq_builder_new(dna);
  GtUchar buffer[EDITSCRIPT_TEST_SEQLEN] = {0};
  GtEditscriptPos get_pos;
  GtBitsequence soll, ist;
  /*
    u: aaacccg-ggttt--acgtacga
       || | || ||| |  |  |  ||
    v: aatc-cggggtatcga--tgtga
    meops(reverse): MmMMMIMMDMmMM
   */
  const char *seq1 = "AAACCCGGGTTTACGTACGA",
             *seq2 = "AATCCGGGGTATCGATGTGA";
  GtEditscript *es = gt_editscript_new(dna);

  gt_error_check(err);

  gt_editscript_pos_reset(&get_pos);
  max = (GtUword) ((1 << es->entry_size) - 1);
  for (i = 0; !had_err && i <= max; ++i) {
    gt_editscript_space_add_next(es, (GtBitsequence) i);
    ist = gt_editscript_space_get_next(es, &get_pos);
    gt_ensure(i == (GtUword) ist);
  }
  if (!had_err) {
    GtUword elements = 0;
    gt_editscript_pos_reset(&get_pos);
    gt_editscript_reset(es);
    soll = (GtBitsequence) ((1 << es->entry_size) - 1);
    gt_editscript_space_add_length(es, soll);
    ist = gt_editscript_space_get_length(es, &get_pos, &elements);
    gt_ensure(soll == ist);
    for (i = 0; !had_err && i <= max; ++i) {
      soll = (GtBitsequence) (i | (i << es->entry_size));
      gt_editscript_space_add_length(es, soll);
      ist = gt_editscript_space_get_length(es, &get_pos, &elements);
      gt_ensure(soll == ist);
    }
    gt_ensure(elements == es->space_elems);
    gt_editscript_pos_reset(&es->fillpos);
  }

  if (!had_err) {
    gt_editscript_reset(es);
    gt_encseq_builder_add_cstr(esb, seq1, EDITSCRIPT_TEST_SEQLEN, "u");
    u = gt_encseq_builder_build(esb, err);
    gt_encseq_builder_reset(esb);
    gt_encseq_builder_add_cstr(esb, seq2, EDITSCRIPT_TEST_SEQLEN, "v");
    v = gt_encseq_builder_build(esb, err);
    gt_encseq_builder_delete(esb);

    gt_editscript_t_add_m();
    gt_editscript_t_add_m();
    gt_editscript_t_add_mm(17UL);
    gt_editscript_t_add_mm(16UL);
    gt_editscript_t_add_m();
    gt_editscript_t_add_d();
    gt_editscript_t_add_d();
    gt_editscript_t_add_m();
    gt_editscript_t_add_i(13UL);
    gt_editscript_t_add_i(12UL);
    gt_editscript_t_add_m();
    gt_editscript_t_add_mm(10UL);
    gt_editscript_t_add_m();
    gt_editscript_t_add_m();
    gt_editscript_t_add_m();
    gt_editscript_t_add_i(6UL);
    gt_editscript_t_add_m();
    gt_editscript_t_add_m();
    gt_editscript_t_add_d();
    gt_editscript_t_add_m();
    gt_editscript_t_add_mm(2UL);
    gt_editscript_t_add_m();
    gt_editscript_t_add_m();

    length = gt_editscript_get_sequence(es, u, 0, GT_READMODE_FORWARD, buffer);
    gt_ensure(length == EDITSCRIPT_TEST_SEQLEN);
  }
  for (i = 0; !had_err && i < EDITSCRIPT_TEST_SEQLEN; ++i) {
    gt_ensure(buffer[i] == gt_encseq_get_encoded_char(v, (GtUword) i,
                                                      GT_READMODE_FORWARD));
  }
  gt_editscript_delete(es);
  if (!had_err) {
    es = gt_editscript_new_with_sequences(v, meops, 0, GT_READMODE_FORWARD);
    length = gt_editscript_get_sequence(es, u, 0, GT_READMODE_FORWARD, buffer);
    gt_ensure(length == EDITSCRIPT_TEST_SEQLEN);
    for (i = 0; !had_err && i < EDITSCRIPT_TEST_SEQLEN; ++i) {
      gt_ensure(buffer[i] == gt_encseq_get_encoded_char(v, (GtUword) i,
                                                        GT_READMODE_FORWARD));
    }
  }
  if (!had_err) {
    GtUword match, mismatch, del, ins;
    gt_editscript_get_stats(es, &match, &mismatch, &ins, &del);
    gt_ensure(match == 13UL);
    gt_ensure(mismatch == 4UL);
    gt_ensure(ins == 3UL);
    gt_ensure(del == 3UL);
  }
  gt_editscript_delete(es);
  gt_alphabet_delete(dna);
  gt_encseq_delete(u);
  gt_encseq_delete(v);
  gt_multieoplist_delete(meops);
  return had_err;
}
