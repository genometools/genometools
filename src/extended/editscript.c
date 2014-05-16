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
#include <string.h>

#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/encseq_api.h"
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
#include "match/xdrop.h"

#define GT_EDITSCRIPT_SHOW_BITS(es)                 \
do {                                                \
  char bitbuffer[65] = {0};                         \
  gt_bitsequence_tostring(bitbuffer, elem);         \
  gt_log_log("%s", bitbuffer/*  + (64-es->entry_size) */);\
} while (false)

typedef struct GtEditscriptPos {
  unsigned int cur_word,
               bitsleft;
}GtEditscriptPos;

struct GtEditscript {
  GtBitsequence      *space;
  unsigned int        size,
                      num_elems;
  unsigned short int  trailing_matches;
  unsigned char       del,
                      entry_size;
};

struct GtEditscriptBuilder {
  GtEditscript    *es;
  GtEditscriptPos  fillpos;
  unsigned char    last_op;
};

static inline void gt_editscript_pos_reset(GtEditscriptPos *pos)
{
  pos->cur_word = 0;
  pos->bitsleft = (unsigned int) GT_INTWORDSIZE;
}

#define GT_EDITSCRIPT_FULLMASK(es) ((GtBitsequence) \
                                ((GtUword) 1 << es->entry_size) - (GtUword) 1)
#define GT_EDITSCRIPT_FIRSTMASK(es) ((GtBitsequence) \
                                 ((GtUword) 1 << (es->entry_size - \
                                                  (GtUword) 1)))
#define GT_EDITSCRIPT_MISDEL_SYM(es) ((GtBitsequence) es->del + 1)
#define GT_EDITSCRIPT_INS_SYM(es) ((GtBitsequence) es->del + 2)

GtEditscript *gt_editscript_new(GtAlphabet *alphabet)
{
  unsigned int alphabet_size;
  GtEditscript *es;
  es = gt_calloc((size_t) 1, sizeof (*es));
  alphabet_size = gt_alphabet_size(alphabet);
  es->entry_size =
    gt_determinebitspervalue((GtUword) alphabet_size + 3);
  gt_assert(es->entry_size <= (unsigned char) (sizeof (unsigned char) * 8));
  es->size = 0U;
  es->space = NULL;
  es->del = (unsigned char) alphabet_size;
  es->num_elems = 0;
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
                                                GtEditscriptPos *fillpos,
                                                GtBitsequence elem)
{
  unsigned int cur_word = fillpos->cur_word,
               remaining = fillpos->bitsleft,
               bits2store = (unsigned int) es->entry_size;
  gt_assert(elem <= GT_EDITSCRIPT_FULLMASK(es));
  if (es->size == 0) {
    es->size = 1U;
    es->space = gt_malloc(es->size * sizeof (*es->space));
    es->space[0] = 0;
  }
  es->num_elems++;
  if (remaining > bits2store) {
    es->space[cur_word] |= elem << (remaining - bits2store);
    remaining -= bits2store;
  }
  else {
    if (cur_word + 1 >= es->size) {
      es->size = GT_MULT2(es->size);
      es->space = gt_realloc(es->space, es->size * sizeof (*(es->space)));
    }
    gt_assert(es->space != NULL);
    es->space[cur_word] |= elem >> (bits2store - remaining);
    cur_word++;
    es->space[cur_word] = 0;
    bits2store -= remaining;
    remaining = (unsigned int) GT_INTWORDSIZE;
    /* this has to be checked because (uint64_t << 64) is not defined */
    if (bits2store != 0) {
      es->space[cur_word] |= elem << (remaining - bits2store);
      remaining -= bits2store;
    }
  }
  fillpos->cur_word = cur_word;
  fillpos->bitsleft = remaining;
}

static inline GtBitsequence gt_editscript_space_get_next(const GtEditscript *es,
                                                         GtEditscriptPos *pos)
{
  GtBitsequence elem = 0;
  unsigned int shift;

  gt_assert(es->size != 0);
  if (pos->bitsleft == 0) {
    pos->cur_word++;
    pos->bitsleft = (unsigned int) GT_INTWORDSIZE;
  }
  shift = pos->bitsleft;

  elem = es->space[pos->cur_word] << (GT_INTWORDSIZE - shift);
  if (shift < es->entry_size) {
    pos->cur_word++;
    elem |= es->space[pos->cur_word] >> (GtBitsequence) shift;
    pos->bitsleft = (unsigned int) (GT_INTWORDSIZE - (es->entry_size - shift));
  }
  else {
    pos->bitsleft -= es->entry_size;
  }
  elem >>= GT_INTWORDSIZE - es->entry_size;
  return elem;
}

static inline void gt_editscript_space_add_length(GtEditscript *es,
                                                  GtEditscriptPos *fillpos,
                                                  GtBitsequence value)
{
  /* scheme: [1111][x][x][x][x] first store consecutive ones for each element
     needed, (like elias gamma) then store the value seperated to elements */
  unsigned int num_elems = 0,
               shift;
  GtBitsequence tmp = value;

  /* needs just one element not starting with 1, store it */
  if (tmp <= GT_EDITSCRIPT_FIRSTMASK(es) - 1) {
    gt_editscript_space_add_next(es, fillpos, tmp);
  }
  else {
    while (tmp != 0) {
      num_elems++;
      tmp >>= es->entry_size;
    }
    /* number of one bits corresponds to num of elements needed */
    tmp = (GtBitsequence) (1 << num_elems) - 1;
    while ((tmp & GT_EDITSCRIPT_FIRSTMASK(es)) != 0) {
      /* add ~0 elements */
      gt_editscript_space_add_next(es, fillpos, GT_EDITSCRIPT_FULLMASK(es));
      tmp >>= es->entry_size;
    }
    if (tmp != 0) {
      /* move bits to front */
      while ((tmp & GT_EDITSCRIPT_FIRSTMASK(es)) == 0) {
        tmp <<= 1;
      }
    }
    /* either add remaining bits or one zero element to seperated */
    gt_editscript_space_add_next(es, fillpos, tmp);

    /* store actual value */
    shift = (num_elems * es->entry_size);
    while (num_elems != 0) {
      num_elems--;
      shift -= es->entry_size;
      gt_editscript_space_add_next(es, fillpos,
                                   (value >> (GtBitsequence) shift) &
                                     GT_EDITSCRIPT_FULLMASK(es));
    }
  }
}

static inline GtBitsequence
gt_editscript_space_get_length(const GtEditscript *es,
                               GtEditscriptPos *pos,
                               unsigned int *elems_used)
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
  if ((elem & GT_EDITSCRIPT_FIRSTMASK(es)) == 0) {
    return elem;
  }
  while (elem == GT_EDITSCRIPT_FULLMASK(es)) {
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

void gt_editscript_reset(GtEditscript *editscript)
{
  gt_assert(editscript);
  if (editscript->space != NULL) {
    editscript->space[0] = 0;
  }
  editscript->num_elems = 0;
}

GtEditscript *gt_editscript_new_with_sequences(const GtEncseq *encseq,
                                               GtMultieoplist *multieops,
                                               GtUword start,
                                               GtReadmode dir)
{
  GtUword vlen, idx, meopidx, meoplen;
  GtEditscript *es;
  GtEditscriptBuilder es_b;
  GtMultieop meop;
  GtUchar cchar;

  gt_assert(encseq != NULL && multieops != NULL);
  es = gt_editscript_new(gt_encseq_alphabet(encseq));
  gt_editscript_builder_reset(&es_b, es);

  vlen = 0;

  meoplen = gt_multieoplist_get_num_entries(multieops);

  for (idx = meoplen; idx != 0; idx--) {
    meop = gt_multieoplist_get_entry(multieops, idx - 1);
    for (meopidx = 0; meopidx < meop.steps; meopidx++) {
      switch (meop.type) {
        case Match:
          vlen++;
          gt_editscript_builder_add_match(&es_b);
          break;
        case Replacement:
        case Mismatch:
          cchar = gt_encseq_get_encoded_char(encseq, start + vlen, dir);
          vlen++;
          gt_editscript_builder_add_mismatch(&es_b, cchar);
          break;
        case Insertion:
          cchar = gt_encseq_get_encoded_char(encseq, start + vlen, dir);
          vlen++;
          gt_editscript_builder_add_insertion(&es_b, cchar);
          break;
        case Deletion:
          gt_editscript_builder_add_deletion(&es_b);
          break;
        default:
          break;
      }
    }
  }
  gt_assert(vlen == gt_multieoplist_get_repins_length(multieops));
  gt_assert(vlen <= gt_encseq_total_length(encseq));
  if (es->num_elems != 0) {
    es->space = gt_realloc(es->space,
                           sizeof (*(es->space)) * (es_b.fillpos.cur_word + 1));
    es->size = es_b.fillpos.cur_word + 1;
  }
  else {
    es->size = 0;
    gt_free(es->space);
    es->space = NULL;
  }
  return(es);
}

void gt_editscript_get_stats(const GtEditscript *editscript,
                             GtUword *match,
                             GtUword *mismatch,
                             GtUword *insertion,
                             GtUword *deletion)
{
  unsigned int elems_served = 0;
  bool misdel = true; /* shut up scan-build */
  GtBitsequence elem;
  GtEditscriptPos getpos;

  gt_assert(editscript != NULL);

  gt_editscript_pos_reset(&getpos);
  *match = *mismatch = *insertion = *deletion = 0;
  if (editscript->num_elems != 0) {
    elem = gt_editscript_space_get_next(editscript, &getpos);
    gt_assert(editscript->del < (unsigned char) elem);
    elems_served++;
    if (elem == GT_EDITSCRIPT_MISDEL_SYM(editscript)) {
      misdel = true;
    }
    else if (elem == GT_EDITSCRIPT_INS_SYM(editscript)) {
      misdel = false;
    }
    (*match) += gt_editscript_space_get_length(editscript,
                                               &getpos,
                                               &elems_served);
    while (elems_served < editscript->num_elems) {
      elem = gt_editscript_space_get_next(editscript, &getpos);
      elems_served++;
      if (elem == GT_EDITSCRIPT_MISDEL_SYM(editscript)) {
        misdel = true;
        (*match) += gt_editscript_space_get_length(editscript,
                                                   &getpos,
                                                   &elems_served);
      }
      else if (elem == GT_EDITSCRIPT_INS_SYM(editscript)) {
        misdel = false;
        (*match) += gt_editscript_space_get_length(editscript,
                                                   &getpos,
                                                   &elems_served);
      }
      else {
        if (misdel) {
          if (elem == (GtBitsequence) editscript->del)
            (*deletion)++;
          else
            (*mismatch)++;
        }
        else
          (*insertion)++;
      }
    }
    gt_assert(elems_served == editscript->num_elems);
  }
  (*match) += editscript->trailing_matches;
}

GtUword gt_editscript_get_target_len(const GtEditscript *editscript)
{
  GtEditscriptPos pos;
  GtBitsequence elem;
  GtUword length = 0;
  unsigned int served = 0;

  gt_assert(editscript != NULL);
  gt_editscript_pos_reset(&pos);
  if (editscript->num_elems != 0) {
    (void) gt_editscript_space_get_next(editscript, &pos);
    served++;
    while (served < editscript->num_elems) {
      length += gt_editscript_space_get_length(editscript, &pos, &served);
      while (served < editscript->num_elems) {
        elem = gt_editscript_space_get_next(editscript, &pos);
        served++;
        if (elem > (GtBitsequence) editscript->del)
          break;
        if (elem < (GtBitsequence) editscript->del)
          length++;
      }
    }
  }
  return length + editscript->trailing_matches;
}

GtUword gt_editscript_get_sequence(const GtEditscript *editscript,
                                   const GtEncseq *encseq,
                                   GtUword start,
                                   GtReadmode dir,
                                   GtUchar *buffer)
{
  GtEditscriptPos pos;
  GtBitsequence elem = 0;
  GtUword vidx = 0,
          uidx;
  unsigned int j,elems_served = 0,
               matchcount;

  gt_editscript_pos_reset(&pos);

  gt_assert(encseq != NULL && editscript != NULL);

  uidx = start;

  if (editscript->num_elems != 0) {
    bool mismatch_or_deletion;
    elem = gt_editscript_space_get_next(editscript, &pos);
    elems_served++;
    while (elems_served < editscript->num_elems) {
      gt_assert(elem <= GT_EDITSCRIPT_INS_SYM(editscript));
      matchcount =
        (unsigned int) gt_editscript_space_get_length(editscript, &pos,
                                                      &elems_served);
      for (j = 0; j < matchcount; ++j) {
        buffer[vidx] = gt_encseq_get_encoded_char(encseq, uidx, dir);
        vidx++;
        uidx++;
      }
      mismatch_or_deletion = elem == GT_EDITSCRIPT_MISDEL_SYM(editscript);
      while (elems_served < editscript->num_elems) {
        elem = gt_editscript_space_get_next(editscript, &pos);
        elems_served++;
        if (elem > (GtBitsequence) editscript->del) {
          break;
        }
        if (elem < (GtBitsequence) editscript->del) {
          buffer[vidx] = (elem == (GtBitsequence) editscript->del - 1) ?
            (GtUchar) WILDCARD :
            (GtUchar) elem;
          vidx++;
        }
        if (mismatch_or_deletion) uidx++;
      }
    }
  }
  for (j = 0; j < (unsigned int) editscript->trailing_matches; ++j) {
    buffer[vidx] = gt_encseq_get_encoded_char(encseq, uidx, dir);
    vidx++;
    uidx++;
  }
  return(vidx);
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
  GtUword bits;
  if (editscript == NULL) {
    editscript = gt_calloc((size_t) 1, sizeof (GtEditscript));
  }

  EDITSCRIPT_IO_ONE(&editscript->entry_size, fp);
  EDITSCRIPT_IO_ONE(&editscript->trailing_matches, fp);
  EDITSCRIPT_IO_ONE(&editscript->del, fp);
  EDITSCRIPT_IO_ONE(&editscript->num_elems, fp);
  if (editscript->num_elems != 0) {
    bits = (GtUword) editscript->num_elems * (GtUword) editscript->entry_size;
    editscript->size = (unsigned int) (bits / (GT_INTWORDSIZE));
    if (bits % GT_INTWORDSIZE) {
      editscript->size++;
    }

    if (editscript->space == NULL) {
      editscript->space =
        gt_malloc(editscript->size * sizeof (*(editscript->space)));
    }
    io_func(editscript->space,
            sizeof (*(editscript->space)),
            (size_t) (editscript->size),
            fp);
  }

  return(editscript);
}

void gt_editscript_show(GtEditscript *editscript, GtAlphabet *alphabet)
{
  GtEditscriptPos pos;
  GtBitsequence elem;
  unsigned int matchcount,
               elems_served = 0;

  gt_editscript_pos_reset(&pos);

  printf("|");
  if (editscript->num_elems != 0) {
    elem = gt_editscript_space_get_next(editscript, &pos);
    elems_served++;
    while (elems_served < editscript->num_elems) {
      gt_assert(elem <= GT_EDITSCRIPT_INS_SYM(editscript));
      matchcount =
        (unsigned int) gt_editscript_space_get_length(editscript, &pos,
                                                      &elems_served);
      if (elem == GT_EDITSCRIPT_INS_SYM(editscript)) {
        if (matchcount != 0) {
          printf("M(%u)|Ins:|", matchcount);
        }
        else {
          printf("Ins:|");
        }
        while (elems_served < editscript->num_elems) {
          elem = gt_editscript_space_get_next(editscript, &pos);
          elems_served++;
          gt_assert(elem != (GtBitsequence) editscript->del);
          elem = (elem == (GtBitsequence) editscript->del - 1) ?
            (GtBitsequence) WILDCARD :
            elem;
          if (elem > (GtBitsequence) editscript->del) {
            break;
          }
          printf("%c|", (int) gt_alphabet_pretty_symbol(alphabet,
                                                        (unsigned int) elem));
        }
      }
      else {
        gt_assert(elem == GT_EDITSCRIPT_MISDEL_SYM(editscript));
        if (matchcount != 0) {
          printf("M(%u)|Misdel:|", matchcount);
        }
        else {
          printf("Misdel:|");
        }
        while (elems_served < editscript->num_elems) {
          elem = gt_editscript_space_get_next(editscript, &pos);
          elems_served++;
          if (elem > (GtBitsequence) editscript->del) {
            break;
          }
          if (elem < (GtBitsequence) editscript->del) {
            elem = (elem == (GtBitsequence) editscript->del - 1) ?
              (GtBitsequence) WILDCARD :
              elem;
            printf("%c|", (int) gt_alphabet_pretty_symbol(alphabet,
                                                          (unsigned int) elem));
          }
          else
            printf("-|");
        }
      }
    }
  }
  if (editscript->trailing_matches != 0) {
    printf("M(%hu)|\n", editscript->trailing_matches);
  }
  else
    printf("\n");
}

GtEditscriptBuilder *gt_editscript_builder_new(GtEditscript *editscript)
{
  GtEditscriptBuilder *es_b = gt_malloc(sizeof (*es_b));
  gt_editscript_pos_reset(&es_b->fillpos);
  es_b->es = editscript;
  es_b->last_op = 0;
  return es_b;
}

void gt_editscript_builder_reset(GtEditscriptBuilder *es_builder,
                                 GtEditscript *editscript)
{
  gt_assert(es_builder);
  gt_assert(editscript);
  gt_editscript_reset(editscript);
  gt_editscript_pos_reset(&es_builder->fillpos);
  es_builder->es = editscript;
  es_builder->last_op = 0;
}

void gt_editscript_builder_add_match(GtEditscriptBuilder *es_builder)
{
  gt_assert(es_builder);
  es_builder->last_op = 0;
  es_builder->es->trailing_matches++;
}

void gt_editscript_builder_add_mismatch(GtEditscriptBuilder *es_builder,
                                        GtUchar c)
{
  GtEditscript *es;
  gt_assert(es_builder != NULL);
  es = es_builder->es;
  gt_assert(c < (GtUchar) es->del || c == (GtUchar) WILDCARD);
  if (c == (GtUchar) WILDCARD) {
    c = (GtUchar) es->del - 1;
  }
  if (es_builder->last_op != (unsigned char) GT_EDITSCRIPT_MISDEL_SYM(es)) {
    es_builder->last_op = (unsigned char) GT_EDITSCRIPT_MISDEL_SYM(es);
    gt_editscript_space_add_next(es, &es_builder->fillpos,
                                 GT_EDITSCRIPT_MISDEL_SYM(es));
    gt_editscript_space_add_length(es, &es_builder->fillpos,
                                   (GtBitsequence) es->trailing_matches);
    es->trailing_matches = 0;
  }
  gt_editscript_space_add_next(es, &es_builder->fillpos, (GtBitsequence) c);
}

void gt_editscript_builder_add_deletion(GtEditscriptBuilder *es_builder)
{
  GtEditscript *es;
  gt_assert(es_builder);
  es = es_builder->es;

  if (es_builder->last_op != (unsigned char) GT_EDITSCRIPT_MISDEL_SYM(es)) {
    es_builder->last_op = (unsigned char) GT_EDITSCRIPT_MISDEL_SYM(es);
    gt_editscript_space_add_next(es, &es_builder->fillpos,
                                 GT_EDITSCRIPT_MISDEL_SYM(es));
    gt_editscript_space_add_length(es, &es_builder->fillpos,
                                   (GtBitsequence) es->trailing_matches);
    es->trailing_matches = 0;
  }
  gt_editscript_space_add_next(es, &es_builder->fillpos,
                               (GtBitsequence) es->del);
}

void gt_editscript_builder_add_insertion(GtEditscriptBuilder *es_builder,
                                         GtUchar c)
{
  GtEditscript *es;
  gt_assert(es_builder);
  es = es_builder->es;

  gt_assert(c < (GtUchar) es->del || c == (GtUchar) WILDCARD);
  if (c == (GtUchar) WILDCARD) {
    c = (GtUchar) es->del - 1;
  }
  if (es_builder->last_op != (unsigned char) GT_EDITSCRIPT_INS_SYM(es)) {
    es_builder->last_op = (unsigned char) GT_EDITSCRIPT_INS_SYM(es);
    gt_editscript_space_add_next(es, &es_builder->fillpos,
                                 GT_EDITSCRIPT_INS_SYM(es));
    gt_editscript_space_add_length(es, &es_builder->fillpos,
                                   (GtBitsequence) es->trailing_matches);
    es->trailing_matches = 0;
  }
  gt_editscript_space_add_next(es, &es_builder->fillpos, (GtBitsequence) c);
}

#define EDITSCRIPT_TEST_SEQLEN 24UL
#define EDITSCRIPT_2_TIMES(X)\
  X;                         \
  X
#define EDITSCRIPT_3_TIMES(X)\
  EDITSCRIPT_2_TIMES(X);     \
  X

int gt_editscript_unit_test(GT_UNUSED GtError *err)
{
  GtEditscriptBuilder es_b;
  GtEditscriptPos get_pos, fill_pos;
  GtUchar buffer[EDITSCRIPT_TEST_SEQLEN] = {0};
  GtMultieoplist *meops = gt_multieoplist_new();
  GtEncseq *v = NULL, *u = NULL;
  GtAlphabet *dna = gt_alphabet_new_dna();
  GtEncseqBuilder *esb = gt_encseq_builder_new(dna);
  GtBitsequence soll, ist;
  GtUword length, i, max;
  int had_err = 0;
  GtUchar csoll;
  /*
u: aaacccg-ggttt--acgtacgnang-a
   || | || ||| |  |  |  | | | |
v: aatc-cggggtatcga--tgtgna-gna
meops(reverse) MIMDMmMmmMDDMIIMmMMMIMMDMmMM
*/
  const char *seq1 = "AAACCCGGGTTTACGTACGNANGA",
             *seq2 = "AATCCGGGGTATCGATGTGNAGNA";
  GtEditscript *es = gt_editscript_new(dna);
  gt_log_log("size of empty: " GT_WU, (GtUword) gt_editscript_size(es));

  gt_error_check(err);

  gt_editscript_pos_reset(&get_pos);
  gt_editscript_pos_reset(&fill_pos);
  max = (GtUword) ((1 << es->entry_size) - 1);
  for (i = 0; !had_err && i <= max; ++i) {
    gt_editscript_space_add_next(es, &fill_pos, (GtBitsequence) i);
    ist = gt_editscript_space_get_next(es, &get_pos);
    gt_log_log("expected: " GT_WU, i);
    gt_log_log(" reality: " GT_WU, (GtUword) ist);
    gt_ensure(i == (GtUword) ist);
  }
  if (!had_err) {
    unsigned int served = 0;
    gt_editscript_pos_reset(&get_pos);
    gt_editscript_pos_reset(&fill_pos);
    gt_editscript_reset(es);
    soll = (GtBitsequence) ((1 << es->entry_size) - 1);
    gt_editscript_space_add_length(es, &fill_pos, soll);
    ist = gt_editscript_space_get_length(es, &get_pos, &served);
    gt_ensure(soll == ist);
    for (i = 0; !had_err && i <= max; ++i) {
      soll = (GtBitsequence) (i | (i << es->entry_size));
      gt_editscript_space_add_length(es, &fill_pos, soll);
      ist = gt_editscript_space_get_length(es, &get_pos, &served);
      gt_ensure(soll == ist);
    }
    gt_ensure(served == es->num_elems);
    gt_editscript_pos_reset(&fill_pos);
  }

  if (!had_err) {
    gt_encseq_builder_add_cstr(esb, seq1, EDITSCRIPT_TEST_SEQLEN, "u");
    u = gt_encseq_builder_build(esb, err);
    gt_encseq_builder_reset(esb);
    gt_encseq_builder_add_cstr(esb, seq2, EDITSCRIPT_TEST_SEQLEN, "v");
    v = gt_encseq_builder_build(esb, err);
    gt_encseq_builder_delete(esb);

    gt_editscript_builder_reset(&es_b, es);
    gt_multieoplist_add_match(meops);
    gt_multieoplist_add_insertion(meops);
    gt_multieoplist_add_match(meops);
    gt_multieoplist_add_deletion(meops);
    EDITSCRIPT_2_TIMES(gt_multieoplist_add_match(meops);
                       gt_multieoplist_add_mismatch(meops));
    gt_multieoplist_add_mismatch(meops);
    gt_multieoplist_add_match(meops);
    EDITSCRIPT_2_TIMES(gt_multieoplist_add_deletion(meops));
    gt_multieoplist_add_match(meops);
    EDITSCRIPT_2_TIMES(gt_multieoplist_add_insertion(meops));
    gt_multieoplist_add_match(meops);
    gt_multieoplist_add_mismatch(meops);
    EDITSCRIPT_3_TIMES(gt_multieoplist_add_match(meops));
    gt_multieoplist_add_insertion(meops);
    EDITSCRIPT_2_TIMES(gt_multieoplist_add_match(meops));
    gt_multieoplist_add_deletion(meops);
    gt_multieoplist_add_match(meops);
    gt_multieoplist_add_mismatch(meops);
    EDITSCRIPT_2_TIMES(gt_multieoplist_add_match(meops));

    EDITSCRIPT_2_TIMES(gt_editscript_builder_add_match(&es_b));
    csoll = gt_encseq_get_encoded_char(v, 2UL, GT_READMODE_FORWARD);
    gt_log_log("add mmchar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_mismatch(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);
    gt_editscript_builder_add_deletion(&es_b);
    EDITSCRIPT_2_TIMES(gt_editscript_builder_add_match(&es_b));
    csoll = gt_encseq_get_encoded_char(v, 6UL, GT_READMODE_FORWARD);
    gt_log_log("add inschar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_insertion(&es_b, csoll);
    EDITSCRIPT_3_TIMES(gt_editscript_builder_add_match(&es_b));
    csoll = gt_encseq_get_encoded_char(v, 10UL, GT_READMODE_FORWARD);
    gt_log_log("add mmchar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_mismatch(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);
    csoll = gt_encseq_get_encoded_char(v, 12UL, GT_READMODE_FORWARD);
    gt_log_log("add inschar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_insertion(&es_b, csoll);
    csoll = gt_encseq_get_encoded_char(v, 13UL, GT_READMODE_FORWARD);
    gt_log_log("add inschar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_insertion(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);
    EDITSCRIPT_2_TIMES(gt_editscript_builder_add_deletion(&es_b));
    gt_editscript_builder_add_match(&es_b);
    csoll = gt_encseq_get_encoded_char(v, 16UL, GT_READMODE_FORWARD);
    gt_log_log("add mmchar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_mismatch(&es_b, csoll);
    csoll = gt_encseq_get_encoded_char(v, 17UL, GT_READMODE_FORWARD);
    gt_log_log("add mmchar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_mismatch(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);
    csoll = gt_encseq_get_encoded_char(v, 19UL, GT_READMODE_FORWARD);
    gt_log_log("add mmchar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_mismatch(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);
    gt_editscript_builder_add_deletion(&es_b);
    gt_editscript_builder_add_match(&es_b);
    csoll = gt_encseq_get_encoded_char(v, 22UL, GT_READMODE_FORWARD);
    gt_log_log("add inschar: %c", gt_alphabet_decode(dna, csoll));
    gt_editscript_builder_add_insertion(&es_b, csoll);
    gt_editscript_builder_add_match(&es_b);

    length = gt_editscript_get_sequence(es, u, 0, GT_READMODE_FORWARD, buffer);
    gt_ensure(length == EDITSCRIPT_TEST_SEQLEN);
    gt_ensure(gt_editscript_get_target_len(es) == length);
    if (gt_log_enabled()) {
      gt_editscript_show(es, dna);
    }
  }
  for (i = 0; !had_err && i < EDITSCRIPT_TEST_SEQLEN; ++i) {
    csoll = gt_encseq_get_encoded_char(v, (GtUword) i,
                                       GT_READMODE_FORWARD);
    gt_log_log("i: "GT_WU " ist: %c, csoll: %c", i,
               gt_alphabet_decode(dna, buffer[i]),
               gt_alphabet_decode(dna, csoll));
    gt_ensure(buffer[i] == csoll);
  }
  if (!had_err) {
    gt_editscript_delete(es);
    es = gt_editscript_new_with_sequences(v, meops, 0, GT_READMODE_FORWARD);
    length = gt_editscript_get_sequence(es, u, 0, GT_READMODE_FORWARD, buffer);
    gt_ensure(length == EDITSCRIPT_TEST_SEQLEN);
    gt_ensure(gt_editscript_get_target_len(es) == length);
    if (gt_log_enabled()) {
      gt_editscript_show(es, dna);
    }
    for (i = 0; !had_err && i < EDITSCRIPT_TEST_SEQLEN; ++i) {
      csoll = gt_encseq_get_encoded_char(v, (GtUword) i,
                                         GT_READMODE_FORWARD);
      gt_log_log("i: "GT_WU " ist: %c, csoll: %c", i,
                 gt_alphabet_decode(dna, buffer[i]),
                 gt_alphabet_decode(dna, csoll));
      gt_ensure(buffer[i] == gt_encseq_get_encoded_char(v, (GtUword) i,
                                                        GT_READMODE_FORWARD));
    }
  }
  if (!had_err) {
    GtUword match, mismatch, del, ins;
    gt_editscript_get_stats(es, &match, &mismatch, &ins, &del);
    gt_ensure(match == 15UL);
    gt_ensure(mismatch == 5UL);
    gt_ensure(ins == 4UL);
    gt_ensure(del == 4UL);
  }
  gt_editscript_delete(es);
  gt_alphabet_delete(dna);
  gt_encseq_delete(u);
  gt_encseq_delete(v);
  gt_multieoplist_delete(meops);

  return had_err;
}
