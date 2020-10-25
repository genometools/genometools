/*
  Copyright (c) 2017  Niklas Niehus <niklas.niehus@live.com>
  Copyright (c) 2017  Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/arraydef_api.h"
#include "core/sequence_buffer_dust.h"
#include "core/sequence_buffer.h"
#include "core/ma_api.h"
#include "core/minmax_api.h"
#include "core/chardef_api.h"

typedef unsigned char GtInl_Queueelem;
#include "match/queue-inline.h"

#define GT_MAXTRIPLETVALUE 63
#define GT_ARRAYEXTENDSIZE 100
#define GT_TREATNASRANDOM 0

typedef struct {
  GtUword begin,
          end;
} GtDustRange;

GT_DECLAREARRAYSTRUCT(GtDustRange); /* GtArrayDustRange */

typedef struct {
  unsigned char val,
                orig;
  GtUword mask_length;
  GtUword next_mask;
  float max_score;
} GtDustBuffer;

struct GtDustMasker {
  GtDustBuffer *buf; /* Buffer for the lookahead corresponding to the
                        windowsize + linker. */
  bool buf_initialized;
  GtUword buf_readpos; /* where to read next value */
  GtUword buf_insertpos; /* where to write next value */
  GtUword buf_remaining;
  GtUword buf_size;

  unsigned char nuc_val1;
  unsigned char nuc_val2; /* last 2 nucleotide values */

  /* following variables named according to Morgulis et al. 2006 */
  unsigned int rv, rw;
  unsigned int cv[GT_MAXTRIPLETVALUE + 1],
               cw[GT_MAXTRIPLETVALUE + 1],
               ctmp[GT_MAXTRIPLETVALUE + 1];
  GtUword L_param;
  GtInl_Queue *w_queue;

  GtUword total_length,   /* Total chars written from file into buffer */
          current_length; /* Chars written in current sequence */
  GtUword current_pos_total; /* Total chars read out of Buffer */
  GtUword current_pos_local; /* Chars read from current sequence */

  GtUword last_seq_start; /* Start position of last sequence in buffer */

  GtUword mask_length; /* How many following chars read need to be masked */
  GtUword next_mask; /* Distance of the next masked region (for linking) */
  bool current_is_masked; /* Last char read was masked */

  bool masking_done; /* completed dust algorithm run */
  GtArrayGtDustRange masked_regions; /* Stores masked regions after 1st run */
  GtDustRange *current_region; /* Current masked region in the array */
  GtUword current_region_index; /* Index of current region in the array */

  /* Variables of the dust algorithm */
  bool echo;
  GtUword windowsize;
  GtUword linker;
  double threshold;
};

GtDustMasker* gt_dust_masker_new(bool echo, GtUword windowsize,
                                 double threshold,
                                 GtUword linker)
{
  GtDustMasker *dust_masker = gt_calloc(1, sizeof *dust_masker);
  dust_masker->echo = echo;
  dust_masker->windowsize = windowsize;
  dust_masker->threshold = threshold;
  dust_masker->linker = linker;

  dust_masker->buf_size = windowsize + linker;
  dust_masker->buf = gt_calloc(dust_masker->buf_size,
                               sizeof *(dust_masker->buf));
  dust_masker->buf_initialized = false;
  dust_masker->buf_readpos = 0;
  dust_masker->buf_insertpos = 0;
  dust_masker->buf_remaining = 0;

  dust_masker->nuc_val1 = 0;
  dust_masker->nuc_val2 = 0;

  dust_masker->rv = 0;
  dust_masker->rw = 0;
  dust_masker->L_param = 0;
  dust_masker->w_queue = gt_inl_queue_new(windowsize);

  dust_masker->total_length = 0;
  dust_masker->current_length = 0;
  dust_masker->current_pos_total = 0;
  dust_masker->current_pos_local = 0;

  dust_masker->last_seq_start = 0;

  GT_INITARRAY(&dust_masker->masked_regions,GtDustRange);
  dust_masker->mask_length = 0;
  dust_masker->next_mask = 0;
  dust_masker->current_is_masked = false;
  dust_masker->current_region_index = 0;

  dust_masker->masking_done = false;

  return dust_masker;
}

void gt_dust_masker_delete(GtDustMasker *dust_masker)
{
  if (dust_masker) {
    gt_free(dust_masker->buf);
    gt_inl_queue_delete(dust_masker->w_queue);
    GT_FREEARRAY(&dust_masker->masked_regions,GtDustRange);
    gt_free(dust_masker);
  }
}

/* Performs (val = val % limit) but assumes that val < 2*limit.
   Faster than a modulo operation. */
static inline void wrap_value_once(GtUword *val, GtUword limit)
{
  if (*val >= limit) {
    *val -= limit;
  }
}

static inline unsigned char nucleotide_value(char c)
{
  switch (c) {
    case ('a'):
    case ('A'): return 0;
    case ('c'):
    case ('C'): return 1;
    case ('g'):
    case ('G'): return 2;
    case ('t'):
    case ('T'): return 3;
#if GT_TREATNASRANDOM
    case ('n'):
    case ('N'): return rand() % 4;
#endif
    default: return 0;
  }
  return 0;
}

static inline void add_triplet_info(unsigned int *r, unsigned int *c,
                                    unsigned char t)
{
  *r = *r + c[t];
  c[t]++;
}

static inline void rem_triplet_info(unsigned int *r, unsigned int *c,
                                    unsigned char t)
{
  c[t]--;
  *r = *r - c[t];
}

static inline void find_perfect(GtDustMasker *dust_masker)
{
  unsigned int r = dust_masker->rv;
  GtUword linker_offset = 0;
  float new_score, max_score = 0.0, score_to_beat = 0.0;
  GtUword idx, best_idx, window_idx, buf_idx, length, step, readpos;
  bool found=false;

  memcpy(dust_masker->ctmp, dust_masker->cv,
         (GT_MAXTRIPLETVALUE + 1) * sizeof *(dust_masker->ctmp));

  if (dust_masker->current_length > dust_masker->windowsize) {
    linker_offset = GT_MIN(dust_masker->linker,
                           dust_masker->current_length
                             - dust_masker->windowsize);
  }

  readpos = dust_masker->buf_readpos;
  if (dust_masker->current_length < dust_masker->buf_size) {
    readpos = dust_masker->last_seq_start;
  }
  length = dust_masker->w_queue->noofelements - dust_masker->L_param - 1;
  for (step = 0; step <= length; step++) {
    unsigned char t;
    idx = length - step;
    window_idx = (readpos + idx + linker_offset);
    wrap_value_once(&window_idx, dust_masker->buf_size);
    score_to_beat = GT_MAX(score_to_beat,
                           dust_masker->buf[window_idx].max_score);
    t = gt_inl_queue_get_at_index(dust_masker->w_queue, idx);
    add_triplet_info(&r, dust_masker->ctmp, t);
    new_score = (float) r/(float)(dust_masker->w_queue->noofelements - idx - 1);
    if (new_score > dust_masker->threshold && new_score >= max_score &&
        new_score >= score_to_beat) {
      found = true;
      max_score = new_score;
      best_idx = idx;
      dust_masker->buf[window_idx].max_score = max_score;
    }
  }

  if (found) {
    buf_idx = (readpos + best_idx + linker_offset);
    wrap_value_once(&buf_idx, dust_masker->buf_size);
    dust_masker->buf[buf_idx].mask_length
      = GT_MAX(dust_masker->w_queue->noofelements + 2 - best_idx,
               dust_masker->buf[buf_idx].mask_length);

    if (dust_masker->linker > 1) {
      GtUword link_length = GT_MIN(best_idx + linker_offset,
                                   dust_masker->linker);
      link_length = GT_MIN(link_length, dust_masker->current_length-1);
      GtUword link_idx = readpos + best_idx + linker_offset - link_length;
      wrap_value_once(&link_idx, dust_masker->buf_size);
      dust_masker->buf[link_idx].next_mask
        = GT_MAX(dust_masker->buf[link_idx].next_mask,link_length);
    }
  }
}

static int gt_dust_masker_shift_window(GtDustMasker* dust_masker,
                                       GtSequenceBuffer* sb,
                                       GtError* err)
{
  int retval;
  GtUchar t_val=0;
  char t_orig=0;

  unsigned char triplet_val, s, nuc_val;

  /* read next character into buffer */
  retval = gt_sequence_buffer_next_with_original_raw(sb, &t_val, &t_orig, err);
  if (retval == -1) {
    return -1;
  } else {
    if (retval == 0) {
      return 0;
    }
  }
  dust_masker->buf_remaining++;
  dust_masker->current_length++;
  dust_masker->total_length++;
  dust_masker->buf[dust_masker->buf_insertpos].val = t_val;
  dust_masker->buf[dust_masker->buf_insertpos].orig = t_orig;
  dust_masker->buf[dust_masker->buf_insertpos].max_score = 0.0;
  dust_masker->buf[dust_masker->buf_insertpos].mask_length = 0;
  dust_masker->buf[dust_masker->buf_insertpos].next_mask = 0;
  dust_masker->buf_insertpos++;
  wrap_value_once(&dust_masker->buf_insertpos, dust_masker->buf_size);

  if (!dust_masker->masking_done) {
    if (t_val != GT_SEPARATOR) {
      nuc_val = nucleotide_value(t_orig);
      triplet_val = dust_masker->nuc_val1 * 16 +
                    dust_masker->nuc_val2 * 4 + nuc_val;
      dust_masker->nuc_val1 = dust_masker->nuc_val2;
      dust_masker->nuc_val2 = nuc_val;

      if (dust_masker->current_length <= 2) {
        return 1;
      }

      /* next part according to SHIFT_WINDOW-procedure in Morgulis et al 2006 */
      if (dust_masker->w_queue->noofelements >= dust_masker->windowsize - 2) {
        s = (char) gt_inl_queue_get(dust_masker->w_queue);
        rem_triplet_info(&dust_masker->rw, dust_masker->cw, s);
        if (dust_masker->L_param > dust_masker->w_queue->noofelements) {
          dust_masker->L_param--;
          rem_triplet_info(&dust_masker->rv, dust_masker->cv, s);
        }
      }
      gt_inl_queue_add(dust_masker->w_queue, triplet_val, false);
      dust_masker->L_param++;
      add_triplet_info(&dust_masker->rw, dust_masker->cw, triplet_val);
      add_triplet_info(&dust_masker->rv, dust_masker->cv, triplet_val);
      if (dust_masker->cv[triplet_val] > (2 * dust_masker->threshold)) {
        do {
          s = gt_inl_queue_get_at_index(dust_masker->w_queue,
                                        dust_masker->w_queue->noofelements
                                          - dust_masker->L_param);
          rem_triplet_info(&dust_masker->rv, dust_masker->cv, s);
          dust_masker->L_param--;
        } while (s != triplet_val);
      }
      if (dust_masker->rw > (float) dust_masker->L_param
          * dust_masker->threshold) {
        find_perfect(dust_masker);
      }
    } else {
      /* Reset variables for next sequence in multiseq fasta file. */
      dust_masker->last_seq_start = dust_masker->buf_insertpos;
      dust_masker->nuc_val1 = 0;
      dust_masker->nuc_val2 = 0;
      dust_masker->rv = 0;
      dust_masker->rw = 0;
      dust_masker->L_param = 0;
      dust_masker->current_length = 0;
      memset(dust_masker->cv, 0,
             sizeof *(dust_masker->cv) * (GT_MAXTRIPLETVALUE + 1));
      memset(dust_masker->cw, 0,
             sizeof *(dust_masker->cw) * (GT_MAXTRIPLETVALUE + 1));
      gt_inl_queue_delete(dust_masker->w_queue);
      dust_masker->w_queue = gt_inl_queue_new(dust_masker->windowsize);
    }
  }
  return 1;
}

int gt_dust_masker_next_with_original(GtDustMasker *dust_masker,
                                      GtSequenceBuffer* sb,
                                      GtUchar *val, char *orig, GtError* err)
{
  int retval;

  if (dust_masker->masking_done) {
    retval = gt_sequence_buffer_next_with_original_raw(sb, val, orig, err);
    if (retval == -1) {
      return -1;
    }
    if (retval == 0) {
      if (dust_masker->echo) {
        if (dust_masker->current_pos_local % 60 != 0) {
          printf("\n");
        }
      }

      dust_masker->current_region_index = 0;
      dust_masker->current_pos_total = 0;
      dust_masker->current_pos_local = 0;
      return 0;
    }
    if (dust_masker->current_region_index <
        dust_masker->masked_regions.nextfreeGtDustRange) {
      GtDustRange range;
      range = dust_masker->masked_regions.spaceGtDustRange[
                         dust_masker->current_region_index];
      if (dust_masker->current_pos_total >= range.begin) {
        if (dust_masker->current_pos_total <= range.end &&
            *val != GT_SEPARATOR) {
          if (*orig >= 'A' && *orig <= 'Z') {
            *orig += 32;
          }
          *val = GT_WILDCARD;
        } else {
          dust_masker->current_region_index++;
        }
      }
    }
    dust_masker->current_pos_total++;
    dust_masker->current_pos_local++;
    if (dust_masker->echo) {
      if (*val == GT_SEPARATOR) {
        if (dust_masker->current_pos_local % 60 != 1) {
          printf("\n");
        }
        dust_masker->current_pos_local = 0;
      } else {
        printf("%c",*orig);
        if (dust_masker->current_pos_local % 60 == 0) {
          printf("\n");
        }
      }
    }
  } else {
    if (!dust_masker->buf_initialized) {
      GtUword idx;
      for (idx = 0; idx < dust_masker->buf_size; idx++) {
        retval = gt_dust_masker_shift_window(dust_masker, sb, err);
        if (retval == -1) {
          return -1;
        } else {
          if (retval == 0) {
            break;
          }
        }
      }
      dust_masker->buf_initialized = true;
    }
    if (dust_masker->buf_remaining > 0) {
      dust_masker->buf_remaining--;
      dust_masker->mask_length
        = GT_MAX(dust_masker->mask_length,
                 dust_masker->buf[dust_masker->buf_readpos].mask_length);

      if (dust_masker->linker > 1) {
        dust_masker->next_mask
          = GT_MAX(dust_masker->next_mask,
                   dust_masker->buf[dust_masker->buf_readpos].next_mask);
        if (dust_masker->mask_length > 0) {
          dust_masker->mask_length = GT_MAX(dust_masker->mask_length,
                                            dust_masker->next_mask);
        }
        if (dust_masker->next_mask > 0) {
          dust_masker->next_mask--;
        }
      }

      *val = dust_masker->buf[dust_masker->buf_readpos].val;
      *orig = dust_masker->buf[dust_masker->buf_readpos].orig;
      if (dust_masker->mask_length > 0 && *val != GT_SEPARATOR) {
        if (*orig >= 'A' && *orig <= 'Z') {
          *orig += 32;
        }
        *val = GT_WILDCARD;
      }

      /* Insert into masked_regions */
      if (dust_masker->mask_length > 0) {
        if (dust_masker->current_is_masked == false) {
          GT_GETNEXTFREEINARRAY(dust_masker->current_region,
                                &dust_masker->masked_regions,
                                GtDustRange, GT_ARRAYEXTENDSIZE);
          dust_masker->current_region->begin = dust_masker->current_pos_total;
        }
        dust_masker->current_region->end = dust_masker->current_pos_total;
        dust_masker->current_is_masked = true;
        dust_masker->mask_length--;
      } else {
        dust_masker->current_is_masked = false;
      }
      dust_masker->buf_readpos++;
      wrap_value_once(&dust_masker->buf_readpos, dust_masker->buf_size);
    } else {
      dust_masker->masking_done = true;
      dust_masker->current_pos_total = 0;
      return 0;
    }
    dust_masker->current_pos_total++;
    retval = gt_dust_masker_shift_window(dust_masker, sb, err);
    if (retval == -1) {
      return -1;
    }
  }
  return 1;
}
