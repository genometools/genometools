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

#define MAXTRIPLETVALUE 63
#define ARRAYEXTENDSIZE 100
#define TREATNASRANDOM 0

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
  unsigned int cv[MAXTRIPLETVALUE + 1],
               cw[MAXTRIPLETVALUE + 1],
               ctmp[MAXTRIPLETVALUE + 1];
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

GtDustMasker* gt_dustmasker_new(bool echo, GtUword windowsize, double threshold,
                                GtUword linker)
{
  GtDustMasker *dm = gt_calloc(1, sizeof *dm);
  dm->echo = echo;
  dm->windowsize = windowsize;
  dm->threshold = threshold;
  dm->linker = linker;

  dm->buf_size = windowsize + linker;
  dm->buf = gt_calloc(dm->buf_size, sizeof *(dm->buf));
  dm->buf_initialized = false;
  dm->buf_readpos = 0;
  dm->buf_insertpos = 0;
  dm->buf_remaining = 0;

  dm->nuc_val1 = 0;
  dm->nuc_val2 = 0;

  dm->rv = 0;
  dm->rw = 0;
  dm->L_param = 0;
  dm->w_queue = gt_inl_queue_new(windowsize);

  dm->total_length = 0;
  dm->current_length = 0;
  dm->current_pos_total = 0;
  dm->current_pos_local = 0;

  dm->last_seq_start = 0;

  GT_INITARRAY(&dm->masked_regions,GtDustRange);
  dm->mask_length = 0;
  dm->next_mask = 0;
  dm->current_is_masked = false;
  dm->current_region_index = 0;

  dm->masking_done = false;

  return dm;
}

void gt_dustmasker_delete(GtDustMasker *dm)
{
  if (dm)
  {
    gt_free(dm->buf);
    gt_inl_queue_delete(dm->w_queue);
    GT_FREEARRAY(&dm->masked_regions,GtDustRange);
    gt_free(dm);
  }
}

/* Performs (val = val % limit) but assumes that val < 2*limit.
   Faster than a modulo operation. */
static inline void wrap_value_once(GtUword *val, GtUword limit)
{
  if (*val >= limit)
    *val -= limit;
}

static inline unsigned char nucleotide_value(char c)
{
  switch (c)
  {
    case ('a'):
    case ('A'): return 0;
    case ('c'):
    case ('C'): return 1;
    case ('g'):
    case ('G'): return 2;
    case ('t'):
    case ('T'): return 3;
#if TREATNASRANDOM
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

static inline void find_perfect(GtDustMasker *dm)
{
  unsigned int r = dm->rv;
  GtUword linker_offset = 0;
  float new_score, max_score = 0.0, score_to_beat = 0.0;
  GtUword idx, best_idx, window_idx, buf_idx, length, step, readpos;
  bool found=false;

  memcpy(dm->ctmp, dm->cv, (MAXTRIPLETVALUE + 1) * sizeof *(dm->ctmp));

  if (dm->current_length > dm->windowsize)
    linker_offset = GT_MIN(dm->linker, dm->current_length - dm->windowsize);

  readpos = dm->buf_readpos;
  if (dm->current_length < dm->buf_size)
    readpos = dm->last_seq_start;

  length = dm->w_queue->noofelements - dm->L_param - 1;
  for (step = 0; step <= length; step++)
  {
    idx = length - step;
    window_idx = (readpos + idx + linker_offset);
    wrap_value_once(&window_idx, dm->buf_size);
    score_to_beat = GT_MAX(score_to_beat, dm->buf[window_idx].max_score);
    unsigned char t = gt_inl_queue_get_at_index(dm->w_queue, idx);
    add_triplet_info(&r, dm->ctmp, t);
    new_score = (float)r / (float)(dm->w_queue->noofelements - idx - 1);
    if (new_score > dm->threshold && new_score >= max_score &&
        new_score >= score_to_beat)
    {
      found = true;
      max_score = new_score;
      best_idx = idx;
      dm->buf[window_idx].max_score = max_score;
    }
  }

  if (found)
  {
    buf_idx = (readpos + best_idx + linker_offset);
    wrap_value_once(&buf_idx, dm->buf_size);
    dm->buf[buf_idx].mask_length
      = GT_MAX(dm->w_queue->noofelements + 2 - best_idx,
               dm->buf[buf_idx].mask_length);

    if (dm->linker > 1)
    {
      GtUword link_length = GT_MIN(best_idx + linker_offset, dm->linker);
      link_length = GT_MIN(link_length, dm->current_length-1);
      GtUword link_idx = readpos + best_idx + linker_offset - link_length;
      wrap_value_once(&link_idx, dm->buf_size);
      dm->buf[link_idx].next_mask = GT_MAX(dm->buf[link_idx].next_mask,
                                           link_length);
    }
  }
}

static int gt_dustmasker_shift_window(GtDustMasker* dm, GtSequenceBuffer* sb,
                                      GtError* err)
{
  int retval;
  GtUchar t_val=0;
  char t_orig=0;

  unsigned char triplet_val, s, nuc_val;

  /* read next character into buffer */
  retval = gt_sequence_buffer_next_with_original_raw(sb, &t_val, &t_orig, err);
  if (retval == -1)
    return -1;
  else if (retval == 0)
    return 0;

  dm->buf_remaining++;
  dm->current_length++;
  dm->total_length++;
  dm->buf[dm->buf_insertpos].val = t_val;
  dm->buf[dm->buf_insertpos].orig = t_orig;
  dm->buf[dm->buf_insertpos].max_score = 0.0;
  dm->buf[dm->buf_insertpos].mask_length = 0;
  dm->buf[dm->buf_insertpos].next_mask = 0;
  dm->buf_insertpos++;
  wrap_value_once(&dm->buf_insertpos, dm->buf_size);

  if (!dm->masking_done)
  {
    if (t_val != GT_SEPARATOR)
    {
      nuc_val = nucleotide_value(t_orig);
      triplet_val = dm->nuc_val1 * 16 + dm->nuc_val2 * 4 + nuc_val;
      dm->nuc_val1 = dm->nuc_val2;
      dm->nuc_val2 = nuc_val;

      if (dm->current_length <= 2)
        return 1;

      /* next part according to SHIFT_WINDOW-procedure in Morgulis et al 2006 */
      if (dm->w_queue->noofelements >= dm->windowsize - 2)
      {
        s = (char)gt_inl_queue_get(dm->w_queue);
        rem_triplet_info(&dm->rw, dm->cw, s);
        if (dm->L_param > dm->w_queue->noofelements)
        {
          dm->L_param--;
          rem_triplet_info(&dm->rv, dm->cv, s);
        }
      }
      gt_inl_queue_add(dm->w_queue, triplet_val, false);
      dm->L_param++;
      add_triplet_info(&dm->rw, dm->cw, triplet_val);
      add_triplet_info(&dm->rv, dm->cv, triplet_val);
      if (dm->cv[triplet_val] > (2 * dm->threshold))
      {
        do {
          s = gt_inl_queue_get_at_index(dm->w_queue,
                                      dm->w_queue->noofelements - dm->L_param);
          rem_triplet_info(&dm->rv, dm->cv, s);
          dm->L_param--;
        } while (s != triplet_val);
      }
      if (dm->rw > (float)dm->L_param * dm->threshold)
        find_perfect(dm);
    }
    else {
      /* Reset variables for next sequence in multiseq fasta file. */
      dm->last_seq_start = dm->buf_insertpos;
      dm->nuc_val1 = 0;
      dm->nuc_val2 = 0;
      dm->rv = 0;
      dm->rw = 0;
      dm->L_param = 0;
      dm->current_length = 0;
      memset(dm->cv, 0, sizeof *(dm->cv) * (MAXTRIPLETVALUE + 1));
      memset(dm->cw, 0, sizeof *(dm->cw) * (MAXTRIPLETVALUE + 1));
      gt_inl_queue_delete(dm->w_queue);
      dm->w_queue = gt_inl_queue_new(dm->windowsize);
    }
  }
  return 1;
}

int gt_dustmasker_next_with_original(GtDustMasker *dm, GtSequenceBuffer* sb,
                                     GtUchar *val, char *orig, GtError* err)
{
  int retval;
  if (dm->masking_done)
  {
    retval = gt_sequence_buffer_next_with_original_raw(sb, val, orig, err);
    if (retval == -1)
    {
      return -1;
    } else if (retval == 0)
    {
      if (dm->echo)
        if (dm->current_pos_local % 60 != 0)
          printf("\n");

      dm->current_region_index = 0;
      dm->current_pos_total = 0;
      dm->current_pos_local = 0;
      return 0;
    }
    else if (dm->current_region_index < dm->masked_regions.nextfreeGtDustRange)
    {
      GtDustRange range;
      range = dm->masked_regions.spaceGtDustRange[dm->current_region_index];
      if (dm->current_pos_total >= range.begin)
      {
        if (dm->current_pos_total <= range.end && *val != GT_SEPARATOR)
        {
          if (*orig >= 'A' && *orig <= 'Z')
            *orig += 32;
          *val = GT_WILDCARD;
        } else
        {
          dm->current_region_index++;
        }
      }
    }
    dm->current_pos_total++;
    dm->current_pos_local++;
    if (dm->echo)
    {
      if (*val == GT_SEPARATOR)
      {
        if (dm->current_pos_local % 60 != 1)
          printf("\n");
        dm->current_pos_local = 0;
      } else
      {
        printf("%c",*orig);
        if (dm->current_pos_local % 60 == 0)
          printf("\n");
      }
    }
  } else
  {
    if (!dm->buf_initialized)
    {
      GtUword idx;
      for (idx = 0; idx < dm->buf_size; idx++)
      {
        retval = gt_dustmasker_shift_window(dm, sb, err);
        if (retval == -1)
          return -1;
        else if (retval == 0)
          break;
      }
      dm->buf_initialized = true;
    }
    if (dm->buf_remaining > 0)
    {
      dm->buf_remaining--;
      dm->mask_length = GT_MAX(dm->mask_length,
                               dm->buf[dm->buf_readpos].mask_length);

      if (dm->linker > 1)
      {
        dm->next_mask = GT_MAX(dm->next_mask,
                               dm->buf[dm->buf_readpos].next_mask);
        if (dm->mask_length > 0)
          dm->mask_length = GT_MAX(dm->mask_length, dm->next_mask);
        if (dm->next_mask > 0)
          dm->next_mask--;
      }

      *val = dm->buf[dm->buf_readpos].val;
      *orig = dm->buf[dm->buf_readpos].orig;
      if (dm->mask_length > 0 && *val != GT_SEPARATOR)
      {
        if (*orig >= 'A' && *orig <= 'Z')
          *orig += 32;
        *val = GT_WILDCARD;
      }

      /* Insert into masked_regions */
      if (dm->mask_length > 0)
      {
        if (dm->current_is_masked == false)
        {
          GT_GETNEXTFREEINARRAY(dm->current_region, &dm->masked_regions,
                                GtDustRange, ARRAYEXTENDSIZE);
          dm->current_region->begin = dm->current_pos_total;
        }
        dm->current_region->end = dm->current_pos_total;
        dm->current_is_masked = true;
        dm->mask_length--;
      } else
      {
        dm->current_is_masked = false;
      }
      dm->buf_readpos++;
      wrap_value_once(&dm->buf_readpos, dm->buf_size);
    } else
    {
      dm->masking_done = true;
      dm->current_pos_total = 0;
      return 0;
    }
    dm->current_pos_total++;
    retval = gt_dustmasker_shift_window(dm, sb, err);
    if (retval == -1)
    {
      return -1;
    }
  }
  return 1;
}
