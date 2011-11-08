/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "extended/mutate.h"

#define MUTATED_DESC_PRIMER " [mutated with rate "

static char* mutate_description(const char *description, unsigned int rate)
{
  unsigned long mutated_description_len;
  char *mutated_description;
  GT_UNUSED int rval;
  gt_assert(description);
  gt_assert(rate <= 100);
  mutated_description_len = strlen(description) + strlen(MUTATED_DESC_PRIMER)
                            + 3  /* for the rate */
                            + 1  /* terminal ']' */
                            + 1; /* terminal '\n' */
  mutated_description = gt_malloc(sizeof (char) * mutated_description_len);
  rval = snprintf(mutated_description, mutated_description_len, "%s%s%u]",
                  description, MUTATED_DESC_PRIMER, rate);
  gt_assert(rval < mutated_description_len);
  return mutated_description;
}

static char random_character(GtAlphabet *alphabet, bool upper_case)
{
  /* we do not want to get wildcard characters */
  char random_char = gt_alphabet_decode(alphabet,
                           gt_rand_max(gt_alphabet_num_of_chars(alphabet) - 1));
  if (upper_case)
    return toupper(random_char);
  return tolower(random_char);
}

static char* mutate_seq(const char *seq, unsigned long len,
                        GtAlphabet *alphabet, unsigned int rate)
{
  unsigned long i, j, allocated, substitution_events = 0, insertion_events = 0,
                deletion_events = 0, total_events = 0;
  double rand_prob, mutate_prob;
  char *mutated_seq;
  bool was_upper;
  gt_assert(seq && alphabet);
  gt_assert(rate <= 100);
  mutate_prob = (double) rate / 100.0;
  allocated = len * 2; /* XXX: possibly reduce this memory consumption */
  mutated_seq = gt_malloc(sizeof (char) * allocated);
  for (i = 0, j = 0; i < len; i++) {
    if (isupper(seq[i]))
      was_upper = true;
    else
      was_upper = false;
    if (gt_rand_0_to_1() <= mutate_prob) {
      /* mutate */
      rand_prob = gt_rand_0_to_1();
      if (rand_prob <= 0.8) {
        /* substitution (80% probability) */
        mutated_seq[j++] = random_character(alphabet, was_upper);
        substitution_events++;
      }
      else if (rand_prob <= 0.9) {
        /* insertion (10% probability) */
        mutated_seq[j++] = seq[i]; /* keep orig. character */
        mutated_seq[j++] = random_character(alphabet, was_upper);
        insertion_events++;
      }
      else {
      /* deletion (10% probability) */
        deletion_events++;
      }
      total_events++;
    }
    else
      mutated_seq[j++] = seq[i]; /* keep original character */
  }
  mutated_seq[j] = '\0'; /* terminate */
  gt_log_log("total number of mutation events: %lu", total_events);
  gt_log_log("number of substitution events: %lu", substitution_events);
  gt_log_log("number of insertion events: %lu", insertion_events);
  gt_log_log("number of deletion events: %lu", deletion_events);
  return mutated_seq;
}

GtSeq* gt_mutate_seq(const char *description, const char *orig_seq,
                     unsigned long len, GtAlphabet *alphabet, unsigned int rate)
{
  char *mutated_description, *mutated_seq;
  GtSeq *seq;
  gt_assert(description && orig_seq && alphabet);
  gt_assert(rate <= 100);
  mutated_description = mutate_description(description, rate);
  mutated_seq = mutate_seq(orig_seq, len, alphabet, rate);
  seq = gt_seq_new_own(mutated_seq, strlen(mutated_seq), alphabet);
  gt_seq_set_description_own(seq, mutated_description);
  return seq;
}
