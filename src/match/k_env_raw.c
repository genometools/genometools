/*
  Copyright (c) 2017 Julian Elvers <julian.elvers@gmail.com>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include <string.h>
#include "core/array2dim_api.h"
#include "core/ma_api.h"
#include "core/score_matrix.h"
#include "match/turnwheels.h"
#include "match/initbasepower.h"
#include "k_env_raw.h"

/* Building the k_environment for a given sequence and contains
   helper methods and structures to do so */

/* Macro to define the size an array of <KenvRawElement> will be enlarged */
#define GT_KENV_RAW_REALLOC_SIZE(SIZE) (256 + (SIZE) * 0.2)

/***** GtKenvScoreHandling *****/
typedef struct
{
  int score;
  unsigned int charCode;
} GtKenvScoreCharCodePair;

typedef struct
{
  GtKenvScoreCharCodePair **score_pairs;
  int *max_scores;
} GtKenvScoreStructure;

static GtKenvScoreStructure *gt_kenv_score_structure_new(unsigned dim)
{
  GtKenvScoreStructure *scores = gt_malloc(1 * sizeof (*scores));
  gt_array2dim_calloc(scores->score_pairs, dim, dim);
  scores->max_scores = gt_calloc(dim, sizeof (*scores->max_scores));

  return scores;
}

static void gt_kenv_score_structure_delete(GtKenvScoreStructure *scores)
{
  if (scores != NULL)
  {
    if (scores->score_pairs != NULL)
      gt_array2dim_delete(scores->score_pairs);
    if (scores->max_scores != NULL)
      gt_free(scores->max_scores);
    gt_free(scores);
  }
}

static void gt_kenv_preprocess_scores(const GtKenvGenerator *kenv_gen);

/***** GtKenvGenerator *****/
struct GtKenvGenerator{
  unsigned int q;
  unsigned int th;
  unsigned int next_seqnum;
  unsigned int alph_size;
  GtCodetype *alph_powers;
  const GtUchar* alph_symbolmap;
  GtKenvScoreStructure *scores;
  char *buffer_Char;
  GtUchar *buffer_GtUchar;
  int *buffer_Int;
  bool allow_x;
  GtKenvRawElement *spaceGtKenvRawElement;
  GtUword allocatedGtKenvRawElement;
  GtUword nextfreeGtKenvRawElement;
  GtUword *kmer_offsets;
  GtUword max_code_length;
  bool preprocess_env;
  GtScoreMatrix* score_matrix;
};

/* Header declaration so it can be used from the constructor */
static void gt_kenv_generator_build_kmer_index(GtKenvGenerator *kenv_gen);

static void gt_kenv_gen_finalize_offsets(GtKenvGenerator *kenv_gen,
                                         GtUword offset_size)
{
  GtUword idx, curr_size, next_size;
  gt_assert(kenv_gen);
  curr_size = kenv_gen->kmer_offsets[0];
  kenv_gen->max_code_length = kenv_gen->kmer_offsets[offset_size - 1];
  for (idx = 0; idx < offset_size; idx++)
  {
    kenv_gen->kmer_offsets[idx] -= curr_size;
    if (idx < (offset_size - 1))
    {
      next_size = kenv_gen->kmer_offsets[idx + 1];
      kenv_gen->kmer_offsets[idx + 1] += curr_size;
      curr_size = next_size;
    }
    if (idx > 0)
    {
      kenv_gen->kmer_offsets[idx] += kenv_gen->kmer_offsets[idx - 1];
    }
  }
}

GtKenvGenerator *gt_kenv_generator_new(const GtKenvAlphabet *kenv_alphabet,
                                       unsigned int q_val, unsigned int th_val,
                                       bool allow_x_val, bool preprocess,
                                       GtError *err)
{
  const char* blosum_file = "testdata/BLOSUM62";
  const char* dna_file = "DNA_SCORES";
  GtKenvGenerator *kenv_gen = gt_malloc(sizeof (*kenv_gen));
  bool use_dna = gt_kenv_alphabet_is_dna(kenv_alphabet);
  GtUword offset_size;

  kenv_gen->q = q_val;
  kenv_gen->th = th_val;
  kenv_gen->next_seqnum = 0;
  kenv_gen->allow_x = allow_x_val;
  kenv_gen->alph_size = gt_kenv_alphabet_get_alph_size(kenv_alphabet);
  kenv_gen->alph_powers = gt_initbasepower(kenv_gen->alph_size, q_val);
  kenv_gen->alph_symbolmap = gt_kenv_alphabet_get_characters(kenv_alphabet);

  kenv_gen->score_matrix
    = gt_score_matrix_new_read(use_dna ? dna_file
                                       : blosum_file,
                               gt_kenv_alphabet_get_alphabet(kenv_alphabet),
                               err);
  if (kenv_gen->score_matrix == NULL)
  {
    gt_free(kenv_gen);
    return NULL;
  }
  kenv_gen->scores = gt_kenv_score_structure_new(kenv_gen->alph_size);
  gt_kenv_preprocess_scores(kenv_gen);
  kenv_gen->buffer_Char = gt_malloc(sizeof (*kenv_gen->buffer_Char) *
                                    (kenv_gen->q + 1));
  kenv_gen->buffer_GtUchar = gt_calloc(kenv_gen->q + 1,
                                        sizeof (*kenv_gen->buffer_GtUchar));
  kenv_gen->buffer_Int = gt_calloc((kenv_gen->q + 1),
                                  sizeof (*kenv_gen->buffer_Int));

  kenv_gen->spaceGtKenvRawElement = NULL;
  kenv_gen->allocatedGtKenvRawElement = 0;
  kenv_gen->nextfreeGtKenvRawElement = 0;
  offset_size = kenv_gen->alph_powers[q_val];
  kenv_gen->kmer_offsets = NULL;
  kenv_gen->preprocess_env = preprocess;
  kenv_gen->max_code_length = 0;

  if (preprocess)
  {
    /* TODO: check, whether we need to allocate this always?? */
    kenv_gen->kmer_offsets = gt_calloc(offset_size,
                                      sizeof (*kenv_gen->kmer_offsets));
    gt_kenv_generator_build_kmer_index(kenv_gen);

    gt_kenv_gen_finalize_offsets(kenv_gen, offset_size);
  }
  return kenv_gen;
}

void gt_kenv_generator_delete(GtKenvGenerator *kenv_gen)
{
  if (kenv_gen != NULL)
  {
    if (kenv_gen->alph_powers != NULL)
      gt_free(kenv_gen->alph_powers);
    if (kenv_gen->scores != NULL)
      gt_kenv_score_structure_delete(kenv_gen->scores);
    if (kenv_gen->buffer_Char != NULL)
      gt_free(kenv_gen->buffer_Char);
    if (kenv_gen->buffer_GtUchar != NULL)
      gt_free(kenv_gen->buffer_GtUchar);
    if (kenv_gen->buffer_Int != NULL)
      gt_free(kenv_gen->buffer_Int);
    if (kenv_gen->spaceGtKenvRawElement != NULL)
      gt_free(kenv_gen->spaceGtKenvRawElement);
    if (kenv_gen->kmer_offsets != NULL)
      gt_free(kenv_gen->kmer_offsets);
    gt_score_matrix_delete(kenv_gen->score_matrix);
    gt_free(kenv_gen);
  }
}

int gt_kenv_generator_get_max_score(const GtKenvGenerator *kenv_gen)
{
  unsigned idx;
  int max_score = INT_MIN, curr_score;

  gt_assert(kenv_gen && kenv_gen->scores);

  for (idx = 0; idx < kenv_gen->alph_size; idx++)
  {
    curr_score = kenv_gen->scores->max_scores[idx];
    if (curr_score > max_score)
      max_score = curr_score;
  }
  return max_score;
}

static int kenv_compare_charcodepairs(const void *p1, const void *p2)
{
  const GtKenvScoreCharCodePair *elem1 = p1, *elem2 = p2;

  if (elem1->score < elem2->score)
      return 1;
  else if (elem1->score > elem2->score)
      return -1;
  else
      return 0;
}

/* Preprocesses the scores from the given score_matrix into the structure
   inside the GtKenvGenerator */
static void gt_kenv_preprocess_scores(const GtKenvGenerator *kenv_gen)
{
  unsigned int row, col, score_dim = kenv_gen->alph_size;

  gt_assert(kenv_gen->scores);

  for (row = 0; row < score_dim; row++)
  {
    for (col = 0; col < score_dim; col++)
    {
      kenv_gen->scores->score_pairs[row][col].score
        = gt_score_matrix_get_score(kenv_gen->score_matrix, row, col);
      kenv_gen->scores->score_pairs[row][col].charCode = col;
    }
    /* Sort each column by descending scores */
    qsort(kenv_gen->scores->score_pairs[row], score_dim,
          sizeof (**(kenv_gen->scores->score_pairs)),
          kenv_compare_charcodepairs);
    kenv_gen->scores->max_scores[row] =
                          kenv_gen->scores->score_pairs[row][0].score;
  }
}

GtCodetype gt_kenv_get_kmer_code(const GtKenvGenerator *kenv_gen,
                                 const GtUchar *word)
{
  unsigned int idx, weight, q_val;
  GtCodetype kmer_code = 0;

  q_val = gt_kenv_generator_get_q(kenv_gen);
  for (idx = 0; idx < q_val; idx++)
  {
    weight = q_val - 1 - idx;
    kmer_code += kenv_gen->alph_powers[weight] * word[idx];
  }

  return kmer_code;
}

/* Will fill the int buffer inside the kenv gen */
static void gt_kenv_generator_calc_pos_thresholds(GtKenvGenerator *kenv_gen)
{
  unsigned int kmer_pos;

  gt_assert(kenv_gen && kenv_gen->buffer_Int && kenv_gen->buffer_GtUchar);
  kenv_gen->buffer_Int[kenv_gen->q - 1] = kenv_gen->th;
  for (kmer_pos = kenv_gen->q - 1; kmer_pos > 0; kmer_pos--)
  {
    /* Define threshold using the maximal score for this character */
    kenv_gen->buffer_Int[kmer_pos - 1] = kenv_gen->buffer_Int[kmer_pos] -
        kenv_gen->scores->max_scores[kenv_gen->buffer_GtUchar[kmer_pos]];
  }
}

/* Fills both the kmer buffer and the threshold buffer inside kenv_gen */
void gt_kenv_generator_fill_buffers_from_sequence(GtKenvGenerator *kenv_gen,
                                                  const GtUchar *kmer)
{
  unsigned int idx;

  gt_assert(kenv_gen && kenv_gen->buffer_Int && kenv_gen->buffer_GtUchar);
  for (idx = 0; idx < kenv_gen->q; idx++)
  {
    kenv_gen->buffer_GtUchar[idx] = kmer[idx];
  }
  gt_kenv_generator_calc_pos_thresholds(kenv_gen);
}

unsigned int gt_kenv_generator_get_q(const GtKenvGenerator *kenv_gen)
{
  gt_assert(kenv_gen);

  return kenv_gen->q;
}

unsigned int gt_kenv_generator_get_th(const GtKenvGenerator *kenv_gen)
{
  gt_assert(kenv_gen);

  return kenv_gen->th;
}

bool gt_kenv_generator_get_allow_x(const GtKenvGenerator *kenv_gen)
{
  gt_assert(kenv_gen);

  return kenv_gen->allow_x;
}

const int *gt_kenv_generator_get_max_score_each_char(
                                                const GtKenvGenerator *kenv_gen)
{
  gt_assert(kenv_gen && kenv_gen->scores->max_scores);

  return kenv_gen->scores->max_scores;
}

/* Fills the array inside the kenv_generator with the GtKenvRawElements */
static void gt_kenv_raw_element_build(GtKenvGenerator *kenv_gen,
                                      GtCodetype kmer, unsigned int curr_len,
                                      int kmer_score, GtCodetype refcode)
{
  /* Check for full word, then we found a valid element for the environment */
  if (curr_len == kenv_gen->q)
  {
    GtUword new_size, next_free;
    /* Reallocate memory if necessary */
    next_free = kenv_gen->nextfreeGtKenvRawElement;
    if (next_free == kenv_gen->allocatedGtKenvRawElement)
    {
      new_size = kenv_gen->allocatedGtKenvRawElement +
                GT_KENV_RAW_REALLOC_SIZE(kenv_gen->allocatedGtKenvRawElement);
      kenv_gen->spaceGtKenvRawElement =
                  gt_realloc(kenv_gen->spaceGtKenvRawElement,
                      new_size * sizeof (*(kenv_gen->spaceGtKenvRawElement)));
      kenv_gen->allocatedGtKenvRawElement = new_size;
    }
    kenv_gen->spaceGtKenvRawElement[next_free].kmer_neighbor = kmer;
    kenv_gen->spaceGtKenvRawElement[next_free].score = kmer_score;
    kenv_gen->nextfreeGtKenvRawElement++;
    if (kenv_gen->kmer_offsets != NULL)
      kenv_gen->kmer_offsets[refcode]++;
  }
  /* Else try to append valid chars from the sorted score matrix */
  else
  {
    unsigned int idx;
    GtUchar curr_code = kenv_gen->buffer_GtUchar[curr_len];
    GtKenvScoreCharCodePair *curr_col_ptr;

    gt_assert(curr_code < kenv_gen->alph_size);
    curr_col_ptr = kenv_gen->scores->score_pairs[curr_code];
    for (idx = 0; idx < kenv_gen->alph_size; idx++)
    {
      int new_score = curr_col_ptr[idx].score;
      if ((kmer_score + new_score) >= kenv_gen->buffer_Int[curr_len])
      {
        GtCodetype new_kmer_code = kmer + curr_col_ptr[idx].charCode *
            kenv_gen->alph_powers[kenv_gen->q - 1 - curr_len];

        gt_kenv_raw_element_build(kenv_gen, new_kmer_code, curr_len + 1,
                                  kmer_score + new_score, refcode);
      }
      else
      {
        /* Since the scores are sorted, this kmer prefix won't lead to
           the score threshold as determined by the pos_th. */
        break;
      }
    }
  }
}

static void gt_kenv_raw_kmer_from_code(GtKenvGenerator *kenv_gen, GtUword code)
{
  unsigned int idx;

  gt_assert(kenv_gen->buffer_Char);
  for (idx = 0; idx < kenv_gen->q; idx++)
  {
    unsigned int temp = kenv_gen->alph_powers[kenv_gen->q - 1 - idx],
                 curr_code = code / temp;
    code = code - curr_code * temp;
    gt_assert(curr_code < kenv_gen->alph_size);
    kenv_gen->buffer_Char[idx] = kenv_gen->alph_symbolmap[curr_code];
    kenv_gen->buffer_GtUchar[idx] = curr_code;
  }
  kenv_gen->buffer_Char[kenv_gen->q] = '\0';
}

static void gt_kenv_generator_build_kmer_index(GtKenvGenerator *kenv_gen)
{
  Turningwheel *tw;
  const unsigned int *curr_kmer;
  unsigned int idx, minchanged;
  GtCodetype kmer_code = 0;

  gt_assert(kenv_gen && kenv_gen->buffer_GtUchar);
  tw = gt_turningwheel_new(kenv_gen->q, kenv_gen->alph_size);

  do
  {
    /* This gives a reference to the current kmer. We need to determine
       pos_thresholds for this and then give it to the recurrent function
       to build the respective k_environment. */
    curr_kmer = gt_turningwheel_state(tw);
    minchanged = gt_turningwheel_minchanged(tw);
    for (idx = minchanged; idx < kenv_gen->q; idx++)
    {
      kenv_gen->buffer_GtUchar[idx] = curr_kmer[idx];
    }
    gt_kenv_generator_calc_pos_thresholds(kenv_gen);
    gt_kenv_raw_element_build(kenv_gen, 0, 0, 0, kmer_code);

    kmer_code++;  /* Used to reference the kmer in the index later on */
  } while (gt_turningwheel_next(tw));

  gt_turningwheel_delete(tw);
}

const GtKenvRawElement *gt_kenv_raw_element_list(GtUword *length,
                                GtCodetype refcode,
                                GtKenvGenerator *kenv_gen)
{
  if (!kenv_gen->preprocess_env)
  {
    /* TODO: hard workaround, not really suitable this way */
    kenv_gen->nextfreeGtKenvRawElement = 0;
    /* Write the kmer into the buffer in kenv_gen */
    gt_kenv_raw_kmer_from_code(kenv_gen, refcode);
    /* Generate the positional thresholds for the score based on the kmer */
    gt_kenv_generator_calc_pos_thresholds(kenv_gen);
    gt_kenv_raw_element_build(kenv_gen, 0, 0, 0, refcode);
    *length = kenv_gen->nextfreeGtKenvRawElement;
    return kenv_gen->spaceGtKenvRawElement;
  }
  else
  {
    if (refcode == (kenv_gen->alph_powers[kenv_gen->q] - 1))
    {
      *length = kenv_gen->max_code_length;
    }
    else
    {
      *length = kenv_gen->kmer_offsets[refcode + 1] -
                                      kenv_gen->kmer_offsets[refcode];
    }
    return kenv_gen->spaceGtKenvRawElement + kenv_gen->kmer_offsets[refcode];
  }
}

int gt_kenv_eval_score(const GtKenvGenerator *kenv_gen,
                       const GtUchar *a_encoded,
                       const GtUchar *b_encoded,
                       GtUword length)
{
  int score_sum = 0;
  GtUword idx;

  for (idx = 0; idx < length; idx++)
  {
    score_sum += gt_score_matrix_get_score(kenv_gen->score_matrix,
                                           a_encoded[idx],b_encoded[idx]);
  }
  return score_sum;
}
