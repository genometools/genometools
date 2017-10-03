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

#include "core/ma_api.h"
#include "match/initbasepower.h"
#include "k_env_alphabet.h"

/***** GtKenvAlphabet *****/
struct GtKenvAlphabet
{
  GtAlphabet *alphabet;
  unsigned int alph_size;
  GtCodetype *alph_powers;
  bool is_dna;
};

GtKenvAlphabet *gt_kenv_alphabet_new(bool use_dna, unsigned int q_val)
{
  GtKenvAlphabet *kenv_alph = gt_malloc(sizeof(*kenv_alph));

  if (use_dna)
    kenv_alph->alphabet = gt_alphabet_new_dna();
  else
    kenv_alph->alphabet = gt_alphabet_new_protein();

  kenv_alph->alph_size = gt_alphabet_num_of_chars(kenv_alph->alphabet);
  kenv_alph->alph_powers = gt_initbasepower(kenv_alph->alph_size, q_val);
  kenv_alph->is_dna = use_dna;
  return kenv_alph;
}

void gt_kenv_alphabet_delete(GtKenvAlphabet *kenv_alph)
{
  if (kenv_alph != NULL)
  {
    if (kenv_alph->alphabet != NULL)
      gt_alphabet_delete(kenv_alph->alphabet);
    if (kenv_alph->alph_powers != NULL)
      gt_free(kenv_alph->alph_powers);
    gt_free(kenv_alph);
  }
}

GtAlphabet *gt_kenv_alphabet_get_alphabet(const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph && kenv_alph && kenv_alph->alphabet);
  return kenv_alph->alphabet;
}

const GtCodetype *gt_kenv_alphabet_get_alph_powers(
                                            const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph && kenv_alph->alph_powers);
  return kenv_alph->alph_powers;
}

unsigned int gt_kenv_alphabet_get_alph_size(const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph);
  return kenv_alph->alph_size;
}

bool gt_kenv_alphabet_is_dna(const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph);
  return kenv_alph->is_dna;
}

const GtUchar* gt_kenv_alphabet_get_symbolmap(const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph && kenv_alph->alphabet);
  return gt_alphabet_symbolmap(kenv_alph->alphabet);
}

const GtUchar* gt_kenv_alphabet_get_characters(const GtKenvAlphabet *kenv_alph)
{
  gt_assert(kenv_alph);
  return gt_alphabet_characters(kenv_alph->alphabet);
}
