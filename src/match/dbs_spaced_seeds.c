/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <inttypes.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "core/codetype.h"
#include "core/unused_api.h"
#include "match/dbs_spaced_seeds.h"

int gt_spaced_seed_span(GtCodetype spaced_seed)
{
  int span = 0;

  for (/* Nothing */; spaced_seed > 0; spaced_seed >>= 1)
  {
    span++;
  }
  return span;
}

int gt_spaced_seed_weight(GtCodetype spaced_seed)
{
  int weight = 0;

  for (/* Nothing */; spaced_seed > 0; spaced_seed >>= 1)
  {
    if (spaced_seed & (GtCodetype) 1)
    {
      weight++;
    }
  }
  return weight;
}

static GtCodetype gt_spaced_seed_spec_tab[] = {
  23075ULL /* 7, 15 */,
  29331ULL /* 8 */,
  27975ULL /* 9 */,
  27823ULL /* 10 */,
  30135ULL /* 11 */,
  30575ULL /* 12 */,
  32495ULL /* 13 */,
  32511ULL /* 14 */,
  39559ULL /* 8, 16 */,
  54039ULL /* 9 */,
  55511ULL /* 10 */,
  59767ULL /* 11 */,
  56687ULL /* 12 */,
  63215ULL /* 13 */,
  64479ULL /* 14 */,
  65471ULL /* 15 */,
  100891ULL /* 8, 17 */,
  108075ULL /* 9 */,
  111271ULL /* 10 */,
  119415ULL /* 11 */,
  125751ULL /* 12 */,
  122287ULL /* 13 */,
  128879ULL /* 14 */,
  128959ULL /* 15 */,
  130943ULL /* 16 */,
  217383ULL /* 9, 18 */,
  234071ULL /* 10 */,
  238903ULL /* 11 */,
  240951ULL /* 12 */,
  251503ULL /* 13 */,
  256887ULL /* 14 */,
  259823ULL /* 15 */,
  261087ULL /* 16 */,
  262015ULL /* 17 */,
  412715ULL /* 9, 19 */,
  469271ULL /* 10 */,
  469399ULL /* 11 */,
  469615ULL /* 12 */,
  486575ULL /* 13 */,
  504751ULL /* 14 */,
  513775ULL /* 15 */,
  507359ULL /* 16 */,
  520127ULL /* 17 */,
  523263ULL /* 18 */,
  860951ULL /* 10, 20 */,
  893607ULL /* 11 */,
  995927ULL /* 12 */,
  963375ULL /* 13 */,
  1009327ULL /* 14 */,
  1029039ULL /* 15 */,
  1027551ULL /* 16 */,
  1031647ULL /* 17 */,
  1040255ULL /* 18 */,
  1048319ULL /* 19 */,
  1902795ULL /* 10, 21 */,
  1739175ULL /* 11 */,
  1880663ULL /* 12 */,
  1992015ULL /* 13 */,
  1952559ULL /* 14 */,
  1955487ULL /* 15 */,
  2055031ULL /* 16 */,
  2060015ULL /* 17 */,
  2080223ULL /* 18 */,
  2080511ULL /* 19 */,
  2095103ULL /* 20 */,
  3754263ULL /* 11, 22 */,
  3969703ULL /* 12 */,
  3970407ULL /* 13 */,
  3847375ULL /* 14 */,
  3905119ULL /* 15 */,
  3909487ULL /* 16 */,
  4110063ULL /* 17 */,
  4126447ULL /* 18 */,
  4176863ULL /* 19 */,
  4177791ULL /* 20 */,
  4193791ULL /* 21 */,
  7508247ULL /* 11, 23 */,
  7490215ULL /* 12 */,
  7950951ULL /* 13 */,
  7956055ULL /* 14 */,
  7951983ULL /* 15 */,
  8074607ULL /* 16 */,
  8219887ULL /* 17 */,
  8220399ULL /* 18 */,
  8240607ULL /* 19 */,
  8320991ULL /* 20 */,
  8355583ULL /* 21 */,
  8387583ULL /* 22 */,
  14848567ULL /* 12, 24 */,
  15280743ULL /* 13 */,
  15911479ULL /* 14 */,
  15912111ULL /* 15 */,
  16149199ULL /* 16 */,
  16174767ULL /* 17 */,
  16469743ULL /* 18 */,
  16217535ULL /* 19 */,
  16629215ULL /* 20 */,
  16644031ULL /* 21 */,
  16760703ULL /* 22 */,
  16776191ULL /* 23 */,
  28387495ULL /* 12, 25 */,
  31755435ULL /* 13 */,
  30624311ULL /* 14 */,
  31019727ULL /* 15 */,
  30775663ULL /* 16 */,
  32872879ULL /* 17 */,
  32303839ULL /* 18 */,
  32988063ULL /* 19 */,
  33222127ULL /* 20 */,
  32996319ULL /* 21 */,
  33283967ULL /* 22 */,
  33488639ULL /* 23 */,
  33546239ULL /* 24 */,
  61019287ULL /* 13, 26 */,
  62007631ULL /* 14 */,
  62178639ULL /* 15 */,
  64578391ULL /* 16 */,
  64330095ULL /* 17 */,
  65755551ULL /* 18 */,
  65756383ULL /* 19 */,
  65894255ULL /* 20 */,
  66022335ULL /* 21 */,
  66026431ULL /* 22 */,
  66576127ULL /* 23 */,
  66977279ULL /* 24 */,
  67092479ULL /* 25 */,
  126495003ULL /* 13, 27 */,
  122309719ULL /* 14 */,
  124131927ULL /* 15 */,
  124308175ULL /* 16 */,
  124954271ULL /* 17 */,
  129160607ULL /* 18 */,
  128896367ULL /* 19 */,
  131786479ULL /* 20 */,
  131784159ULL /* 21 */,
  131984863ULL /* 22 */,
  133151711ULL /* 23 */,
  133685183ULL /* 24 */,
  133954559ULL /* 25 */,
  134201343ULL /* 26 */,
  244880463ULL /* 14, 28 */,
  254945615ULL /* 15 */,
  256519375ULL /* 16 */,
  255145071ULL /* 17 */,
  249914783ULL /* 18 */,
  262878623ULL /* 19 */,
  263615855ULL /* 20 */,
  263579375ULL /* 21 */,
  263634399ULL /* 22 */,
  264092639ULL /* 23 */,
  266303423ULL /* 24 */,
  267378559ULL /* 25 */,
  267909119ULL /* 26 */,
  268402687ULL /* 27 */,
  508768943ULL /* 17, 29 */,
  513435311ULL /* 18 */,
  499881567ULL /* 19 */,
  525769951ULL /* 20 */,
  527260911ULL /* 21 */,
  527674815ULL /* 22 */,
  527920095ULL /* 23 */,
  532134879ULL /* 24 */,
  534640575ULL /* 25 */,
  534740735ULL /* 26 */,
  536345599ULL /* 27 */,
  1051087767ULL /* 18, 30 */,
  1028869743ULL /* 19 */,
  1047213423ULL /* 20 */,
  1035629407ULL /* 21 */,
  1054521823ULL /* 22 */,
  1055717055ULL /* 23 */,
  1055878079ULL /* 24 */,
  1056373695ULL /* 25 */,
  1065220031ULL /* 26 */,
  1069514239ULL /* 27 */,
  1073479167ULL /* 28 */,
  2040932015ULL /* 18, 31 */,
  2057774495ULL /* 19 */,
  2067064495ULL /* 20 */,
  2103078127ULL /* 21 */,
  2071258815ULL /* 22 */,
  2104348351ULL /* 23 */,
  2111548911ULL /* 24 */,
  2126216159ULL /* 25 */,
  2130115519ULL /* 26 */,
  2138828735ULL /* 27 */,
  2143223551ULL /* 28 */,
  2147220991ULL /* 29 */,
  4207733599ULL /* 22, 32 */,
  4208813935ULL /* 23 */,
  4218133983ULL /* 24 */,
  4225429215ULL /* 25 */,
  4223523807ULL /* 26 */,
  4226775999ULL /* 27 */,
  4260872063ULL /* 28 */,
  4286545663ULL /* 29 */,
  4292868095ULL /* 30 */
};

static int gt_spaced_seed_span_start_tab[] = {
  0, 8, 16, 25, 34, 44, 54, 65, 76, 88, 100, 113, 126, 140, 154, 165, 176, 188
};

static int gt_spaced_seed_first_weight_tab[] = {
  7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 17, 18, 18, 22
};

void gt_spaced_seed_weight_range(int *min_weight,int *max_weight, int span)
{
  const int size_table = sizeof gt_spaced_seed_spec_tab/
                         sizeof gt_spaced_seed_spec_tab[0];
  const int num_spans = sizeof gt_spaced_seed_span_start_tab/
                        sizeof gt_spaced_seed_span_start_tab[0];
  int span_offset, end_of_span;

  gt_assert(GT_SPACED_SEED_FIRST_SPAN <= span &&
            span - GT_SPACED_SEED_FIRST_SPAN <= num_spans - 1);
  *min_weight
    = gt_spaced_seed_first_weight_tab[span - GT_SPACED_SEED_FIRST_SPAN];
  if (span - GT_SPACED_SEED_FIRST_SPAN == num_spans - 1)
  {
    end_of_span = size_table;
  } else
  {
    end_of_span
      = gt_spaced_seed_span_start_tab[span + 1 - GT_SPACED_SEED_FIRST_SPAN];
  }
  span_offset = gt_spaced_seed_span_start_tab[span - GT_SPACED_SEED_FIRST_SPAN];
  *max_weight = *min_weight + end_of_span - span_offset - 1;
}

typedef struct
{
  GtCodetype extract;
  int shiftright;
} GtSpacedSeedSpecValue;

struct GtSpacedSeedSpec
{
  GtSpacedSeedSpecValue *spec_tab;
  size_t num_specs;
};

GtSpacedSeedSpec *gt_spaced_seed_spec_new(GtCodetype spacedseed)
{
  uint8_t blocks_length[32] = {0}, shiftleft = 0, shiftright = 0;
  GtCodetype ss_copy, last = (GtCodetype) 1, GT_UNUSED from_blocks = 0;
  GtUword idx, block_num = 0, spec_counter = 0;
  GtSpacedSeedSpec *seed_spec;

  gt_assert(spacedseed & ((GtCodetype) 1));
  blocks_length[0] = 1;
  for (ss_copy = spacedseed >> 1; ss_copy > 0; ss_copy >>= 1)
  {
    GtCodetype current = ss_copy & (GtCodetype) 1;

    if (current != last)
    {
      block_num++;
      last = current;
    }
    gt_assert(block_num < sizeof blocks_length/sizeof blocks_length[0]);
    blocks_length[block_num]++;
  }
  block_num++;
  gt_assert(block_num % 2 == 1);
  seed_spec = gt_malloc(sizeof *seed_spec);
  seed_spec->num_specs = 1 + block_num/2;
  seed_spec->spec_tab = gt_malloc(seed_spec->num_specs *
                                  sizeof *seed_spec->spec_tab);
  gt_assert(seed_spec->spec_tab != NULL);
  last = (GtCodetype) 1;
  for (idx = 0; idx < block_num; idx++)
  {
    const uint8_t width = blocks_length[idx];
    if (last == (GtCodetype) 1)
    {
      from_blocks |= (((((GtCodetype) 1) << width) - 1) << shiftleft);
      last = 0;
      gt_assert(spec_counter < seed_spec->num_specs);
      seed_spec->spec_tab[spec_counter].extract
        = ((((GtCodetype) 1) << (2 * width)) - 1) << (2 * shiftleft);
      seed_spec->spec_tab[spec_counter++].shiftright = 2 * shiftright;
    } else
    {
      last = (GtCodetype) 1;
      shiftright += width;
    }
    shiftleft += width;
  }
  gt_assert(spec_counter == seed_spec->num_specs && from_blocks == spacedseed);
  return seed_spec;
}

void gt_spaced_seed_spec_delete(GtSpacedSeedSpec *seed_spec)
{
  if (seed_spec != NULL)
  {
    gt_free(seed_spec->spec_tab);
    gt_free(seed_spec);
  }
}

static int gt_spaced_seed_tab_num_extract(int weight,int span)
{
  int first_weight, span_offset;

  gt_assert(GT_SPACED_SEED_FIRST_SPAN <= span);
  span_offset = gt_spaced_seed_span_start_tab[span - GT_SPACED_SEED_FIRST_SPAN];
  first_weight
    = gt_spaced_seed_first_weight_tab[span - GT_SPACED_SEED_FIRST_SPAN];
  gt_assert(first_weight <= weight);
  return span_offset + weight - first_weight;
}

GtSpacedSeedSpec *gt_spaced_seed_spec_new_from_ws(int weight,int span)
{
  int min_weight, max_weight, seed_num;;

  gt_spaced_seed_weight_range(&min_weight,&max_weight, span);
  gt_assert(min_weight <= weight && weight <= max_weight);
  seed_num = gt_spaced_seed_tab_num_extract(weight,span);
  return gt_spaced_seed_spec_new(gt_spaced_seed_spec_tab[seed_num]);
}

GtCodetype gt_spaced_seed_extract_generic(const GtSpacedSeedSpec *seed_spec,
                                          GtCodetype kmer)
{
  GtCodetype ext = 0;
  const GtSpacedSeedSpecValue *seed_spec_ptr;

  gt_assert(seed_spec != NULL);
  for (seed_spec_ptr = seed_spec->spec_tab;
       seed_spec_ptr < seed_spec->spec_tab + seed_spec->num_specs;
       seed_spec_ptr++)
  {
    ext |= ((kmer & seed_spec_ptr->extract) >> seed_spec_ptr->shiftright);
  }
  return ext;
}
