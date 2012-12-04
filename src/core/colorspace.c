/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/colorspace.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"

typedef struct GtColorSpaceDecoder GtColorSpaceDecoder;

struct color_state {
  char self;
  struct color_state *links[5];
};

struct GtColorSpaceDecoder {
  struct color_state *current,
                     state_A,
                     state_C,
                     state_G,
                     state_T,
                     state_N;
};

static GtColorSpaceDecoder *gt_colorspace_decoder_new(void)
{
  GtColorSpaceDecoder *cd = gt_malloc(sizeof(*cd));
  cd->current = NULL;
  cd->state_A.self = 'a';
  cd->state_C.self = 'c';
  cd->state_G.self = 'g';
  cd->state_T.self = 't';
  cd->state_N.self = 'n';
  cd->state_A.links[0] = &(cd->state_A);
  cd->state_A.links[1] = &(cd->state_C);
  cd->state_A.links[2] = &(cd->state_G);
  cd->state_A.links[3] = &(cd->state_T);
  cd->state_A.links[4] = &(cd->state_N);
  cd->state_C.links[0] = &(cd->state_C);
  cd->state_C.links[1] = &(cd->state_A);
  cd->state_C.links[2] = &(cd->state_T);
  cd->state_C.links[3] = &(cd->state_G);
  cd->state_C.links[4] = &(cd->state_N);
  cd->state_G.links[0] = &(cd->state_G);
  cd->state_G.links[1] = &(cd->state_T);
  cd->state_G.links[2] = &(cd->state_A);
  cd->state_G.links[3] = &(cd->state_C);
  cd->state_G.links[4] = &(cd->state_N);
  cd->state_T.links[0] = &(cd->state_T);
  cd->state_T.links[1] = &(cd->state_G);
  cd->state_T.links[2] = &(cd->state_C);
  cd->state_T.links[3] = &(cd->state_A);
  cd->state_T.links[4] = &(cd->state_N);
  cd->state_N.links[0] =
  cd->state_N.links[1] =
  cd->state_N.links[2] =
  cd->state_N.links[3] =
  cd->state_N.links[4] = &(cd->state_N);
  return cd;
}

static void gt_colorspace_decoder_delete(GtColorSpaceDecoder *cd)
{
  if (!cd) return;
  gt_free(cd);
}

static int set_first_state(GtColorSpaceDecoder *cd, char start, GtError *err)
{
  int had_err = 0;

  gt_assert(cd->current == NULL);

  switch (start)
  {
    case 'a':
    case 'A':
      cd->current = &(cd->state_A);
      break;
    case 'c':
    case 'C':
      cd->current = &(cd->state_C);
      break;
    case 'g':
    case 'G':
      cd->current = &(cd->state_G);
      break;
    case 't':
    case 'T':
      cd->current = &(cd->state_T);
      break;
    default:
      gt_error_set(err, "encountered wrong start character while encoding "
                        "color space string: %c!\n", start);
      had_err = -1;
  }
  return had_err;
}

static int set_next_state(GtColorSpaceDecoder *cd, char next, GtError *err)
{
  int had_err = 0;

  gt_assert(cd->current != NULL);

  switch (next)
  {
    case '0':
    case '1':
    case '2':
    case '3':
      cd->current = cd->current->links[next - 48];
      break;
    case '.':
    case '4':
    case '5':
      cd->current = cd->current->links[4];
      break;
    default:
      gt_error_set(err, "encountered wrong character while encoding color "
                        "space string: %c\n", next);
      had_err = -1;
  }
  return had_err;
}

int gt_colorspace_decode_string(GtStr *color_string, GtStr *result,
                                GtError *err)
{
  int had_err = 0;
  unsigned long str_len, idx;
  char *input;
  GtColorSpaceDecoder *cd;

  cd = gt_colorspace_decoder_new();

  gt_assert(cd);
  gt_assert(color_string && result);
  gt_assert(gt_str_length(result) == 0);

  input = gt_str_get(color_string);
  str_len = gt_str_length(color_string);

  had_err = set_first_state(cd, input[0], err);
  if (!had_err)
  {
    gt_str_append_char(result, cd->current->self);

    for (idx = 1UL; !had_err && idx < str_len; idx++)
    {
      had_err = set_next_state(cd, input[idx], err);
      gt_str_append_char(result, cd->current->self);
    }
    gt_assert(gt_str_length(result) == str_len);
  }
  gt_colorspace_decoder_delete(cd);
  return had_err;
}

int gt_colorspace_unit_test(GtError *err)
{
  int had_err = 0,
      test_err = 0;
  GtStr *in_gt = gt_str_new_cstr("A0112233122331211332311223"),
        *con_gt = gt_str_new_cstr("AACAGATACTCGCAGTGCGATGTCGA"),
        *output = gt_str_new();

  had_err = gt_colorspace_decode_string(in_gt,
                                        output,
                                        err);
  gt_ensure(had_err, gt_str_cmp(con_gt, output));

  if (!had_err)
  {
    GtError *my_err = gt_error_new();

    gt_str_delete(in_gt);
    gt_str_reset(output);
    in_gt = gt_str_new_cstr("X0011");

    test_err = gt_colorspace_decode_string(in_gt,
                                           output,
                                           my_err);

    gt_ensure(had_err, test_err);
    gt_ensure(had_err, gt_error_is_set(my_err));

    gt_error_delete(my_err);
    test_err = 0;
  }
  if (!had_err)
  {
    GtError *my_err = gt_error_new();

    gt_str_delete(in_gt);
    gt_str_reset(output);
    in_gt = gt_str_new_cstr("a0011a");

    test_err = gt_colorspace_decode_string(in_gt,
                                           output,
                                           my_err);

    gt_ensure(had_err, test_err);
    gt_ensure(had_err, gt_error_is_set(my_err));

    gt_error_delete(my_err);
    test_err = 0;
  }
  gt_str_delete(in_gt);
  gt_str_delete(con_gt);
  gt_str_delete(output);
  return had_err;
}
