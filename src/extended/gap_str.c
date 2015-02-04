/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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
#include "core/cstr_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/parseutils_api.h"
#include "core/splitter_api.h"
#include "extended/gap_str.h"

struct GtGapStr
{
  GtUword ref_len,
          tar_len,
          ali_len,
          step;
  bool is_protein_alignment;
};

static int gt_gap_str_parse(GtGapStr *gap_str, const char *str, GtError *err)
{
  int had_err = 0;
  char *mystr = gt_cstr_dup(str),
       *token = NULL,
       code;
  GtUword nof_tokens, i;
  GtSplitter *splitter;
  gt_assert(gap_str && str);
  gt_error_check(err);

  splitter = gt_splitter_new();
  gt_splitter_reset(splitter);
  gt_splitter_split(splitter, mystr, strlen(mystr), ' ');
  nof_tokens = gt_splitter_size(splitter);

  for (i = 0; !had_err && i < nof_tokens; i++) {
    int toklen;
    token = gt_splitter_get_token(splitter, i);
    toklen = strlen(token);
    if (toklen < 2) {
      if (token[0] != ' ' && token[0] != '\n') {
        gt_error_set(err, "edit operation too short: \"%s\"", token);
        had_err = -1;
        break;
      }
    }
    if (!had_err) {
      gt_assert(toklen > 1);
      switch (token[0]) {
        case 'M':
        case 'I':
        case 'D':
          break;
        case 'F':
        case 'R':
          if (!gap_str->is_protein_alignment) {
            gt_error_set(err, "invalid edit operation code '%c' (only allowed "
                              "in nucleotide-protein matches)",
                         token[0]);
            had_err = -1;
          }
          break;
        default:
          gt_error_set(err, "invalid edit operation code '%c' encountered",
                       token[0]);
          had_err = -1;
          break;
      }
    }
    if (!had_err) {
      GtUword oplen;
      code = (token++)[0];
      if (0 != gt_parse_uword(&oplen, token)) {
        gt_error_set(err, "cannot parse edit length from string: \"%s\"",
                     token);
        had_err = -1;
        break;
      }
      switch (code) {
        case 'M':
          gap_str->ref_len += gap_str->step * oplen;
          gap_str->ali_len += gap_str->step * oplen;
          gap_str->tar_len += gap_str->step * oplen;
          break;
        case 'I':
          gap_str->ali_len += gap_str->step * oplen;
          gap_str->tar_len += gap_str->step * oplen;
          break;
        case 'D':
          gap_str->ref_len += gap_str->step * oplen;
          gap_str->ali_len += gap_str->step * oplen;
          break;
        case 'F':
          gap_str->ref_len += oplen;
          gap_str->ali_len += oplen;
          break;
        case 'R':
          if (gap_str->ref_len < oplen || gap_str->ali_len < oplen) {
            gt_error_set(err, "reverse frameshift is too long (" GT_WU ")",
                         oplen);
            had_err = -1;
          }
          gap_str->ref_len -= oplen;
          gap_str->ali_len -= oplen;
          break;
        default:
        gt_assert(false);   /* cannot happen */
      }
    }
  }

  gt_free(mystr);
  gt_splitter_delete(splitter);
  return had_err;
}

static GtGapStr *gt_gap_str_new(const char *str, bool is_protein, GtError *err)
{
  GtGapStr *gap_str;
  gap_str = gt_malloc(sizeof (GtGapStr));
  gap_str->is_protein_alignment = is_protein;
  gap_str->step = is_protein ? 3 : 1;
  gap_str->ali_len = 0;
  gap_str->ref_len = 0;
  gap_str->tar_len = 0;
  if (gt_gap_str_parse(gap_str, str, err) != 0) {
    gt_free(gap_str);
    return NULL;
  }
  return gap_str;
}

GtGapStr*  gt_gap_str_new_nucleotide(const char *str, GtError *err)
{
  return gt_gap_str_new(str, false, err);
}

GtGapStr*  gt_gap_str_new_protein(const char *str, GtError *err)
{
  return gt_gap_str_new(str, true, err);
}

GtUword gt_gap_str_length_alignment(const GtGapStr *gap_str)
{
  gt_assert(gap_str);
  return gap_str->ali_len;
}

GtUword gt_gap_str_length_reference(const GtGapStr *gap_str)
{
  gt_assert(gap_str);
  return gap_str->ref_len;
}

GtUword gt_gap_str_length_target(const GtGapStr *gap_str)
{
  gt_assert(gap_str);
  return gap_str->tar_len;
}

void gt_gap_str_delete(GtGapStr *gap_str)
{
  gt_free(gap_str);
}
