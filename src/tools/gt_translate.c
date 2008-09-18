/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/option.h"
#include "core/ma.h"
#include "core/translate.h"
#include "core/versionfunc.h"
#include "tools/gt_translate.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  OPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("DNA_sequence", "Translate DNA sequence in all the "
                         "reading frames.");
  gt_option_parser_set_min_max_args(op, 1, 1);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_translate(int argc, const char **argv, GtError *err)
{
  char *frame1 = NULL,
       *frame2 = NULL,
       *frame3 = NULL;
  unsigned long seqlen;
  const char *seq;
  int parsed_args;
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* translation */
  seq = argv[parsed_args];
  seqlen = strlen(seq);

  if (seqlen >= 3) {
    gt_translate_all_frames(&frame1, &frame2, &frame3, seq, seqlen);
    printf("frame1: %s\nframe2: %s\nframe3: %s\n", frame1, frame2, frame3);
  }

  gt_free(frame1);
  gt_free(frame2);
  gt_free(frame3);

  return 0;
}
