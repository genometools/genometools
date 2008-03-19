/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/bioseq.h"
#include "libgtcore/fasta.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/stringdistri.h"
#include "libgtcore/unused.h"
#include "tools/gt_sequniq.h"

static OptionParser* gt_sequniq_option_parser_new(UNUSED
                                                      void *tool_arguments)
{
  OptionParser *op;
  op = option_parser_new("[option ...] sequence_file [...] ",
                         "Filter out repeated sequences in given in given "
                         "sequence_file(s).");
  option_parser_set_min_args(op, 1);
  return op;
}

static int gt_sequniq_runner(int argc, const char **argv,
                                 UNUSED void *tool_arguments, Error *err)
{
  Bioseq *bs;
  StringDistri *sd;
  unsigned long i, j;
  unsigned long long duplicates = 0, num_of_sequences = 0;
  int had_err = 0;

  error_check(err);
  sd = stringdistri_new();

  for (i = 0; !had_err && i < argc; i++) {
    if (!(bs = bioseq_new(argv[i], err)))
      had_err = -1;
    if (!had_err) {
      for (j = 0; j < bioseq_number_of_sequences(bs); j++) {
        if (!stringdistri_get(sd, bioseq_get_md5_fingerprint(bs, j))) {
          stringdistri_add(sd, bioseq_get_md5_fingerprint(bs, j));
          fasta_show_entry(bioseq_get_description(bs, j),
                           bioseq_get_sequence(bs, j),
                           bioseq_get_sequence_length(bs, j), 0);
        }
        else
          duplicates++;
        num_of_sequences++;
      }
    }
    bioseq_delete(bs);
  }

  /* show statistics */
  if (!had_err) {
    fprintf(stderr, "# %llu out of %llu sequences have been removed (%.3f%%)\n",
            duplicates, num_of_sequences,
            ((double) duplicates / num_of_sequences) * 100.0);
  }

  stringdistri_delete(sd);

  return had_err;
}

Tool* gt_sequniq(void)
{
  return tool_new(NULL,
                  NULL,
                  gt_sequniq_option_parser_new,
                  NULL,
                  gt_sequniq_runner);
}
