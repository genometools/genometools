/*
  Copyright (c) 2006-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/fasta.h"
#include "core/fasta_separator.h"
#include "core/xansi_api.h"

void gt_fasta_show_entry(const char *description, const char *sequence,
                         GtUword sequence_length, GtUword width,
                         GtFile *outfp)
{
  gt_fasta_show_entry_with_suffix(description, sequence, sequence_length, NULL,
                                  width, outfp);
}

void gt_fasta_show_entry_str(const char *description, const char *sequence,
                             GtUword sequence_length, GtUword width,
                             GtStr *outstr)
{
  gt_fasta_show_entry_with_suffix_str(description, sequence, sequence_length,
                                      NULL, width, outstr);
}

void gt_fasta_show_entry_with_suffix(const char *description,
                                     const char *sequence,
                                     GtUword sequence_length,
                                     const char *suffix, GtUword width,
                                     GtFile *outfp)
{
  GtUword i, current_length, suffix_length;
  gt_assert(sequence);
  gt_file_xfputc(GT_FASTA_SEPARATOR, outfp);
  if (description)
    gt_file_xfputs(description, outfp);
  gt_file_xfputc('\n', outfp);
  suffix_length = suffix ? strlen(suffix) : 0;
  for (i = 0, current_length = 0; i < sequence_length + suffix_length;
       i++, current_length++) {
    if (width && current_length == width) {
      gt_file_xfputc('\n', outfp);
      current_length = 0;
    }
    if (i < sequence_length)
      gt_file_xfputc(sequence[i], outfp);
    else
      gt_file_xfputc(suffix[i-sequence_length], outfp);
  }
  gt_file_xfputc('\n', outfp);
}

void gt_fasta_show_entry_with_suffix_str(const char *description,
                                         const char *sequence,
                                         GtUword sequence_length,
                                         const char *suffix, GtUword width,
                                         GtStr *outstr)
{
  GtUword i, current_length, suffix_length;
  gt_assert(sequence && outstr);
  gt_str_append_char(outstr, GT_FASTA_SEPARATOR);
  if (description)
    gt_str_append_cstr(outstr, description);
  gt_str_append_char(outstr, '\n');
  suffix_length = suffix ? strlen(suffix) : 0;
  for (i = 0, current_length = 0; i < sequence_length + suffix_length;
       i++, current_length++) {
    if (width && current_length == width) {
      gt_str_append_char(outstr, '\n');
      current_length = 0;
    }
    if (i < sequence_length)
      gt_str_append_char(outstr, sequence[i]);
    else
      gt_str_append_char(outstr, suffix[i-sequence_length]);
  }
  gt_str_append_char(outstr, '\n');
}
