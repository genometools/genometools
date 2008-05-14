/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>
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

#include "libgtcore/option.h"
#include "libgtmatch/eis-bwtseq-context-param.h"

extern void
registerCtxMapOptions(OptionParser *op, int *ilogOut)
{
  Option *option = option_new_int_min_max(
    "ctxilog", "specify the interval of context sampling as log value\n"
    "parameter i means that each 2^i-th position of source is sampled for "
    "rank\n-1 => chooses default of log(log(sequence length))\n"
    "-2 => generates no map",
    ilogOut, CTX_MAP_ILOG_NOMAP, CTX_MAP_ILOG_AUTOSIZE,
    sizeof (Seqpos) * CHAR_BIT - 1);
  option_parser_add_option(op, option);
}
