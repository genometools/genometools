/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include "libgtmatch/eis-blockcomp-param.h"
#include "libgtmatch/eis-encidxseq.h"

extern void
registerBlockEncOptions(OptionParser *op, struct blockEncParams *paramOutput)
{
  Option *option;

  option = option_new_uint_min("bsize", "specify size of blocks",
                               &paramOutput->blockSize, 8U, 1U);
  option_parser_add_option(op, option);
  option = option_new_uint_min("blbuck", "specify number of blocks per bucket",
                               &paramOutput->bucketBlocks, 8U, 1U);
  option_parser_add_option(op, option);
}
