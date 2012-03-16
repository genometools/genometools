/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "core/assert_api.h"
#include "core/log_api.h"
#include "core/xansi_api.h"
#include "extended/bitoutstream.h"

struct GtBitOutStream {
  unsigned long bits_left;
  GtBitsequence bitseqbuffer;
  FILE *fp;
  long pagesize;
};

GtBitOutStream* gt_bitoutstream_new(FILE *fp)
{
  GtBitOutStream *b;
  b = gt_calloc((size_t) 1, sizeof (GtBitOutStream));
  b->bits_left = GT_INTWORDSIZE;
  b->fp = fp;
  b->bitseqbuffer = 0;
  b->pagesize = sysconf(_SC_PAGESIZE);
  return b;
}

int gt_bitoutstream_append(GtBitOutStream *b,
                           GtBitsequence code,
                           unsigned long bits2write,
                           GT_UNUSED GtError *err)
{
  int had_err = 0;

  if (b->bits_left < bits2write) {
    unsigned int overhang = 0;
    overhang = bits2write - b->bits_left;
    b->bitseqbuffer |= code >> overhang;
    gt_xfwrite(&b->bitseqbuffer,
               sizeof (GtBitsequence),
               1, b->fp);
    b->bitseqbuffer = 0;
    b->bits_left = GT_INTWORDSIZE - overhang;
  }
  else {
    b->bits_left -= bits2write;
  }
  b->bitseqbuffer |= code << b->bits_left;
  return had_err;
}

int gt_bitoutstream_flush(GtBitOutStream *b,
                          GT_UNUSED GtError *err)
{
  gt_xfwrite(&b->bitseqbuffer, sizeof (GtBitsequence), 1, b->fp);

  b->bitseqbuffer = 0;
  b->bits_left = GT_INTWORDSIZE;
  return 0;
}

int gt_bitoutstream_flush_and_jump_to_next_page(GtBitOutStream *b,
                                                GtError *err)
{
  unsigned long fpos;
  gt_bitoutstream_flush(b, err);
  if (!(ftell(b->fp) % b->pagesize))
    return 0;
  fpos = (ftell(b->fp) / b->pagesize + 1) * b->pagesize;
  gt_xfseek(b->fp, fpos, SEEK_SET);
  return 0;
}

off_t gt_bitoutstream_pos(GtBitOutStream *b)
{
  return (off_t) ftell(b->fp);
}

void gt_bitoutstream_delete(GtBitOutStream *b)
{
  gt_free(b);
}
