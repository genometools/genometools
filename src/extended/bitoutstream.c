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
  GtBitOutStream *bitstream;
  gt_assert(fp);

  bitstream = gt_calloc((size_t) 1, sizeof (GtBitOutStream));
  bitstream->bits_left = GT_INTWORDSIZE;
  bitstream->fp = fp;
  bitstream->bitseqbuffer = 0;
  bitstream->pagesize = sysconf(_SC_PAGESIZE);
  return bitstream;
}

void gt_bitoutstream_append(GtBitOutStream *bitstream,
                            GtBitsequence code,
                            unsigned long bits2write)
{
  if (bitstream->bits_left < bits2write) {
    unsigned int overhang = 0;
    overhang = bits2write - bitstream->bits_left;
    bitstream->bitseqbuffer |= code >> overhang;
    gt_xfwrite(&bitstream->bitseqbuffer,
               sizeof (GtBitsequence),
               1, bitstream->fp);
    bitstream->bitseqbuffer = 0;
    bitstream->bits_left = GT_INTWORDSIZE - overhang;
  }
  else {
    bitstream->bits_left -= bits2write;
  }
  bitstream->bitseqbuffer |= code << bitstream->bits_left;
}

void gt_bitoutstream_flush(GtBitOutStream *bitstream)
{
  gt_assert(bitstream);
  gt_xfwrite(&bitstream->bitseqbuffer, sizeof (GtBitsequence),
             1, bitstream->fp);

  bitstream->bitseqbuffer = 0;
  bitstream->bits_left = GT_INTWORDSIZE;
}

void gt_bitoutstream_flush_advance(GtBitOutStream *bitstream)
{
  unsigned long fpos;
  gt_assert(bitstream);
  gt_bitoutstream_flush(bitstream);
  if ((ftell(bitstream->fp) % bitstream->pagesize) != 0) {
    fpos = (ftell(bitstream->fp) / bitstream->pagesize + 1) *
           bitstream->pagesize;
    gt_xfseek(bitstream->fp, fpos, SEEK_SET);
  }
}

off_t gt_bitoutstream_pos(const GtBitOutStream *bitstream)
{
  return (off_t) ftell(bitstream->fp);
}

void gt_bitoutstream_delete(GtBitOutStream *bitstream)
{
  gt_free(bitstream);
}
