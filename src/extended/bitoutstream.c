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

#ifndef S_SPLINT_S
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "core/assert_api.h"
#include "core/log_api.h"
#include "core/xansi_api.h"
#include "extended/bitoutstream.h"

struct GtBitOutStream {
  int bits_left;
  long pagesize;
  unsigned long written_bits;
  GtBitsequence bitseqbuffer;
  FILE *fp;
};

GtBitOutStream* gt_bitoutstream_new(FILE *fp)
{
  GtBitOutStream *bitstream;
  gt_assert(fp);

  bitstream = gt_calloc((size_t) 1, sizeof (GtBitOutStream));
  bitstream->bits_left = GT_INTWORDSIZE;
  bitstream->fp = fp;
  bitstream->bitseqbuffer = 0;
  bitstream->written_bits = 0;
  bitstream->pagesize = sysconf((int) _SC_PAGESIZE);
  return bitstream;
}

void gt_bitoutstream_append(GtBitOutStream *bitstream,
                            GtBitsequence code,
                            unsigned bits_to_write)
{
  if ((unsigned) bitstream->bits_left < bits_to_write) {
    unsigned overhang = bits_to_write - bitstream->bits_left;
    bitstream->bitseqbuffer |= code >> overhang;
    gt_xfwrite(&bitstream->bitseqbuffer,
               sizeof (GtBitsequence),
               (size_t) 1, bitstream->fp);
    bitstream->bitseqbuffer = 0;
    bitstream->bits_left = GT_INTWORDSIZE - overhang;
    bitstream->written_bits += GT_INTWORDSIZE;
  }
  else {
    bitstream->bits_left -= bits_to_write;
  }
  bitstream->bitseqbuffer |= code << bitstream->bits_left;
}

void gt_bitoutstream_append_bittab(GtBitOutStream *bitstream,
                                   GtBittab *tab) {
  unsigned long j,
                size = gt_bittab_size(tab);
  for (j = 0; j < size; j++) {
    if (bitstream->bits_left == 0) {
      gt_xfwrite(&bitstream->bitseqbuffer,
                 sizeof (GtBitsequence),
                 (size_t) 1, bitstream->fp);
      bitstream->bitseqbuffer = 0;
      bitstream->bits_left = GT_INTWORDSIZE;
      bitstream->written_bits += GT_INTWORDSIZE;
    }
    bitstream->bits_left--;
    if (gt_bittab_bit_is_set(tab, j))
      bitstream->bitseqbuffer |= ((GtBitsequence) 1) << (bitstream->bits_left);
  }
}

void gt_bitoutstream_flush(GtBitOutStream *bitstream)
{
  gt_assert(bitstream);
  gt_xfwrite(&bitstream->bitseqbuffer, sizeof (GtBitsequence),
             (size_t) 1, bitstream->fp);
  bitstream->written_bits += (GT_INTWORDSIZE - bitstream->bits_left);

  bitstream->bitseqbuffer = 0;
  bitstream->bits_left = GT_INTWORDSIZE;
}

void gt_bitoutstream_flush_advance(GtBitOutStream *bitstream)
{
  long fpos;
  bool is_not_at_pageborder = (ftell(bitstream->fp) % bitstream->pagesize) != 0;

  gt_assert(bitstream);

  gt_bitoutstream_flush(bitstream);

  if (is_not_at_pageborder) {
    fpos = (ftell(bitstream->fp) / bitstream->pagesize + 1) *
           bitstream->pagesize;
    gt_xfseek(bitstream->fp, fpos, SEEK_SET);
  }
}

long gt_bitoutstream_pos(const GtBitOutStream *bitstream)
{
  return ftell(bitstream->fp);
}

void gt_bitoutstream_delete(GtBitOutStream *bitstream)
{
  if (bitstream != NULL)
    gt_log_log("written %lu bits", bitstream->written_bits);
  gt_free(bitstream);
}
