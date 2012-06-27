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
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/safearith.h"
#include "core/xansi_api.h"
#include "extended/bitinstream.h"

struct GtBitInStream {
  bool last_chunk;
  char *path;
  GtBitsequence *bitseqbuffer;
  int cur_bit;
  unsigned long cur_bitseq,
                bufferlength,
                pages_to_map,
                read_bits;
  size_t cur_filepos,
         filesize;
  long pagesize;
};

GtBitInStream *gt_bitinstream_new(const char* path,
                                  size_t offset,
                                  unsigned long pages_to_map)
{
  GtBitInStream *bitstream = gt_malloc(sizeof (*bitstream));

  bitstream->pagesize =sysconf((int) _SC_PAGESIZE);
  bitstream->last_chunk = false;
  gt_safe_assign(bitstream->filesize, gt_file_estimate_size(path));
  bitstream->path = gt_cstr_dup(path);
  bitstream->pages_to_map = pages_to_map;

  if ((unsigned long) bitstream->filesize <
      (bitstream->pages_to_map * bitstream->pagesize))
    bitstream->pages_to_map =
      (unsigned long) ((bitstream->filesize / bitstream->pagesize) + 1);

  bitstream->bitseqbuffer = NULL;
  bitstream->read_bits = 0;
  gt_bitinstream_reinit(bitstream,
                        offset);

  bitstream->bufferlength = (bitstream->pages_to_map * bitstream->pagesize) /
                            sizeof (*bitstream->bitseqbuffer);
  return bitstream;
}

void gt_bitinstream_reinit(GtBitInStream *bitstream,
                           size_t offset)
{
  size_t mapsize = (size_t) (bitstream->pagesize * bitstream->pages_to_map);

  gt_assert(offset < bitstream->filesize);
  gt_assert((offset % bitstream->pagesize) == 0);

  bitstream->cur_filepos = offset;

  gt_fa_xmunmap(bitstream->bitseqbuffer);

  if (bitstream->cur_filepos + mapsize > bitstream->filesize) {
    mapsize = bitstream->filesize - bitstream->cur_filepos;
    bitstream->bufferlength = (unsigned long)  mapsize /
                                sizeof (*bitstream->bitseqbuffer);
    bitstream->last_chunk = true;
  }
  bitstream->bitseqbuffer =
    gt_fa_xmmap_read_range(bitstream->path,
                           mapsize,
                           offset);

  gt_assert(bitstream->bitseqbuffer != NULL);

  bitstream->cur_bit = 0;
  bitstream->cur_bitseq = 0;
}

int gt_bitinstream_get_next_bit(GtBitInStream *bitstream,
                                bool * bit)
{
  const int eof = 0, more_to_read = 1;
  if (bitstream->cur_bit == GT_INTWORDSIZE) {
    if (bitstream->cur_bitseq < bitstream->bufferlength - 1) {
      bitstream->cur_bit = 0;
      bitstream->cur_bitseq++;
    }
    else {
      if (bitstream->last_chunk) {
        return eof;
      }
      else {
        gt_bitinstream_reinit(bitstream,
                              bitstream->cur_filepos +
                                bitstream->pagesize *
                                bitstream->pages_to_map);
      }
    }
  }
  gt_assert(bitstream->cur_bitseq != bitstream->bufferlength);
  *bit =  GT_ISBITSET(bitstream->bitseqbuffer[bitstream->cur_bitseq],
                      bitstream->cur_bit++) != 0;
  bitstream->read_bits++;
  return more_to_read;
}

void gt_bitinstream_delete(GtBitInStream *bitstream)
{
  if (bitstream != NULL) {
    gt_log_log("read %lu bits", bitstream->read_bits);
    gt_fa_xmunmap(bitstream->bitseqbuffer);
    gt_free(bitstream->path);
    gt_free(bitstream);
  }
}
