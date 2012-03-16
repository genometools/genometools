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
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/xansi_api.h"
#include "extended/bitinstream.h"

struct GtBitInStream {
  unsigned long cur_bit,
                cur_bitseq,
                bufferlength;
  off_t cur_page;
  GtBitsequence *bitseqbuffer;
  off_t filesize;
  char *path;
  long pagesize;
  unsigned long num_of_pages;
};

GtBitInStream *gt_bitinstream_new(char* path,
                                  off_t offset,
                                  unsigned long pages_to_map,
                                  GtError *err)
{
  int had_err = 0;
  GtBitInStream *bitstream = gt_malloc(sizeof (*bitstream));

  bitstream->pagesize =sysconf(_SC_PAGESIZE);
  bitstream->filesize = gt_file_estimate_size(path);
  bitstream->path = strdup(path);
  bitstream->num_of_pages = pages_to_map;

  if (bitstream->filesize < bitstream->num_of_pages * bitstream->pagesize)
    bitstream->num_of_pages = (bitstream->filesize / bitstream->pagesize) + 1;

  had_err = gt_bitinstream_reinit(bitstream,
                                  offset,
                                  err);

  bitstream->bufferlength = (bitstream->num_of_pages * bitstream->pagesize) /
                            sizeof (*bitstream->bitseqbuffer);
  if (!had_err)
    return bitstream;
  else
    return NULL;
}

int gt_bitinstream_reinit(GtBitInStream *bitstream,
                              off_t offset,
                              GT_UNUSED GtError *err)
{
  bitstream->cur_page = offset;
  gt_assert(bitstream->cur_page < bitstream->filesize);
  gt_assert((bitstream->cur_page % bitstream->pagesize) == 0);

  gt_fa_xmunmap(bitstream->bitseqbuffer);

  bitstream->bitseqbuffer = gt_fa_xmmap_read_range(bitstream->path,
                                                   bitstream->pagesize *
                                                     bitstream->num_of_pages,
                                                   bitstream->cur_page);

  bitstream->cur_bit = 0;
  bitstream->cur_bitseq = 0;

  return 0;
}

int gt_bitinstream_get_next_bit(GtBitInStream *bitstream,
                                    bool * bit,
                                    GtError *err)
{
  int had_err = 0;
  if (bitstream->cur_bit == GT_INTWORDSIZE) {
    if (bitstream->cur_bitseq < bitstream->bufferlength - 1) {
      bitstream->cur_bit = 0;
      bitstream->cur_bitseq++;
    }
    else {
      if (bitstream->filesize <=
            bitstream->cur_page + (bitstream->num_of_pages *
              bitstream->pagesize))
        return 0;
      else
        had_err = gt_bitinstream_reinit(bitstream,
                                        bitstream->cur_page +
                                          bitstream->pagesize *
                                          bitstream->num_of_pages,
                                        err);
    }
    if (had_err)
      return had_err;
  }
  *bit =  GT_ISBITSET(bitstream->bitseqbuffer[bitstream->cur_bitseq],
                      bitstream->cur_bit++);
  return 1;
}

void gt_bitinstream_delete(GtBitInStream *bitstream)
{
  if (!bitstream)
    return;
  gt_fa_xmunmap(bitstream->bitseqbuffer);
  gt_free(bitstream);
}
