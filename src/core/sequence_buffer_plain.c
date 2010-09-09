/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#include "core/error.h"
#include "core/sequence_buffer_plain.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_inline.h"
#include "core/unused_api.h"

struct GtSequenceBufferPlain {
  const GtSequenceBuffer parent_instance;
  bool nextfile,
       firstseqinfile;
};

#define gt_sequence_buffer_plain_cast(SB)\
        gt_sequence_buffer_cast(gt_sequence_buffer_plain_class(), SB)

static int gt_sequence_buffer_plain_advance(GtSequenceBuffer *sb, GtError *err)
{
  int currentchar;
  unsigned long currentoutpos = 0, currentfileread = 0;
  GtSequenceBufferMembers *pvt;
  GtSequenceBufferPlain *sbp;

  sbp = gt_sequence_buffer_plain_cast(sb);
  pvt = sb->pvt;

  gt_error_check(err);
  if (pvt->descptr != NULL)
  {
    gt_error_set(err, "no headers in plain sequence file");
    return -1;
  }
  while (true)
  {
    if (currentoutpos >= (unsigned long) OUTBUFSIZE)
    {
      if (pvt->filelengthtab != NULL)
      {
        pvt->filelengthtab[pvt->filenum].length
           += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength
           += (uint64_t) currentfileread;
      }
      break;
    }
    if (sbp->nextfile)
    {
      if (pvt->filelengthtab != NULL)
      {
        pvt->filelengthtab[pvt->filenum].length = 0;
        pvt->filelengthtab[pvt->filenum].effectivelength = 0;
      }
      sbp->nextfile = false;
      sbp->firstseqinfile = true;
      currentfileread = 0;
      pvt->inputstream = gt_file_xopen(gt_str_array_get(pvt->filenametab,
                                                  (unsigned long) pvt->filenum),
                                          "rb");
      pvt->currentinpos = 0;
      pvt->currentfillpos = 0;
    } else
    {
      currentchar = inlinebuf_getchar(sb, pvt->inputstream);
      if (currentchar == EOF)
      {
        gt_file_delete(pvt->inputstream);
        pvt->inputstream = NULL;
        if (pvt->filelengthtab != NULL)
        {
          pvt->filelengthtab[pvt->filenum].length
            += (uint64_t) currentfileread;
          pvt->filelengthtab[pvt->filenum].effectivelength
            += (uint64_t) currentfileread;
        }
        if ((unsigned long) pvt->filenum
                                       == gt_str_array_size(pvt->filenametab)-1)
        {
          pvt->complete = true;
          break;
        }
        pvt->filenum++;
        sbp->nextfile = true;
      } else
      {
        currentfileread++;
        pvt->outbuf[currentoutpos++] = (unsigned char) currentchar;
      }
    }
  }
  if (currentoutpos == 0)
  {
    gt_error_set(err, "no characters in plain file(s) %s ...",
              gt_str_array_get(pvt->filenametab,0));
    return -2;
  }
  pvt->nextfree = currentoutpos;
  return 0;
}

static unsigned long
gt_sequence_buffer_plain_get_file_index(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  return (unsigned long) sb->pvt->filenum;
}

void gt_sequence_buffer_plain_free(GT_UNUSED GtSequenceBuffer *sb)
{
  /* not needed */
}

const GtSequenceBufferClass* gt_sequence_buffer_plain_class(void)
{
  static const GtSequenceBufferClass sbc = { sizeof (GtSequenceBufferPlain),
                                        gt_sequence_buffer_plain_advance,
                                        gt_sequence_buffer_plain_get_file_index,
                                        gt_sequence_buffer_plain_free };
  return &sbc;
}

GtSequenceBuffer* gt_sequence_buffer_plain_new(const GtStrArray *sequences)
{
  GtSequenceBuffer *sb;
  GtSequenceBufferPlain *sbf;
  sb = gt_sequence_buffer_create(gt_sequence_buffer_plain_class());
  sbf = gt_sequence_buffer_plain_cast(sb);
  sb->pvt->filenametab = sequences;
  sb->pvt->filenum = 0;
  sbf->firstseqinfile = true;
  sbf->nextfile = true;
  sb->pvt->nextread = sb->pvt->nextfree = 0;
  sb->pvt->complete = false;
  sb->pvt->lastspeciallength = 0;
  return sb;
}
