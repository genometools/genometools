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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include "core/cstr_api.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_rep.h"
#include "core/sequence_buffer_inline.h"

#define FASTASEPARATOR    '>'
#define NEWLINESYMBOL     '\n'
#define CRSYMBOL          '\r'

struct GtSequenceBufferFasta {
  const GtSequenceBuffer parent_instance;
  GtStr *headerbuffer;
  bool indesc,
       firstseqinfile,
       firstoverallseq,
       nextfile;
};

#define gt_sequence_buffer_fasta_cast(SB)\
        gt_sequence_buffer_cast(gt_sequence_buffer_fasta_class(), SB)

static int gt_sequence_buffer_fasta_advance(GtSequenceBuffer *sb, GtError *err)
{
  int currentchar, ret = 0;
  unsigned long currentoutpos = 0, currentfileadd = 0, currentfileread = 0;
  GtSequenceBufferMembers *pvt;
  GtSequenceBufferFasta *sbf;

  gt_error_check(err);

  sbf = (GtSequenceBufferFasta*) sb;
  pvt = sb->pvt;
  while (true)
  {
    if (currentoutpos >= (unsigned long) OUTBUFSIZE)
    {
      if (pvt->filelengthtab != NULL)
      {
        pvt->filelengthtab[pvt->filenum].length
          += (uint64_t) currentfileread;
        pvt->filelengthtab[pvt->filenum].effectivelength
          += (uint64_t) currentfileadd;
      }
      break;
    }
    if (sbf->nextfile)
    {
      if (pvt->filelengthtab != NULL)
      {
        pvt->filelengthtab[pvt->filenum].length = 0;
        pvt->filelengthtab[pvt->filenum].effectivelength = 0;
      }
      sbf->nextfile = false;
      sbf->indesc = false;
      sbf->firstseqinfile = true;
      currentfileadd = 0;
      currentfileread = 0;
      pvt->linenum = (uint64_t) 1;
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
          pvt->filelengthtab[pvt->filenum].length += currentfileread;
          pvt->filelengthtab[pvt->filenum].effectivelength += currentfileadd;
        }
        if ((unsigned long) pvt->filenum
                                     == gt_str_array_size(pvt->filenametab)-1)
        {
          pvt->complete = true;
          break;
        }
        pvt->filenum++;
        sbf->nextfile = true;
      } else
      {
        currentfileread++;
        if (sbf->indesc)
        {
          if (currentchar == NEWLINESYMBOL)
          {
            pvt->linenum++;
            sbf->indesc = false;
          }
          if (pvt->descptr != NULL)
          {
            if (currentchar == NEWLINESYMBOL)
            {
              gt_desc_buffer_finish(pvt->descptr);
            } else
            {
              if (currentchar != CRSYMBOL)
                gt_desc_buffer_append_char(pvt->descptr, currentchar);
            }
          }
        } else
        {
          if (!isspace((int) currentchar))
          {
            if (currentchar == FASTASEPARATOR)
            {
              if (sbf->firstoverallseq)
              {
                sbf->firstoverallseq = false;
                sbf->firstseqinfile = false;
              } else
              {
                if (sbf->firstseqinfile)
                {
                  sbf->firstseqinfile = false;
                } else
                {
                  currentfileadd++;
                }
                pvt->outbuf[currentoutpos++] = (unsigned char) SEPARATOR;
                pvt->lastspeciallength++;
              }
              sbf->indesc = true;
            } else
            {
              if ((ret = process_char(sb, currentoutpos,
                                      (unsigned char) currentchar, err)))
                return ret;
              currentoutpos++;
              currentfileadd++;
            }
          }
        }
      }
    }
  }
  if (sbf->firstoverallseq)
  {
    gt_error_set(err,"no sequences in multiple fasta file(s) %s ...",
              gt_str_array_get(pvt->filenametab,0));
    return -2;
  }
  pvt->nextfree = currentoutpos;
  return 0;
}

static void gt_sequence_buffer_fasta_free(GtSequenceBuffer *sb)
{
  GtSequenceBufferFasta *sbf = gt_sequence_buffer_fasta_cast(sb);
  gt_file_delete(sb->pvt->inputstream);
  gt_str_delete(sbf->headerbuffer);
}

static unsigned long
gt_sequence_buffer_fasta_get_file_index(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  return (unsigned long) sb->pvt->filenum;
}

bool gt_sequence_buffer_fasta_guess(const char* txt)
{
  return (*txt == FASTASEPARATOR);
}

const GtSequenceBufferClass* gt_sequence_buffer_fasta_class(void)
{
  static const GtSequenceBufferClass sbc = { sizeof (GtSequenceBufferFasta),
                                        gt_sequence_buffer_fasta_advance,
                                        gt_sequence_buffer_fasta_get_file_index,
                                        gt_sequence_buffer_fasta_free };
  return &sbc;
}

GtSequenceBuffer* gt_sequence_buffer_fasta_new(const GtStrArray *sequences)
{
  GtSequenceBuffer *sb;
  GtSequenceBufferFasta *sbf;
  sb = gt_sequence_buffer_create(gt_sequence_buffer_fasta_class());
  sbf = gt_sequence_buffer_fasta_cast(sb);
  sb->pvt->filenametab = sequences;
  sbf->headerbuffer = gt_str_new();
  sb->pvt->filenum = 0;
  sbf->firstoverallseq = true;
  sbf->firstseqinfile = true;
  sbf->nextfile = true;
  sb->pvt->nextread = sb->pvt->nextfree = 0;
  sb->pvt->complete = false;
  sb->pvt->lastspeciallength = 0;
  return sb;
}
