/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef FASTABUFFER_IMP_H
#define FASTABUFFER_IMP_H

#include <sys/types.h>

#define OUTPUTFILEBUFFERSIZE 4096
#define INPUTFILEBUFFERSIZE  4096

struct FastaBuffer
{
  unsigned int filenum;
  uint64_t linenum;
  unsigned long nextread,
                nextfree;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete,
       nextfile;
  Queue *descptr;
  GenFile *inputstream;
  Uchar outputbuffer[OUTPUTFILEBUFFERSIZE],
        inputbuffer[INPUTFILEBUFFERSIZE];
  ssize_t currentinpos, currentfillpos;
  uint64_t lastspeciallength;
  Filelengthvalues *filelengthtab;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  bool plainformat;
  unsigned long *characterdistribution;
  Arraychar headerbuffer;
};

int advanceformatbufferstate(FastaBuffer *fb, Error *);

static inline int fastabuffer_next(FastaBuffer *fb,Uchar *val, Error *e)
{
  if (fb->nextread >= fb->nextfree)
  {
    if (fb->complete)
    {
      return 0;
    }
    if (advanceformatbufferstate(fb, e) != 0)
    {
      return -1;
    }
    fb->nextread = 0;
    if (fb->nextfree == 0)
    {
      return 0;
    }
  }
  *val = fb->outputbuffer[fb->nextread++];
  return 1;
}

#endif
