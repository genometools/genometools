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

#include <inttypes.h>
#include "libgtcore/chardef.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"
#include "format64.h"
#include "iterseq.h"

struct Scansequenceiterator
{
  FastaBuffer *fb;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  Queue *descptr;
  Str *sequencebuffer;
  unsigned long length;
  uint64_t unitnum;
  bool withsequence, exhausted;
};

Scansequenceiterator *newScansequenceiterator(const StrArray *filenametab,
                                              const Uchar *symbolmap,
                                              bool withsequence,
                                              Env *env)
{
  Scansequenceiterator *sseqit;
  env_error_check(env);
  sseqit = env_ma_malloc(env, sizeof (Scansequenceiterator));
  sseqit->sequencebuffer = str_new(env);
  sseqit->length = 0;
  sseqit->descptr = queue_new(env);
  sseqit->fb = fastabuffer_new(filenametab,
                               symbolmap,
                               false,
                               NULL,
                               sseqit->descptr,
                               NULL,
                               env);
  sseqit->exhausted = false;
  sseqit->unitnum = 0;
  sseqit->withsequence = withsequence;
  return sseqit;
}

void freeScansequenceiterator(Scansequenceiterator **sseqit,Env *env)
{
  queue_delete_with_contents((*sseqit)->descptr,env);
  fastabuffer_delete((*sseqit)->fb, env);
  str_delete((*sseqit)->sequencebuffer, env);
  env_ma_free(*sseqit, env);
}

int nextScansequenceiterator(const Uchar **sequence,
                             unsigned long *len,
                             char **desc,
                             Scansequenceiterator *sseqit,
                             Env *env)
{
  Uchar charcode;
  int retval;
  bool haserr = false, foundseq = false;

  if (sseqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = fastabuffer_next(sseqit->fb,&charcode,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      sseqit->exhausted = true;
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sseqit->length == 0 && sseqit->withsequence)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(sseqit->unitnum));
        haserr = true;
        break;
      }
      *desc = queue_get(sseqit->descptr,env);
      *len = sseqit->length;
      if (sseqit->withsequence)
      {
        *sequence = str_get(sseqit->sequencebuffer); /* XXX: ownership prob. */
      }
      str_reset(sseqit->sequencebuffer);
      sseqit->length = 0;
      foundseq = true;
      sseqit->unitnum++;
      break;
    } else
    {
      if (sseqit->withsequence)
      {
        str_append_char(sseqit->sequencebuffer, charcode, env);
      }
      sseqit->length++;
    }
  }
  if (!haserr && sseqit->length > 0)
  {
    *desc = queue_get(sseqit->descptr,env);
    if (sseqit->withsequence)
    {
      *sequence = str_get(sseqit->sequencebuffer); /* XXX: ownership problem */
    }
    *len = sseqit->length;
    foundseq = true;
    str_reset(sseqit->sequencebuffer);
    sseqit->length = 0;
  }
  if (haserr)
  {
    return -1;
  }
  if (foundseq)
  {
    return 1;
  }
  return 0;
}
