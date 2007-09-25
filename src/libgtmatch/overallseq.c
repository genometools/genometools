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
#include "arraydef.h"
#include "format64.h"
#include "overallseq.h"

int overallquerysequences(int(*processsequence)(void *,
                                                uint64_t,
                                                const Uchar *,
                                                unsigned long,
                                                const char *,
                                                Env *),
                          void *info,
                          ArrayUchar *sequencebuffer,
                          const StrArray *filenametab,
                          Queue *descptr,
                          const Uchar *symbolmap,
                          Env *env)
{
  FastaBuffer *fb;
  Uchar charcode;
  int retval;
  uint64_t unitnum = 0;
  bool haserr = false;
  char *desc;

  env_error_check(env);
  fb = fastabuffer_new(filenametab,
                       symbolmap,
                       false,
                       NULL,
                       descptr,
                       NULL,
                       env);
  sequencebuffer->nextfreeUchar = 0;
  while (true)
  {
    retval = fastabuffer_next(fb,&charcode,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sequencebuffer->nextfreeUchar == 0)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(unitnum));
        haserr = true;
        break;
      }
      desc = queue_get(descptr,env);
      if (processsequence(info,
                          unitnum,
                          sequencebuffer->spaceUchar,
                          sequencebuffer->nextfreeUchar,
                          desc,
                          env) != 0)
      {
        haserr = true;
        FREESPACE(desc);
        break;
      }
      FREESPACE(desc);
      sequencebuffer->nextfreeUchar = 0;
      unitnum++;
    } else
    {
      STOREINARRAY(sequencebuffer,Uchar,1024,charcode);
    }
  }
  if (!haserr && sequencebuffer->nextfreeUchar > 0)
  {
    desc = queue_get(descptr,env);
    if (processsequence(info,
                        unitnum,
                        sequencebuffer->spaceUchar,
                        sequencebuffer->nextfreeUchar,
                        desc,
                        env) != 0)
    {
      haserr = true;
    }
    FREESPACE(desc);
    sequencebuffer->nextfreeUchar = 0;
  }
  fastabuffer_delete(fb, env);
  return haserr ? -1 : 0;
}
