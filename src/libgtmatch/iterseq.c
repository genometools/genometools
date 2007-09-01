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
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "arraydef.h"
#include "chardef.h"
#include "fbs-def.h"
#include "seqdesc.h"
#include "gqueue-def.h"
#include "format64.h"

#include "fbsadv.pr"
#include "genericqueue.pr"

#include "readnextUchar.gen"

int overallquerysequences(int(*processsequence)(void *,
                                                uint64_t,
                                                const Uchar *,
                                                unsigned long,
                                                const char *,
                                                Env *),
                          void *info,
                          ArrayUchar *sequencebuffer,
                          const StrArray *filenametab,
                          Sequencedescription *sequencedescription,
                          const Uchar *symbolmap,
                          Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  int retval;
  uint64_t unitnum = 0;
  bool haserr = false;
  char *desc;

  env_error_check(env);
  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        sequencedescription,
                        env);
  sequencebuffer->nextfreeUchar = 0;
  while (true)
  {
    retval = readnextUchar(&charcode,&fbs,env);
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
      desc = dequeuegeneric(sequencedescription->descptr,env);
      if (desc == NULL)
      {
        haserr = true;
        break;
      }
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
    desc = dequeuegeneric(sequencedescription->descptr,env);
    if (desc == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      if (processsequence(info,
                         unitnum,
                         sequencebuffer->spaceUchar,
                         sequencebuffer->nextfreeUchar,
                         desc,
                         env) != 0)
      {
        haserr = true;
      }
    }
    FREESPACE(desc);
    sequencebuffer->nextfreeUchar = 0;
  }
  return haserr ? -1 : 0;
}
