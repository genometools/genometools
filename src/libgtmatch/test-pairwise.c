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

#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "symboldef.h"
#include "spacedef.h"

int runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                               const char *file1,
                               const char *file2,
                               Env *env)
{
  const Uchar *useq = NULL, *vseq = NULL;
  size_t ulen, vlen;
  bool haserr = false;

  useq = (const Uchar *) env_fa_mmap_read(env,file1,&ulen);
  if (useq == NULL)
  {
    env_error_set(env,"cannot map file \"%s\": %s",file1,strerror(errno));
    haserr = true;
  }
  if (!haserr)
  {
    vseq = (const Uchar *) env_fa_mmap_read(env,file2,&vlen);
    if (vseq == NULL)
    {
      env_error_set(env,"cannot map file \"%s\": %s",file2,strerror(errno));
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (checkfunction(true,useq,(unsigned long) ulen,
                           vseq,(unsigned long) vlen) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (checkfunction(false,useq,(unsigned long) ulen,
                            vseq,(unsigned long) vlen) != 0)
    {
      haserr = true;
    }
  }
  env_fa_xmunmap((void *) useq,env);
  env_fa_xmunmap((void *) vseq,env);
  return haserr ? -1 : 0;
}

int runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                           const char *text)
{
  size_t i, len;

  len = strlen(text);
  for (i=(size_t) 1; i<=len/2; i++)
  {
    if (checkfunction(true,(const Uchar *) text,(unsigned long) i,
                           (const Uchar *) (text+i), (unsigned long) len-i)
        != 0)
    {
      return -1;
    }
  }
  return 0;
}

int applycheckfunctiontotext(const Uchar *text,
                             unsigned long textlen,
                             void *info)
{
  unsigned long i;
  Checkcmppairfuntype checkfunction = (Checkcmppairfuntype) info;

  printf("%s\n",(char *) text);
  for (i=0; i<=textlen/2; i++)
  {
    if (checkfunction(true,text,i,text+i,textlen-i) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static int applyall(const char *alpha,unsigned long textlen,void *info,
                    int (*apply)(const Uchar *,unsigned long,void *),
                    Env *env)
{
  unsigned long i, *w, z = textlen-1, asize = (unsigned long) strlen(alpha);
  Uchar *text;
  bool haserr = false;

  ALLOCASSIGNSPACE(w,NULL,unsigned long,textlen+1);
  ALLOCASSIGNSPACE(text,NULL,Uchar,textlen+1);
  for (i=0; i<=textlen; i++)
  {
    w[i] = 0;
  }
  text[textlen] = (Uchar) '\0';

  while (true)
  {
    for (i = 0; i<textlen; i++)
    {
      text[i] = (Uchar) alpha[w[i]];
    }
    if (apply(text,textlen,info) != 0)
    {
      haserr = true;
      break;
    }
    while (true)
    {
      w[z]++;
      if (w[z] == asize)
      {
        w[z] = 0;
        if (z == 0)
        {
          return 0;
        }
        z--;
      } else
      {
        z = textlen-1;
        /*@innerbreak@*/ break;
      }
    }
  }
  FREESPACE(w);
  FREESPACE(text);
  /*@ignore@*/
  return haserr ? -1 : 0;
  /*@end@*/
}

int runcheckfunctiononalphalen(Checkcmppairfuntype checkfunction,
                               const char *charlist,
                               unsigned long len,
                               Env *env)
{
  if (applyall(charlist,
               len,
               (void *) checkfunction,
               applycheckfunctiontotext,
               env) != 0)
  {
    return -1;
  }
  return 0;
}
