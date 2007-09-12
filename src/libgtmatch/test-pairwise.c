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
#include "test-pairwise.h"

#include "greedyedist.pr"
#include "squarededist.pr"

void runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                                const char *file1,
                                const char *file2,
                                Env *env)
{
  const Uchar *useq = NULL, *vseq = NULL;
  size_t ulen, vlen;
  bool forward = true;

  useq = (const Uchar *) env_fa_mmap_read(env,file1,&ulen);
  if (useq == NULL)
  {
    fprintf(stderr,"cannot map file \"%s\": %s\n",file1,strerror(errno));
    exit(EXIT_FAILURE);
  }
  vseq = (const Uchar *) env_fa_mmap_read(env,file2,&vlen);
  if (vseq == NULL)
  {
    fprintf(stderr,"cannot map file \"%s\": %s",file2,strerror(errno));
    exit(EXIT_FAILURE);
  }
  while (true)
  {
    checkfunction(forward,useq,(unsigned long) ulen,
                          vseq,(unsigned long) vlen,env);
    if (!forward)
    {
      break;
    }
    forward = false;
  }
  env_fa_xmunmap((void *) useq,env);
  env_fa_xmunmap((void *) vseq,env);
}

unsigned long runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                                     const char *text,
                                     Env *env)
{
  unsigned long i, len;

  len = (unsigned long) strlen(text);
  for (i=1UL; i<=len/2; i++)
  {
    checkfunction(true,
                  (const Uchar *) text,
                  i,
                  (const Uchar *) (text+i),
                  len-i,
                  env);
  }
  return len/2;
}

unsigned long applycheckfunctiontotext(const Uchar *text,
                                       unsigned long textlen,
                                       void *info,
                                       Env *env)
{
  unsigned long i;
  Checkcmppairfuntype checkfunction = (Checkcmppairfuntype) info;

#ifdef DEBUG
  printf("%s\n",(char *) text);
#endif
  for (i=0; i<=textlen/2; i++)
  {
    checkfunction(true,text,i,text+i,textlen-i,env);
  }
  return textlen/2+1;
}

static unsigned long applyall(const char *alpha,
                              unsigned long textlen,void *info,
                              unsigned long (*apply)(const Uchar *,
                                                     unsigned long,
                                                     void *,Env *),
                              Env *env)
{
  unsigned long i, *w, z = textlen-1,
                testcases = 0,
                asize = (unsigned long) strlen(alpha);
  Uchar *text;
  bool stop = false;

  ALLOCASSIGNSPACE(w,NULL,unsigned long,textlen+1);
  ALLOCASSIGNSPACE(text,NULL,Uchar,textlen+1);
  for (i=0; i<=textlen; i++)
  {
    w[i] = 0;
  }
  text[textlen] = (Uchar) '\0';
  while (!stop)
  {
    for (i = 0; i<textlen; i++)
    {
      text[i] = (Uchar) alpha[w[i]];
    }
    testcases += apply(text,textlen,info,env);
    while (true)
    {
      w[z]++;
      if (w[z] == asize)
      {
        w[z] = 0;
        if (z == 0)
        {
          stop = true;
          break;
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
  return testcases;
}

unsigned long runcheckfunctiononalphalen(Checkcmppairfuntype checkfunction,
                                         const char *charlist,
                                         unsigned long len,
                                         Env *env)
{
  return applyall(charlist,
                  len,
                  (void *) checkfunction,
                  applycheckfunctiontotext,
                  env);
}

void checkgreedyunitedist(/*@unused@*/ bool forward,
                          const Uchar *useq,
                          unsigned long ulen,
                          const Uchar *vseq,
                          unsigned long vlen,
                          Env *env)
{
  unsigned long edist1, edist2;

  edist1 = greedyunitedist(useq,ulen,vseq,vlen,env);
  edist2 = squarededistunit (useq,ulen,vseq,vlen,env);
#ifdef DEBUG
  printf("edist = %lu\n",edist1);
#endif
  if (edist1 != edist2)
  {
    fprintf(stderr,"greedyunitedist = %lu != %lu = squarededistunit\n",
                   edist1,edist2);
    exit(EXIT_FAILURE);
  }
}
