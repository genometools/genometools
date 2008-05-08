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
#include "libgtcore/fa.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "spacedef.h"
#include "test-pairwise.h"

#include "greedyedist.pr"
#include "squarededist.pr"

void runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                                const char *file1,
                                const char *file2)
{
  const Uchar *useq = NULL, *vseq = NULL;
  size_t ulen, vlen;
  bool forward = true;

  useq = (const Uchar *) fa_mmap_read(file1,&ulen);
  if (useq == NULL)
  {
    fprintf(stderr,"cannot map file \"%s\": %s\n",file1,strerror(errno));
    exit(EXIT_FAILURE);
  }
  vseq = (const Uchar *) fa_mmap_read(file2,&vlen);
  if (vseq == NULL)
  {
    fprintf(stderr,"cannot map file \"%s\": %s",file2,strerror(errno));
    exit(EXIT_FAILURE);
  }
  while (true)
  {
    checkfunction(forward,useq,(unsigned long) ulen,
                          vseq,(unsigned long) vlen);
    if (!forward)
    {
      break;
    }
    forward = false;
  }
  fa_xmunmap((void *) useq);
  fa_xmunmap((void *) vseq);
}

unsigned long runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                                     const char *text)
{
  unsigned long i, len;

  len = (unsigned long) strlen(text);
  for (i=1UL; i<=len/2; i++)
  {
    checkfunction(true,
                  (const Uchar *) text,
                  i,
                  (const Uchar *) (text+i),
                  len-i);
  }
  return len/2;
}

unsigned long applycheckfunctiontotext(const Uchar *text,
                                       unsigned long textlen,
                                       void *info)
{
  unsigned long i;
  Checkcmppairfuntype checkfunction = (Checkcmppairfuntype) info;

#ifdef SKDEBUG
  printf("%s\n",(char *) text);
#endif
  for (i=0; i<=textlen/2; i++)
  {
    checkfunction(true,text,i,text+i,textlen-i);
  }
  return textlen/2+1;
}

static unsigned long applyall(const char *alpha,
                              unsigned long textlen,void *info,
                              unsigned long (*apply)(const Uchar *,
                                                     unsigned long,
                                                     void *))
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
    testcases += apply(text,textlen,info);
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
                                         unsigned long len)
{
  return applyall(charlist,
                  len,
                  (void *) checkfunction,
                  applycheckfunctiontotext);
}

void checkgreedyunitedist(UNUSED bool forward,
                          const Uchar *useq,
                          unsigned long ulen,
                          const Uchar *vseq,
                          unsigned long vlen)
{
  unsigned long edist1, edist2;

  edist1 = greedyunitedist(useq,ulen,vseq,vlen);
  edist2 = squarededistunit (useq,ulen,vseq,vlen);
#ifdef SKDEBUG
  printf("edist = %lu\n",edist1);
#endif
  if (edist1 != edist2)
  {
    fprintf(stderr,"greedyunitedist = %lu != %lu = squarededistunit\n",
                   edist1,edist2);
    exit(EXIT_FAILURE);
  }
}
