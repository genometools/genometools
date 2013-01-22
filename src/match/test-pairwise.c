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
#include "core/fa.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "test-pairwise.h"
#include "greedyedist.h"
#include "squarededist.h"

void gt_runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                                const char *file1,
                                const char *file2)
{
  const GtUchar *useq = NULL, *vseq = NULL;
  size_t ulen, vlen;
  bool forward = true;
  GtError *err;

  err = gt_error_new();
  useq = (const GtUchar *) gt_fa_mmap_read(file1,&ulen,err);
  if (useq == NULL)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(err));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  vseq = (const GtUchar *) gt_fa_mmap_read(file2,&vlen,err);
  if (vseq == NULL)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(err));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_error_delete(err);
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
  gt_fa_xmunmap((void *) useq);
  gt_fa_xmunmap((void *) vseq);
}

unsigned long gt_runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                                     const char *text)
{
  unsigned long i, len;

  len = (unsigned long) strlen(text);
  for (i=1UL; i<=len/2; i++)
  {
    checkfunction(true,
                  (const GtUchar *) text,
                  i,
                  (const GtUchar *) (text+i),
                  len-i);
  }
  return len/2;
}

unsigned long gt_applycheckfunctiontotext(const GtUchar *text,
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
                              unsigned long (*apply)(const GtUchar *,
                                                     unsigned long,
                                                     void *))
{
  unsigned long i, *w, z = textlen-1,
                testcases = 0,
                asize = (unsigned long) strlen(alpha);
  GtUchar *text;
  bool stop = false;

  w = gt_malloc(sizeof *w * (textlen+1));
  text = gt_malloc(sizeof *text * (textlen+1));
  for (i=0; i<=textlen; i++)
  {
    w[i] = 0;
  }
  text[textlen] = (GtUchar) '\0';
  while (!stop)
  {
    for (i = 0; i<textlen; i++)
    {
      text[i] = (GtUchar) alpha[w[i]];
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
  gt_free(w);
  gt_free(text);
  return testcases;
}

unsigned long gt_runcheckfunctiononalphalen(Checkcmppairfuntype checkfunction,
                                         const char *charlist,
                                         unsigned long len)
{
  return applyall(charlist,
                  len,
                  (void *) checkfunction,
                  gt_applycheckfunctiontotext);
}

void gt_checkgreedyunitedist(GT_UNUSED bool forward,
                             const GtUchar *useq,
                             unsigned long ulen,
                             const GtUchar *vseq,
                             unsigned long vlen)
{
  unsigned long edist1, edist2;
  GtFrontResource *frontresource = gt_frontresource_new(10UL);
  GtSeqabstract *greedyedistuseq = gt_seqabstract_new_ptr(useq,ulen,0),
                *greedyedistvseq = gt_seqabstract_new_ptr(vseq,vlen,0);

  edist1 = greedyunitedist(frontresource,greedyedistuseq,greedyedistvseq);
  edist2 = gt_squarededistunit (useq,ulen,vseq,vlen);
#ifdef SKDEBUG
  printf("edist = %lu\n",edist1);
#endif
  if (edist1 != edist2)
  {
    fprintf(stderr,"greedyunitedist = %lu != %lu = gt_squarededistunit\n",
                   edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_seqabstract_delete(greedyedistuseq);
  gt_seqabstract_delete(greedyedistvseq);
  gt_frontresource_delete(frontresource);
}
