/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include <string.h>
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/arraydef.h"
#include "sfx-sain-bareesq.h"

struct GtBareEncseq
{
  GtUchar *sequence;
  GtUword totallength, specialcharacters, numofchars, *charcount;
  GtArrayGtSainSpecialrange specialranges;
};

void gt_bare_encseq_delete(GtBareEncseq *bare_encseq)
{
  if (bare_encseq != NULL)
  {
    gt_free(bare_encseq->charcount);
    printf("number of special characters: " GT_WU " (" GT_WU " range(s))\n",
           bare_encseq->specialcharacters,
           bare_encseq->specialranges.nextfreeGtSainSpecialrange);
    GT_FREEARRAY(&bare_encseq->specialranges,GtSainSpecialrange);
    gt_free(bare_encseq);
  }
}

GtBareEncseq *gt_bare_encseq_new(GtUchar *filecontents,size_t numofbytes,
                                 GtError *err)
{
  GtUchar *writeptr = filecontents, *readptr = filecontents, smap[UCHAR_MAX+1];
  const GtUchar undefined = (GtUchar) UCHAR_MAX,
        *endptr = filecontents + numofbytes;
  bool firstline = true, haserr = false;
  size_t idx;
  const char *wildcard_list = "nsywrkvbdhmNSYWRKVBDHM";
  GtUword lastspecialrange_length = 0;
  GtSainSpecialrange *srptr = NULL;
  GtBareEncseq *bare_encseq = gt_malloc(sizeof *bare_encseq);

  for (idx = 0; idx <= (size_t) UCHAR_MAX; idx++)
  {
    smap[idx] = undefined;
  }
  smap['a'] = 0;
  smap['A'] = 0;
  smap['c'] = (GtUchar) 1;
  smap['C'] = (GtUchar) 1;
  smap['g'] = (GtUchar) 2;
  smap['G'] = (GtUchar) 2;
  smap['t'] = (GtUchar) 3;
  smap['T'] = (GtUchar) 3;
  smap['u'] = (GtUchar) 3;
  smap['U'] = (GtUchar) 3;
  bare_encseq->specialcharacters = 0;
  bare_encseq->numofchars = 4UL;
  bare_encseq->charcount = gt_calloc((size_t) bare_encseq->numofchars,
                                     sizeof *bare_encseq->charcount);
  GT_INITARRAY(&bare_encseq->specialranges,GtSainSpecialrange);
  for (idx = 0; idx < strlen(wildcard_list); idx++)
  {
    smap[(int) wildcard_list[idx]] = WILDCARD;
  }
  readptr = filecontents;
  while (!haserr && readptr < endptr)
  {
    if (*readptr == '>')
    {
      if (!firstline)
      {
        if (lastspecialrange_length == 0)
        {
          GT_GETNEXTFREEINARRAY(srptr,&bare_encseq->specialranges,
                                GtSainSpecialrange,128UL);
          srptr->start = (GtUword) (writeptr - filecontents);
        }
        lastspecialrange_length++;
        *writeptr++ = SEPARATOR;
        bare_encseq->specialcharacters++;
      } else
      {
        firstline = false;
      }
      while (readptr < endptr && *readptr != '\n')
      {
        readptr++;
      }
      readptr++;
    } else
    {
      while (readptr < endptr && *readptr != '\n')
      {
        if (!isspace(*readptr))
        {
          GtUchar cc = smap[*readptr];
          if (cc == undefined)
          {
            gt_error_set(err,"illegal input characters %c\n",*readptr);
            haserr = true;
            break;
          }
          if (ISSPECIAL(cc))
          {
            if (lastspecialrange_length == 0)
            {
              GT_GETNEXTFREEINARRAY(srptr,&bare_encseq->specialranges,
                                    GtSainSpecialrange,128UL);
              srptr->start = (GtUword) (writeptr - filecontents);
            }
            lastspecialrange_length++;
            bare_encseq->specialcharacters++;
          } else
          {
            gt_assert((GtUword) cc < bare_encseq->numofchars);
            bare_encseq->charcount[(int) cc]++;
            if (lastspecialrange_length > 0)
            {
              gt_assert(srptr != NULL);
              srptr->length = lastspecialrange_length;
            }
            lastspecialrange_length = 0;
          }
          *writeptr++ = cc;
        }
        readptr++;
      }
      readptr++;
    }
  }
  bare_encseq->sequence = filecontents;
  bare_encseq->totallength = (GtUword) (writeptr - filecontents);
  if (haserr)
  {
    gt_bare_encseq_delete(bare_encseq);
    return NULL;
  }
  return bare_encseq;
}

const GtArrayGtSainSpecialrange *gt_bare_encseq_specialranges(
                                           const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return &bare_encseq->specialranges;
}

const GtUchar *gt_bare_encseq_sequence(const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return bare_encseq->sequence;
}

GtUword gt_bare_encseq_total_length(const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return bare_encseq->totallength;
}

GtUword gt_bare_encseq_numofchars(const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return bare_encseq->numofchars;
}

GtUword gt_bare_encseq_charcount(const GtBareEncseq *bare_encseq,GtUchar idx)
{
  gt_assert(bare_encseq != NULL && (GtUword) idx < bare_encseq->numofchars);
  return bare_encseq->charcount[idx];
}

GtUword gt_bare_encseq_specialcharacters(const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return bare_encseq->specialcharacters;
}
