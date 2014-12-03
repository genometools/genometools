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

#include <ctype.h>
#include <limits.h>
#include <string.h>

#include "core/alphabet_api.h"
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/readmode.h"
#include "core/types_api.h"
#include "match/bare-encseq.h"

typedef struct
{
  GtUword start,
          length;
} GtBareSpecialrange;

GT_DECLAREARRAYSTRUCT(GtBareSpecialrange);

struct GtBareEncseq
{
  GtUchar *sequence;
  GtUword totallength, specialcharacters, numofchars, *charcount;
  GtArrayGtBareSpecialrange specialranges;
};

void gt_bare_encseq_delete(GtBareEncseq *bare_encseq)
{
  if (bare_encseq != NULL)
  {
    gt_free(bare_encseq->charcount);
    GT_FREEARRAY(&bare_encseq->specialranges,GtBareSpecialrange);
    gt_free(bare_encseq);
  }
}

GtBareEncseq *gt_bare_encseq_new(GtUchar *sequence,GtUword len,
                                 GtUword numofchars)
{
  GtBareEncseq *bare_encseq = gt_malloc(sizeof *bare_encseq);
  const GtUchar *readptr;
  GtBareSpecialrange *srptr = NULL;
  GtUword lastspecialrange_length = 0;

  bare_encseq->specialcharacters = 0;
  bare_encseq->numofchars = numofchars;
  bare_encseq->charcount = gt_calloc((size_t) bare_encseq->numofchars,
                                     sizeof *bare_encseq->charcount);
  GT_INITARRAY(&bare_encseq->specialranges,GtBareSpecialrange);
  for (readptr = sequence; readptr < sequence + len; readptr++)
  {
    GtUchar cc = *readptr;
    if (ISSPECIAL(cc))
    {
      if (lastspecialrange_length == 0)
      {
        GT_GETNEXTFREEINARRAY(srptr,&bare_encseq->specialranges,
                              GtBareSpecialrange,128UL);
        srptr->start = (GtUword) (readptr - sequence);
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
  }
  if (lastspecialrange_length > 0)
  {
    gt_assert(srptr != NULL);
    srptr->length = lastspecialrange_length;
  }
  bare_encseq->sequence = sequence;
  bare_encseq->totallength = len;
  return bare_encseq;
}

GtBareEncseq *gt_bare_encseq_parse_new(GtUchar *filecontents,size_t numofbytes,
                                       const GtAlphabet *alphabet,
                                       GtError *err)
{
  GtUchar *writeptr = filecontents, *readptr = filecontents;
  const GtUchar *endptr = filecontents + numofbytes;
  bool firstline = true, haserr = false;
  GtUword lastspecialrange_length = 0;
  GtBareSpecialrange *srptr = NULL;
  GtBareEncseq *bare_encseq = gt_malloc(sizeof *bare_encseq);
  const GtUchar *smap = gt_alphabet_symbolmap(alphabet);

  bare_encseq->specialcharacters = 0;
  bare_encseq->numofchars = (GtUword) gt_alphabet_num_of_chars(alphabet);
  bare_encseq->charcount = gt_calloc((size_t) bare_encseq->numofchars,
                                     sizeof *bare_encseq->charcount);
  GT_INITARRAY(&bare_encseq->specialranges,GtBareSpecialrange);
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
                                GtBareSpecialrange,128UL);
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
          if (cc == UNDEFCHAR)
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
                                    GtBareSpecialrange,128UL);
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
  if (lastspecialrange_length > 0)
  {
    gt_assert(srptr != NULL);
    srptr->length = lastspecialrange_length;
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

const GtUchar *gt_bare_encseq_sequence(const GtBareEncseq *bare_encseq)
{
  gt_assert(bare_encseq != NULL);
  return bare_encseq->sequence;
}

GtUchar gt_bare_encseq_get_encoded_char(const GtBareEncseq *bare_encseq,
                                        GtUword position)
{
  return bare_encseq->sequence[position];
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

struct GtBareSpecialrangeiterator
{
  const GtBareSpecialrange *startptr, *endptr, *current;
  bool moveforward;
};

GtBareSpecialrangeiterator* gt_bare_encseq_specialrangeiterator_new(
                       const GtBareEncseq *bare_encseq,
                       bool moveforward)
{
  gt_assert(bare_encseq != NULL);
  if (bare_encseq->specialranges.nextfreeGtBareSpecialrange == 0)
  {
    return NULL;
  } else
  {
    GtBareSpecialrangeiterator *sri = gt_malloc(sizeof *sri);

    sri->startptr = bare_encseq->specialranges.spaceGtBareSpecialrange;
    sri->endptr = bare_encseq->specialranges.spaceGtBareSpecialrange +
                  bare_encseq->specialranges.nextfreeGtBareSpecialrange;
    sri->moveforward = moveforward;
    sri->current = moveforward ? sri->startptr : sri->endptr - 1;
    return sri;
  }
}

void gt_bare_encseq_specialrangeiterator_delete(GtBareSpecialrangeiterator *sri)
{
  if (sri != NULL)
  {
    gt_free(sri);
  }
}

bool gt_bare_encseq_specialrangeiterator_next(GtBareSpecialrangeiterator *sri,
                                              GtRange *range)
{
  if (sri->current != NULL)
  {
    range->start = sri->current->start;
    range->end = sri->current->start + sri->current->length;
    if (sri->moveforward)
    {
      if (sri->current < sri->endptr - 1)
      {
        sri->current++;
      } else
      {
        sri->current = NULL;
      }
    } else
    {
      if (sri->current > sri->startptr)
      {
        sri->current--;
      } else
      {
        sri->current = NULL;
      }
    }
    return true;
  }
  return false;
}

void bare_encseq_convert(GtBareEncseq *bare_encseq,bool forward,bool direct)
{
  GtUchar *leftptr, *rightptr;

  if (forward)
  {
    gt_assert(!direct);
    for (leftptr = bare_encseq->sequence;
         leftptr < bare_encseq->sequence + bare_encseq->totallength;
         leftptr++)
    {
      if (ISNOTSPECIAL(*leftptr))
      {
        *leftptr = GT_COMPLEMENTBASE(*leftptr);
      }
    }
  } else
  {
    if (direct)
    {
      for (leftptr = bare_encseq->sequence,
           rightptr = bare_encseq->sequence + bare_encseq->totallength - 1;
           leftptr < rightptr; leftptr++, rightptr--)
      {
        GtUchar tmp = *leftptr;
        *leftptr = *rightptr;
        *rightptr = tmp;
      }
    } else
    {
      for (leftptr = bare_encseq->sequence,
           rightptr = bare_encseq->sequence + bare_encseq->totallength - 1;
           leftptr <= rightptr; leftptr++, rightptr--)
      {
        GtUchar tmp = *leftptr;
        *leftptr = ISSPECIAL(*rightptr) ? *rightptr
                                        : GT_COMPLEMENTBASE(*rightptr);
        *rightptr = ISSPECIAL(tmp) ? tmp
                                   : GT_COMPLEMENTBASE(tmp);
      }
    }
  }
}
