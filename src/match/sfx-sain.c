/*
  Copyright (c) 2012-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/timer_api.h"
#include "core/mathsupport.h"
#include "sfx-lwcheck.h"
#include "bare-encseq.h"
#include "sfx-sain.h"
#include "sfx-linlcp.h"

#define GT_SAIN_SHOWTIMER(DESC)\
        if (timer != NULL && gt_logger_enabled(logger))\
        {\
          gt_timer_show_progress(timer,DESC,stdout);\
        }

/* GT_READMODE_FORWARD <=> GT_READMODE_REVERSE
   GT_READMODE_COMPL <=> GT_READMODE_REVCOMPL
*/

#define gt_readmode_inverse_direction(RM)\
        ((int) (RM)) % 2 == 0 ? (GtReadmode) ((int) (RM) + 1)\
                              : (GtReadmode) ((int) (RM) - 1)

typedef enum
{
  GT_SAIN_PLAINSEQ,
  GT_SAIN_INTSEQ,
  GT_SAIN_ENCSEQ,
  GT_SAIN_BARE_ENCSEQ
} GtSainSeqtype;

typedef signed int GtSsainindextype;

typedef struct
{
  GtUword totallength,
          numofchars;
  GtUsainindextype currentround,
                  *bucketsize,
                  *bucketfillptr,
                  *sstarfirstcharcount,
                  *roundtable;
  union
  {
    const GtEncseq *encseq;
    const GtUchar *plainseq;
    const GtUsainindextype *array;
  } seq;
  GtReadmode readmode; /* only relevant for encseq and bare_encseq */
  const GtBareEncseq *bare_encseq;
  GtSainSeqtype seqtype;
  bool bucketfillptrpoints2suftab,
       bucketsizepoints2suftab,
       roundtablepoints2suftab;
} GtSainseq;

static bool gt_sain_decideforfastmethod(GtUword maxvalue,
                                        GtUword len,
                                        GtUword numofchars)
{
  return maxvalue < (GtUword) GT_FIRSTTWOBITS && len > 1024UL
         && len >= GT_MULT2(numofchars)
         ? true : false;
}

static void gt_sain_allocate_tmpspace(GtSainseq *sainseq,
                                      GtUword maxvalue,
                                      GtUword len)
{
  sainseq->bucketsize
    = (GtUsainindextype *)
      gt_calloc((size_t) sainseq->numofchars,sizeof (*sainseq->bucketsize));
  sainseq->bucketfillptr
    = (GtUsainindextype *)
      gt_malloc(sizeof (*sainseq->bucketfillptr) * sainseq->numofchars);
  if (gt_sain_decideforfastmethod(maxvalue,len,sainseq->numofchars))
  {
    sainseq->roundtable = (GtUsainindextype *)
                          gt_malloc(sizeof (*sainseq->roundtable) *
                                    (size_t) GT_MULT2(sainseq->numofchars));
  } else
  {
    sainseq->roundtable = NULL;
  }
  sainseq->sstarfirstcharcount
    = (GtUsainindextype *)
      gt_calloc((size_t) sainseq->numofchars,
                sizeof (*sainseq->sstarfirstcharcount));
  sainseq->bucketfillptrpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  sainseq->roundtablepoints2suftab = false;
}

static GtSainseq *gt_sainseq_new_from_encseq(const GtEncseq *encseq,
                                             GtReadmode readmode)
{
  GtUword idx;
  GtSainseq *sainseq = (GtSainseq *) gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_ENCSEQ;
  sainseq->seq.encseq = encseq;
  sainseq->bare_encseq = NULL;
  sainseq->readmode = readmode;
  sainseq->totallength = gt_encseq_total_length(encseq);
  sainseq->numofchars = (GtUword) gt_encseq_alphabetnumofchars(encseq);

  gt_sain_allocate_tmpspace(sainseq,sainseq->totallength+GT_COMPAREOFFSET,
                                    sainseq->totallength);
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      sainseq->bucketsize[idx]
        = (GtUsainindextype)
          gt_encseq_charcount(encseq,(GtUchar) GT_COMPLEMENTBASE(idx));
    } else
    {
      sainseq->bucketsize[idx]
        = (GtUsainindextype) gt_encseq_charcount(encseq,(GtUchar) idx);
    }
  }
  return sainseq;
}

static GtSainseq *gt_sainseq_new_from_plainseq(const GtUchar *plainseq,
                                               GtUword len)
{
  const GtUchar *cptr;
  GtSainseq *sainseq = (GtSainseq *) gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_PLAINSEQ;
  sainseq->seq.plainseq = plainseq;
  sainseq->totallength = len;
  sainseq->numofchars = UCHAR_MAX+1;
  sainseq->bare_encseq = NULL;
  sainseq->readmode = GT_READMODE_FORWARD;
  gt_sain_allocate_tmpspace(sainseq,len+1,len);
  for (cptr = sainseq->seq.plainseq; cptr < sainseq->seq.plainseq + len; cptr++)
  {
    sainseq->bucketsize[*cptr]++;
  }
  return sainseq;
}

static GtSainseq *gt_sainseq_new_from_bare_encseq(const GtBareEncseq
                                                     *bare_encseq,
                                                  GtReadmode readmode)
{
  GtUword idx;
  GtSainseq *sainseq = (GtSainseq *) gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_BARE_ENCSEQ;
  sainseq->seq.plainseq = gt_bare_encseq_sequence(bare_encseq);
  sainseq->totallength = gt_bare_encseq_total_length(bare_encseq);
  sainseq->numofchars = gt_bare_encseq_numofchars(bare_encseq);
  sainseq->bare_encseq = bare_encseq;
  sainseq->readmode = readmode;
  gt_sain_allocate_tmpspace(sainseq,sainseq->totallength+GT_COMPAREOFFSET,
                            sainseq->totallength);
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      sainseq->bucketsize[idx]
        = (GtUsainindextype) gt_bare_encseq_charcount(bare_encseq,
                                             (GtUchar) GT_COMPLEMENTBASE(idx));
    } else
    {
      sainseq->bucketsize[idx]
        = (GtUsainindextype) gt_bare_encseq_charcount(bare_encseq,
                                                      (GtUchar) idx);
    }
  }
  return sainseq;
}

static GtSainseq *gt_sainseq_new_from_array(GtUsainindextype *arr,
                                            GtUword len,
                                            GtUword numofchars,
                                            GtUsainindextype *suftab,
                                            GtUsainindextype firstusable,
                                            GtUword suftabentries)
{
  GtUword charidx;
  GtUsainindextype *cptr;
  GtSainseq *sainseq = (GtSainseq *) gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_INTSEQ;
  sainseq->seq.array = arr;
  sainseq->totallength = len;
  sainseq->bare_encseq = NULL;
  sainseq->readmode = GT_READMODE_FORWARD;
  sainseq->numofchars = numofchars;
  gt_assert((GtUword) firstusable < suftabentries);
  if (suftabentries - firstusable >= numofchars)
  {
    /*printf("bucketsize: reclaim suftab["GT_WU","GT_WU"]\n",
            suftabentries - numofchars,suftabentries-1);*/
    sainseq->bucketsize = suftab + suftabentries - numofchars;
    sainseq->bucketsizepoints2suftab = true;
  } else
  {
    /*printf("bucketsize requires "GT_WU" entries and only "GT_WU" are left\n",
              numofchars,suftabentries - firstusable);*/
    sainseq->bucketsizepoints2suftab = false;
    sainseq->bucketsize
      = (GtUsainindextype *)
        gt_malloc(sizeof (*sainseq->bucketsize) * numofchars);
  }
  if (suftabentries - firstusable >= GT_MULT2(numofchars))
  {
    /*printf("bucketfillptr: reclaim suftab["GT_WU","GT_WU"]\n",
            suftabentries - GT_MULT2(numofchars),
            suftabentries - numofchars -1);*/
    sainseq->bucketfillptr = suftab + suftabentries - GT_MULT2(numofchars);
    sainseq->bucketfillptrpoints2suftab = true;
  } else
  {
    /*
    printf("bucketfillptr requires "GT_WU" entries and only "
           GT_WU" are left\n", numofchars,suftabentries - firstusable);*/
    sainseq->bucketfillptrpoints2suftab = false;
    sainseq->bucketfillptr
      = (GtUsainindextype *)
        gt_malloc(sizeof (*sainseq->bucketfillptr) * numofchars);
  }
  if (gt_sain_decideforfastmethod(len+1,len,numofchars))
  {
    if (suftabentries - firstusable >= GT_MULT4(numofchars))
    {
      /*printf("roundtable: reclaim suftab["GT_WU","GT_WU"]\n",
              suftabentries - GT_MULT4(numofchars),
              suftabentries - GT_MULT2(numofchars) - 1);*/
      sainseq->roundtable = suftab + suftabentries - GT_MULT4(numofchars);
      sainseq->roundtablepoints2suftab = true;
    } else
    {
      sainseq->roundtablepoints2suftab = false;
      sainseq->roundtable
        = (GtUsainindextype *)
          gt_malloc(sizeof (*sainseq->roundtable) * GT_MULT2(numofchars));
    }
  } else
  {
    sainseq->roundtablepoints2suftab = false;
    sainseq->roundtable = NULL;
  }
  sainseq->sstarfirstcharcount = NULL;
  for (charidx = 0; charidx < sainseq->numofchars; charidx++)
  {
    sainseq->bucketsize[charidx] = 0;
  }
  for (cptr = arr; cptr < arr + sainseq->totallength; cptr++)
  {
    gt_assert((GtUword) *cptr < numofchars);
    sainseq->bucketsize[*cptr]++;
  }
  return sainseq;
}

static void gt_sainseq_delete(GtSainseq *sainseq)
{
  if (sainseq != NULL)
  {
    if (!sainseq->bucketfillptrpoints2suftab)
    {
      gt_free(sainseq->bucketfillptr);
    }
    if (!sainseq->bucketsizepoints2suftab)
    {
      gt_free(sainseq->bucketsize);
    }
    if (!sainseq->roundtablepoints2suftab)
    {
      gt_free(sainseq->roundtable);
    }
    if (sainseq->seqtype != GT_SAIN_INTSEQ)
    {
      gt_free(sainseq->sstarfirstcharcount);
    }
    gt_free(sainseq);
  }
}

#ifdef GT_SAIN_WITHCOUNTS
static GtUword gt_sain_countcharaccess = 0;
#endif

static GtUword gt_sainseq_getchar(const GtSainseq *sainseq,
                                  GtUword position)
{
  gt_assert(position < sainseq->totallength);
#ifdef GT_SAIN_WITHCOUNTS
  gt_sain_countcharaccess++;
#endif
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      return (GtUword) sainseq->seq.plainseq[position];
    case GT_SAIN_INTSEQ:
      return (GtUword) sainseq->seq.array[position];
    case GT_SAIN_ENCSEQ:
      {
        GtUchar cc = gt_encseq_get_encoded_char(sainseq->seq.encseq,
                                                position,
                                                sainseq->readmode);
        return ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (GtUword) cc;
      }
    case GT_SAIN_BARE_ENCSEQ:
      {
        GtUchar cc = sainseq->seq.plainseq[position];
        return ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (GtUword) cc;
      }
  }
#ifndef S_SPLINT_S
  return 0;
#endif
}

static void gt_sain_endbuckets(GtSainseq *sainseq)
{
  GtUword charidx;
  GtUsainindextype previous;

  previous = sainseq->bucketfillptr[0] = sainseq->bucketsize[0];
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    previous += sainseq->bucketsize[charidx];
    sainseq->bucketfillptr[charidx] = previous;
  }
}

static void gt_sain_startbuckets(GtSainseq *sainseq)
{
  GtUword charidx;
  GtUsainindextype previous;

  previous = sainseq->bucketfillptr[0] = 0;
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    previous += sainseq->bucketsize[charidx-1];
    sainseq->bucketfillptr[charidx] = previous;
  }
}

typedef struct
{
  GtUword buf_size, numofchars, cachesize;
  GtUsainindextype *values, *fillptr, *suftab;
  uint16_t *nextidx;
  int log_bufsize;
  size_t size;
} GtSainbuffer;

static GtSainbuffer *gt_sainbuffer_new(GtUsainindextype *suftab,
                                       GtUsainindextype *fillptr,
                                       GtUword numofchars,
                                       GtUword totallength,
                                       GtLogger *logger)
{
  if (numofchars <= UCHAR_MAX+1)
  {
    GtSainbuffer *buf = (GtSainbuffer *) gt_malloc(sizeof (*buf));

    buf->size = sizeof (*buf);
    buf->fillptr = fillptr;
    buf->suftab = suftab;
    buf->log_bufsize = 18 - (sizeof (GtUsainindextype) == (size_t) 4 ? 1 : 2) -
                             (int) gt_determinebitspervalue(numofchars);
    buf->buf_size = 1UL << buf->log_bufsize;
    buf->numofchars = numofchars;
    gt_assert(buf->buf_size <= UINT16_MAX);
    buf->cachesize = numofchars << buf->log_bufsize;
    buf->size += sizeof (*buf->values) * buf->cachesize;
    buf->size += sizeof (*buf->nextidx) * numofchars;
    if ((GtUword) buf->size * 10UL >= totallength)
    {
      gt_free(buf);
      buf = NULL;
    } else
    {
      buf->values = (GtUsainindextype *)
                    gt_malloc(sizeof (*buf->values) * buf->cachesize);
      buf->nextidx = gt_calloc((size_t) numofchars,sizeof (*buf->nextidx));
      gt_logger_log(logger,"used buffer of "GT_WU" bytes (bufsize="GT_WU")",
                         (GtUword) buf->size,buf->buf_size);
    }
    return buf;
  } else
  {
    return NULL;
  }
}

static void gt_sainbuffer_update(GtSainbuffer *buf,
                                 GtUword charidx,
                                 GtUsainindextype value)
{
  const GtUword offset = charidx << buf->log_bufsize;

  buf->values[offset + (GtUword) buf->nextidx[charidx]] = value;
  if ((GtUword) buf->nextidx[charidx] < buf->buf_size - 1)
  {
    buf->nextidx[charidx]++;
  } else
  {
    GtUsainindextype *writeptr = buf->suftab + buf->fillptr[charidx] - 1,
                     *readptr = buf->values + offset;
    const GtUsainindextype *readend = readptr + buf->buf_size;

    while (readptr < readend)
    {
      *(writeptr--) = *(readptr++);
    }
    buf->nextidx[charidx] = 0;
    buf->fillptr[charidx] -= buf->buf_size;
  }
}

static void gt_sainbuffer_flushall(GtSainbuffer *buf)
{
  GtUword charidx;

  if (buf == NULL)
  {
    return;
  }
  for (charidx = 0; charidx < buf->numofchars; charidx++)
  {
    const GtUsainindextype bufleft = (GtUsainindextype) buf->nextidx[charidx];

    if (bufleft > 0)
    {
      GtUsainindextype *writeptr = buf->suftab + buf->fillptr[charidx] - 1,
                       *readptr = buf->values + (charidx << buf->log_bufsize);
      const GtUsainindextype *readend = readptr + bufleft;

      while (readptr < readend)
      {
        *(writeptr--) = *(readptr++);
      }
      buf->nextidx[charidx] = 0;
      buf->fillptr[charidx] -= bufleft;
    }
  }
}

static void gt_sainbuffer_delete(GtSainbuffer *buf)
{
  if (buf != NULL)
  {
    gt_free(buf->values);
    gt_free(buf->nextidx);
    gt_free(buf);
  }
}

#define GT_SAINUPDATEBUCKETPTR(CURRENTCC)\
        if (bucketptr != NULL)\
        {\
          if ((CURRENTCC) != lastupdatecc)\
          {\
            fillptr[lastupdatecc] = (GtUsainindextype) (bucketptr - suftab);\
            bucketptr = suftab + fillptr[lastupdatecc = CURRENTCC];\
          }\
        } else\
        {\
          bucketptr = suftab + fillptr[lastupdatecc = CURRENTCC];\
        }

static void gt_sain_special_singleSinduction1(GtSainseq *sainseq,
                                              GtSsainindextype *suftab,
                                              GtSsainindextype position)
{
  GtUword currentcc = gt_sainseq_getchar(sainseq,(GtUword) position);

  if (currentcc < sainseq->numofchars)
  {
    GtUword leftcontextcc;
    GtUsainindextype putidx = --sainseq->bucketfillptr[currentcc];

    gt_assert(position > 0);
    leftcontextcc = gt_sainseq_getchar(sainseq,(GtUword) --position);
    if (sainseq->roundtable != NULL)
    {
      GtUword t = (currentcc << 1) |
                  (leftcontextcc > currentcc ? 1UL : 0);

      gt_assert (sainseq->roundtable[t] <= sainseq->currentround);
      if (sainseq->roundtable[t] < sainseq->currentround)
      {
        sainseq->roundtable[t] = sainseq->currentround;
      }
      position += sainseq->totallength;
    }
    suftab[putidx] = (leftcontextcc > currentcc) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab["GT_WU"]="GT_WD"\n",putidx,suftab[putidx]);
#endif
  }
}

static void gt_sain_induceStypes1fromspecialranges(GtSainseq *sainseq,
                                                   GtSsainindextype *suftab)
{
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    if (gt_encseq_has_specialranges(sainseq->seq.encseq))
    {
      GtSpecialrangeiterator *sri;
      GtRange range;
      sri = gt_specialrangeiterator_new(sainseq->seq.encseq,
                                        GT_ISDIRREVERSE(sainseq->readmode)
                                          ? true : false);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        if (GT_ISDIRREVERSE(sainseq->readmode))
        {
          gt_range_reverse(sainseq->totallength,&range);
        }
        if (range.start > 1UL)
        {
          gt_sain_special_singleSinduction1(sainseq,
                                            suftab,
                                            (GtSsainindextype)
                                            (range.start - 1));
        }
      }
      gt_specialrangeiterator_delete(sri);
    }
  } else
  {
    gt_assert(sainseq->seqtype == GT_SAIN_BARE_ENCSEQ &&
              sainseq->bare_encseq != NULL);
    if (gt_bare_encseq_specialcharacters(sainseq->bare_encseq) > 0)
    {
      GtRange range;
      GtBareSpecialrangeiterator *sri
        = gt_bare_encseq_specialrangeiterator_new(sainseq->bare_encseq,
                                        GT_ISDIRREVERSE(sainseq->readmode)
                                          ? true : false);

      gt_assert(sri != NULL);
      while (gt_bare_encseq_specialrangeiterator_next(sri,&range))
      {
        if (GT_ISDIRREVERSE(sainseq->readmode))
        {
          gt_range_reverse(sainseq->totallength,&range);
        }
        if (range.start > 1UL)
        {
          gt_sain_special_singleSinduction1(sainseq,
                                            suftab,
                                            (GtSsainindextype)
                                            (range.start - 1));
        }
      }
      gt_bare_encseq_specialrangeiterator_delete(sri);
    }
  }
}

static void gt_sain_special_singleSinduction2(const GtSainseq *sainseq,
                                              GtSsainindextype *suftab,
                                              GtSsainindextype position,
                                              GT_UNUSED GtUword
                                                nonspecialentries)
{
  GtUword currentcc;

  currentcc = gt_sainseq_getchar(sainseq,(GtUword) --position);
  if (currentcc < sainseq->numofchars)
  {
    GtUsainindextype putidx = --sainseq->bucketfillptr[currentcc];

    gt_assert((GtUword) putidx < nonspecialentries);
    suftab[putidx] = (position == 0 ||
                      gt_sainseq_getchar(sainseq,
                                         (GtUword)
                                         (position-1)) > currentcc)
                      ? ~position : position;
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab["GT_WU"]="GT_WD"\n",putidx,suftab[putidx]);
#endif
  }
}

static void gt_sain_induceStypes2fromspecialranges(const GtSainseq *sainseq,
                                                   GtSsainindextype *suftab,
                                                   GtUword nonspecialentries)
{
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    if (gt_encseq_has_specialranges(sainseq->seq.encseq))
    {
      GtSpecialrangeiterator *sri;
      GtRange range;
      sri = gt_specialrangeiterator_new(sainseq->seq.encseq,
                                        GT_ISDIRREVERSE(sainseq->readmode)
                                          ? true : false);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        if (GT_ISDIRREVERSE(sainseq->readmode))
        {
          gt_range_reverse(sainseq->totallength,&range);
        }
        if (range.start > 0)
        {
          gt_sain_special_singleSinduction2(sainseq,
                                            suftab,
                                            (GtSsainindextype) range.start,
                                            nonspecialentries);
        }
      }
      gt_specialrangeiterator_delete(sri);
    }
  } else
  {
    gt_assert(sainseq->seqtype == GT_SAIN_BARE_ENCSEQ &&
              sainseq->bare_encseq != NULL);
    if (gt_bare_encseq_specialcharacters(sainseq->bare_encseq) > 0)
    {
      GtRange range;
      GtBareSpecialrangeiterator *sri
        = gt_bare_encseq_specialrangeiterator_new(sainseq->bare_encseq,
                                        GT_ISDIRREVERSE(sainseq->readmode)
                                          ? true : false);

      gt_assert(sri != NULL);
      while (gt_bare_encseq_specialrangeiterator_next(sri,&range))
      {
        if (GT_ISDIRREVERSE(sainseq->readmode))
        {
          gt_range_reverse(sainseq->totallength,&range);
        }
        if (range.start > 0)
        {
          gt_sain_special_singleSinduction2(sainseq,
                                            suftab,
                                            (GtSsainindextype) range.start,
                                            nonspecialentries);
        }
      }
      gt_bare_encseq_specialrangeiterator_delete(sri);
    }
  }
}

#include "match/sfx-sain.inc"

static GtUword gt_sain_insertSstarsuffixes(GtSainseq *sainseq,
                                           GtUsainindextype *suftab,
                                           GtLogger *logger)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      return gt_sain_PLAINSEQ_insertSstarsuffixes(sainseq,
                                                  sainseq->seq.plainseq,
                                                  suftab,logger);
    case GT_SAIN_ENCSEQ:
      return gt_sain_ENCSEQ_insertSstarsuffixes(sainseq,
                                                sainseq->seq.encseq,
                                                suftab,logger);
    case GT_SAIN_INTSEQ:
      return gt_sain_INTSEQ_insertSstarsuffixes(sainseq,
                                                sainseq->seq.array,
                                                suftab,logger);
    case GT_SAIN_BARE_ENCSEQ:
      return gt_sain_BARE_ENCSEQ_insertSstarsuffixes(sainseq,
                                                     sainseq->seq.plainseq,
                                                     suftab,logger);
  }
  /*@ignore@*/
  return 0;
  /*@end@*/
}

static void gt_sain_incrementfirstSstar(GtSainseq *sainseq,
                                        GtUsainindextype *suftab)
{
  GtUword charidx;
  GtUsainindextype sum = 0;

  for (charidx = 0; charidx < sainseq->numofchars; charidx++)
  {
    sum += sainseq->bucketsize[charidx];
    gt_assert(sainseq->bucketfillptr[charidx] <= sum);
    if (sainseq->bucketfillptr[charidx] < sum)
    {
      suftab[sainseq->bucketfillptr[charidx]] += sainseq->totallength;
    }
    sainseq->roundtable[charidx] = 0;
    sainseq->roundtable[charidx+sainseq->numofchars] = 0;
  }
}

static void gt_sain_induceLtypesuffixes1(GtSainseq *sainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
    case GT_SAIN_BARE_ENCSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_PLAINSEQ_induceLtypesuffixes1
        : gt_sain_PLAINSEQ_fast_induceLtypesuffixes1)
           (sainseq,sainseq->seq.plainseq,suftab,nonspecialentries);
      break;
    case GT_SAIN_ENCSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_ENCSEQ_induceLtypesuffixes1
        : gt_sain_ENCSEQ_fast_induceLtypesuffixes1)
           (sainseq,sainseq->seq.encseq,suftab,nonspecialentries);
      break;
    case GT_SAIN_INTSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_INTSEQ_induceLtypesuffixes1
        : gt_sain_INTSEQ_fast_induceLtypesuffixes1)
           (sainseq,sainseq->seq.array,suftab,nonspecialentries);
      break;
  }
}

static void gt_sain_adjustsuftab(GtUword totallength,GtSsainindextype *suftab,
                                 GtUword nonspecialentries)
{
  GtSsainindextype *suftabptr;

  for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
       suftabptr--)
  {
    if (*suftabptr > 0 && *suftabptr < (GtSsainindextype) totallength)
    {
      GtSsainindextype *nextgteq;

      *suftabptr += (GtSsainindextype) totallength;
      nextgteq = suftabptr - 1;
      while (nextgteq >= suftab && *nextgteq < (GtSsainindextype) totallength)
      {
        nextgteq--;
      }
      if (*nextgteq >= (GtSsainindextype) totallength)
      {
        *nextgteq -= (GtSsainindextype) totallength;
      }
      suftabptr = nextgteq;
    }
  }
}

static void gt_sain_induceStypesuffixes1(GtSainseq *sainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
    case GT_SAIN_BARE_ENCSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_PLAINSEQ_induceStypesuffixes1
        : gt_sain_PLAINSEQ_fast_induceStypesuffixes1)
           (sainseq,sainseq->seq.plainseq,suftab,nonspecialentries);
      break;
    case GT_SAIN_ENCSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_ENCSEQ_induceStypesuffixes1
        : gt_sain_ENCSEQ_fast_induceStypesuffixes1)
           (sainseq,sainseq->seq.encseq,suftab,nonspecialentries);
      break;
    case GT_SAIN_INTSEQ:
      (sainseq->roundtable == NULL
        ? gt_sain_INTSEQ_induceStypesuffixes1
        : gt_sain_INTSEQ_fast_induceStypesuffixes1)
           (sainseq,sainseq->seq.array,suftab,nonspecialentries);
      break;
  }
}

static void gt_sain_moveSstar2front(GtUword countSstartype,
                                    GtSsainindextype *suftab,
                                    GT_UNUSED GtUword nonspecialentries)
{
  GtUword readidx, writeidx = 0;
  GtSsainindextype position;

  for (readidx = 0; (position = suftab[readidx]) < 0; readidx++)
  {
    position = ~position;
    suftab[readidx] = position;
    gt_assert ((readidx + 1) < nonspecialentries);
  }
  if (readidx < countSstartype)
  {
    for (writeidx = readidx, readidx++; /* Nothing */; readidx++)
    {
      gt_assert (readidx < nonspecialentries);
      if ((position = suftab[readidx]) < 0)
      {
        position = ~position;
        gt_assert(writeidx < readidx);
        suftab[writeidx++] = position;
        suftab[readidx] = 0;
        if (writeidx == countSstartype)
        {
          break;
        }
      } else
      {
        suftab[readidx] = 0;
      }
    }
  } else
  {
    writeidx = readidx;
  }
  gt_assert(writeidx == countSstartype);
}

static GtUword gt_sain_fast_moveSstar2front(
                                      GtUword totallength,
                                      GtUword countSstartype,
                                      GtSsainindextype *suftab,
                                      GT_UNUSED GtUword nonspecialentries)
{
  GtUword readidx, namecount = 0, writeidx = 0;
  GtSsainindextype position;

  for (readidx = 0; (position = suftab[readidx]) < 0; readidx++)
  {
    position = ~position;
    if (position >= (GtSsainindextype) totallength)
    {
      namecount++;
    }
    suftab[readidx] = position;
    gt_assert ((readidx + 1) < nonspecialentries);
  }
  if (readidx < countSstartype)
  {
    for (writeidx = readidx, readidx++; /* Nothing */; readidx++)
    {
      gt_assert (readidx < nonspecialentries);
      if ((position = suftab[readidx]) < 0)
      {
        position = ~position;
        if (position >= (GtSsainindextype) totallength)
        {
          namecount++;
        }
        gt_assert(writeidx < readidx);
        suftab[writeidx++] = position;
        suftab[readidx] = 0;
        if (writeidx == countSstartype)
        {
          break;
        }
      } else
      {
        suftab[readidx] = 0;
      }
    }
  } else
  {
    writeidx = readidx;
  }
  gt_assert(writeidx == countSstartype);
  return namecount;
}

static void gt_sain_fast_assignSstarnames(GtUword totallength,
                                          GtUword countSstartype,
                                          GtUsainindextype *suftab,
                                          GtUword numberofnames,
                                          GtUword nonspecialentries)
{
  GtUsainindextype *suftabptr, *secondhalf = suftab + countSstartype;

  if ((GtUword) numberofnames < countSstartype)
  {
    GtUword currentname = numberofnames + 1;

    for (suftabptr = suftab + nonspecialentries - 1; suftabptr >= suftab;
         suftabptr--)
    {
      GtUsainindextype position;

      if ((position = *suftabptr) >= (GtUsainindextype) totallength)
      {
        position -= totallength;
        gt_assert(currentname > 0);
        currentname--;
      }
      if (currentname <= numberofnames)
      {
        secondhalf[GT_DIV2(position)] = (GtUsainindextype) currentname;
      }
    }
  } else
  {
    for (suftabptr = suftab; suftabptr < suftab + nonspecialentries;
         suftabptr++)
    {
      if (*suftabptr >= (GtUsainindextype) totallength)
      {
        *suftabptr -= totallength;
      }
    }
  }
}

static void gt_sain_induceLtypesuffixes2(const GtSainseq *sainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
    case GT_SAIN_BARE_ENCSEQ:
      gt_sain_PLAINSEQ_induceLtypesuffixes2(sainseq,sainseq->seq.plainseq,
                                            suftab,nonspecialentries);
      break;
    case GT_SAIN_ENCSEQ:
      gt_sain_ENCSEQ_induceLtypesuffixes2(sainseq,sainseq->seq.encseq,suftab,
                                          nonspecialentries);
      break;
    case GT_SAIN_INTSEQ:
      gt_sain_INTSEQ_induceLtypesuffixes2(sainseq,sainseq->seq.array,
                                          suftab,nonspecialentries);
      break;
  }
}

static void gt_sain_induceStypesuffixes2(const GtSainseq *sainseq,
                                         GtSsainindextype *suftab,
                                         GtUword nonspecialentries)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
    case GT_SAIN_BARE_ENCSEQ:
      gt_sain_PLAINSEQ_induceStypesuffixes2(sainseq,sainseq->seq.plainseq,
                                            suftab,nonspecialentries);
      break;
    case GT_SAIN_ENCSEQ:
      gt_sain_ENCSEQ_induceStypesuffixes2(sainseq,sainseq->seq.encseq,suftab,
                                          nonspecialentries);
      break;
    case GT_SAIN_INTSEQ:
      gt_sain_INTSEQ_induceStypesuffixes2(sainseq,sainseq->seq.array,
                                          suftab,nonspecialentries);
      break;
  }
}

static int gt_sain_compare_suffixes(const GtSainseq *sainseq,
                                    GtUword start1,
                                    GtUword start2)
{
  gt_assert(start1 <= sainseq->totallength &&
            start2 <= sainseq->totallength &&
            start1 != start2);

  while (true)
  {
    GtUword cc1, cc2;

    if (start1 == sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sainseq_getchar(sainseq,start1);
    cc2 = gt_sainseq_getchar(sainseq,start2);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    start1++;
    start2++;
  }
}

static void gt_sain_setundefined(bool fwd,GtUsainindextype *suftab,
                                 GtUword start, GtUword end)
{
  GtUsainindextype *ptr;

  if (fwd)
  {
    for (ptr = suftab + start; ptr <= suftab + end; ptr++)
    {
      *ptr = 0;
    }
  } else
  {
    for (ptr = suftab + end; ptr >= suftab + start; ptr--)
    {
      *ptr = 0;
    }
  }
}

static void gt_sain_movenames2front(GtUsainindextype *suftab,
                                    GtUword numberofsuffixes,
                                    GtUword totallength)
{
  GtUsainindextype *rptr, *wptr;
  const GtUsainindextype *maxrptr;

  gt_assert(totallength > 0);
  maxrptr = suftab + numberofsuffixes + GT_DIV2(totallength-1);
  for (rptr = wptr = suftab + numberofsuffixes; rptr <= maxrptr; rptr++)
  {
    if (*rptr > 0)
    {
      *(wptr++) = *rptr - 1UL; /* As we have used names with
                                  offset 1 to distinguish them
                                  from the undefined values
                                  signified by 0 */
    }
  }
  gt_assert(wptr == suftab + GT_MULT2(numberofsuffixes));
}

static void gt_sain_checkorder(const GtSainseq *sainseq,
                               const GtUsainindextype *suftab,
                               GtUword start,
                               GtUword end)
{
  GtUword idx;

  for (idx = start+1; idx <= end; idx++)
  {
    int cmp = gt_sain_compare_suffixes(sainseq,(GtUword) suftab[idx-1],
                                       (GtUword) suftab[idx]);

    if (cmp != -1)
    {
      fprintf(stderr,
              "%s: check interval ["GT_WU","GT_WU"] at idx="GT_WU": suffix "
              ""GT_WU" >= "GT_WU"\n",__func__,start,end,idx,
              (GtUword) suftab[idx-1],
              (GtUword) suftab[idx]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void gt_sain_expandorder2original(GtSainseq *sainseq,
                                         GtUword numberofsuffixes,
                                         GtUsainindextype *suftab)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      gt_sain_PLAINSEQ_expandorder2original(sainseq,sainseq->seq.plainseq,
                                            numberofsuffixes,suftab);
      break;
    case GT_SAIN_ENCSEQ:
      gt_sain_ENCSEQ_expandorder2original(sainseq,sainseq->seq.encseq,
                                          numberofsuffixes,suftab);
      break;
    case GT_SAIN_INTSEQ:
      gt_sain_INTSEQ_expandorder2original(sainseq,sainseq->seq.array,
                                          numberofsuffixes,suftab);
      break;
    case GT_SAIN_BARE_ENCSEQ:
      gt_sain_BARE_ENCSEQ_expandorder2original(sainseq,sainseq->seq.plainseq,
                                               numberofsuffixes,suftab);
      break;
  }
}

static void gt_sain_assignSstarlength(GtSainseq *sainseq,
                                      GtUsainindextype *lentab)
{
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      gt_sain_PLAINSEQ_assignSstarlength(sainseq,sainseq->seq.plainseq,lentab);
      break;
    case GT_SAIN_ENCSEQ:
      gt_sain_ENCSEQ_assignSstarlength(sainseq,sainseq->seq.encseq,lentab);
      break;
    case GT_SAIN_INTSEQ:
      gt_sain_INTSEQ_assignSstarlength(sainseq,sainseq->seq.array,lentab);
      break;
    case GT_SAIN_BARE_ENCSEQ:
      gt_sain_BARE_ENCSEQ_assignSstarlength(sainseq,sainseq->seq.plainseq,
                                            lentab);
      break;
  }
}

static GtUword gt_sain_assignSstarnames(const GtSainseq *sainseq,
                                        GtUword countSstartype,
                                        GtUsainindextype *suftab)
{
  GtUsainindextype *suftabptr, *secondhalf = suftab + countSstartype,
                   previouspos;
  GtUword previouslen, currentname = 1UL;

  previouspos = suftab[0];
  previouslen = (GtUword) secondhalf[GT_DIV2(previouspos)];
  secondhalf[GT_DIV2(previouspos)] = (GtUsainindextype) currentname;
  for (suftabptr = suftab + 1UL; suftabptr < suftab + countSstartype;
       suftabptr++)
  {
    int cmp = 0;
    GtUsainindextype position = *suftabptr;
    GtUword currentlen = 0;

    currentlen = (GtUword) secondhalf[GT_DIV2(position)];
    if (previouslen == currentlen)
    {
      switch (sainseq->seqtype)
      {
        case GT_SAIN_PLAINSEQ:
          cmp = gt_sain_PLAINSEQ_compare_Sstarstrings(sainseq,
                                                      sainseq->seq.plainseq,
                                                      (GtUword) previouspos,
                                                      (GtUword) position,
                                                      currentlen);
          break;
        case GT_SAIN_ENCSEQ:
          cmp = gt_sain_ENCSEQ_compare_Sstarstrings(sainseq,
                                                    sainseq->seq.encseq,
                                                    (GtUword) previouspos,
                                                    (GtUword) position,
                                                    currentlen);
          break;
        case GT_SAIN_INTSEQ:
          cmp = gt_sain_INTSEQ_compare_Sstarstrings(sainseq,
                                                    sainseq->seq.array,
                                                    (GtUword) previouspos,
                                                    (GtUword) position,
                                                    currentlen);
          break;
        case GT_SAIN_BARE_ENCSEQ:
          cmp = gt_sain_BARE_ENCSEQ_compare_Sstarstrings(sainseq,
                                                         sainseq->seq.plainseq,
                                                         (GtUword) previouspos,
                                                         (GtUword) position,
                                                         currentlen);
          break;
      }
      gt_assert(cmp != 1);
    } else
    {
      cmp = -1;
    }
    if (cmp == -1)
    {
      currentname++;
    }
    /* write the names in order of positions. As the positions of
       the Sstar suffixes differ by at least 2, the used address
       is unique */
    previouslen = currentlen;
    secondhalf[GT_DIV2(position)] = (GtUsainindextype) currentname;
    previouspos = position;
  }
  return currentname;
}

static void gt_sain_determineSstarfirstchardist(GtSainseq *sainseq)
{
  const GtUsainindextype *seqptr;
  GtUword nextcc = GT_UNIQUEINT(sainseq->totallength);
  bool nextisStype = true;

  gt_assert(sainseq->seqtype == GT_SAIN_INTSEQ);
  for (seqptr = sainseq->seq.array + sainseq->totallength - 1;
       seqptr >= sainseq->seq.array; seqptr--)
  {
    GtUword currentcc = (GtUword) *seqptr;
    bool currentisStype = (currentcc < nextcc ||
                           (currentcc == nextcc && nextisStype)) ? true : false;
    if (!currentisStype && nextisStype)
    {
      sainseq->sstarfirstcharcount[nextcc]++;
    }
    sainseq->bucketsize[currentcc]++;
    nextisStype = currentisStype;
    nextcc = currentcc;
  }
}

static void gt_sain_insertsortedSstarsuffixes(const GtSainseq *sainseq,
                                              GtUsainindextype *suftab,
                                              GtUword readidx,
                                              GtUword nonspecialentries)
{
  GtUword cc, fillidx = nonspecialentries;

  for (cc = sainseq->numofchars - 1; /* Nothing */; cc--)
  {
    if (sainseq->sstarfirstcharcount[cc] > 0)
    {
      GtUword putidx = fillidx - 1;

      gt_assert(readidx <= putidx);
      if (readidx < putidx)
      {
        GtUsainindextype offset;

        for (offset = 0; offset < sainseq->sstarfirstcharcount[cc]; offset++)
        {
          suftab[putidx - offset] = suftab[readidx - offset];
          suftab[readidx - offset] = 0;
#ifdef SAINSHOWSTATE
          printf("insertsorted: suftab["GT_WU"]="GT_WU"\n",putidx-offset,
                                                   suftab[putidx-offset]);
          printf("insertsorted: suftab["GT_WU"]=undef\n",readidx-offset);
#endif
        }
      }
    }
    gt_assert(fillidx >= (GtUword) sainseq->bucketsize[cc]);
    fillidx -= sainseq->bucketsize[cc];
    gt_assert(sainseq->bucketsize[cc] >=
              sainseq->sstarfirstcharcount[cc]);
    if (sainseq->bucketsize[cc] > sainseq->sstarfirstcharcount[cc])
    {
      gt_sain_setundefined(false,suftab,fillidx,
                           fillidx + sainseq->bucketsize[cc] -
                           sainseq->sstarfirstcharcount[cc] - 1);
    }
    readidx -= (GtUword) sainseq->sstarfirstcharcount[cc];
    if (cc == 0)
    {
      break;
    }
  }
}

static void gt_sain_filltailsuffixes(GtUsainindextype *suftabtail,
                                     const GtSainseq *sainseq,
                                     GtReadmode readmode)
{
  if (sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    GtUword specialcharacters
              = gt_encseq_specialcharacters(sainseq->seq.encseq),
            totallength = gt_encseq_total_length(sainseq->seq.encseq);

    if (gt_encseq_has_specialranges(sainseq->seq.encseq))
    {
      GtSpecialrangeiterator *sri;
      GtRange range;
      GtUword countspecial = 0;

      sri = gt_specialrangeiterator_new(sainseq->seq.encseq,
                                        GT_ISDIRREVERSE(readmode) ? false
                                                                  : true);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        GtUword idx;

        if (GT_ISDIRREVERSE(readmode))
        {
          gt_range_reverse(totallength,&range);
        }
        for (idx = range.start; idx < range.end; idx++)
        {
          gt_assert(countspecial < specialcharacters && idx < totallength);
          suftabtail[countspecial++] = (GtUsainindextype) idx;
        }
      }
      gt_assert(countspecial == specialcharacters);
      gt_specialrangeiterator_delete(sri);
    }
    suftabtail[specialcharacters] = (GtUsainindextype) totallength;
  } else
  {
    GtUword specialcharacters, totallength;

    gt_assert(sainseq->seqtype == GT_SAIN_BARE_ENCSEQ &&
              sainseq->bare_encseq != NULL);
    specialcharacters = gt_bare_encseq_specialcharacters(sainseq->bare_encseq);
    totallength = gt_bare_encseq_total_length(sainseq->bare_encseq);
    if (specialcharacters > 0)
    {
      GtRange range;
      GtUword countspecial = 0;
      GtBareSpecialrangeiterator *sri
        = gt_bare_encseq_specialrangeiterator_new(sainseq->bare_encseq,
                                        GT_ISDIRREVERSE(readmode) ? false
                                                                  : true);

      gt_assert(sri != NULL);
      while (gt_bare_encseq_specialrangeiterator_next(sri,&range))
      {
        GtUword idx;

        if (GT_ISDIRREVERSE(readmode))
        {
          gt_range_reverse(totallength,&range);
        }
        for (idx = range.start; idx < range.end; idx++)
        {
          gt_assert(countspecial < specialcharacters && idx < totallength);
          suftabtail[countspecial++] = (GtUsainindextype) idx;
        }
      }
      gt_assert(countspecial == specialcharacters);
      gt_bare_encseq_specialrangeiterator_delete(sri);
    }
    suftabtail[specialcharacters] = (GtUsainindextype) totallength;
  }
}

static GtUchar sain_accesschar_bare_encseq(const void *encseq,GtUword position,
                                           GT_UNUSED GtReadmode readmode)
{
  return gt_bare_encseq_get_encoded_char((const GtBareEncseq *) encseq,
                                         position);
}

GtUword sain_charcount_bare_encseq(const void *encseq,GtUchar idx)
{
  return gt_bare_encseq_charcount((const GtBareEncseq *) encseq,idx);
}

static GtUchar sain_accesschar_encseq(const void *encseq,GtUword position,
                                      GtReadmode readmode)
{
  return gt_encseq_get_encoded_char((const GtEncseq *) encseq,
                                    position,
                                    readmode);
}

GtUword sain_charcount_encseq(const void *encseq,GtUchar idx)
{
  return gt_encseq_charcount((const GtEncseq *) encseq, idx);
}

static GtUchar sain_accesschar_plainseq(const void *encseq,GtUword position,
                                        GT_UNUSED GtReadmode readmode)
{
  return ((const GtSainseq *) encseq)->seq.plainseq[position];
}

static GtUword sain_charcount_plainseq(const void *encseq,GtUchar idx)
{
  return (GtUword) ((const GtSainseq *) encseq)->bucketsize[idx];
}

static void gt_sain_rec_sortsuffixes(unsigned int level,
                                     GtSainseq *sainseq,
                                     GtUsainindextype *suftab,
                                     GtUsainindextype firstusable,
                                     GtUword nonspecialentries,
                                     GtUword suftabentries,
                                     bool intermediatecheck,
                                     bool finalcheck,
                                     GtLogger *logger,
                                     GtTimer *timer)
{
  GtUword countSstartype;

  GT_SAIN_SHOWTIMER("insert Sstar suffixes");
  countSstartype = gt_sain_insertSstarsuffixes(sainseq,suftab,logger);
  gt_logger_log(logger,"level %u: sort sequence of length "GT_WU" over "
                       ""GT_WU" symbols (%.2f)",
          level,sainseq->totallength,sainseq->numofchars,
          (double) sainseq->numofchars/sainseq->totallength);
  gt_logger_log(logger,"Sstar-type: "GT_WU" (%.2f)",countSstartype,
                 (double) countSstartype/sainseq->totallength);
  if (countSstartype > 0)
  {
    GtUword numberofnames;

    if (sainseq->roundtable != NULL)
    {
      gt_sain_incrementfirstSstar(sainseq,suftab);
    }
    gt_sain_startbuckets(sainseq);
    GT_SAIN_SHOWTIMER(sainseq->roundtable == NULL
                      ? "induce L suffixes" : "fast induce L suffixes");
    gt_sain_induceLtypesuffixes1(sainseq,(GtSsainindextype *) suftab,
                                 nonspecialentries);
    if (sainseq->roundtable != NULL)
    {
      gt_sain_adjustsuftab(sainseq->totallength,(GtSsainindextype *) suftab,
                           nonspecialentries);
    }
    gt_sain_endbuckets(sainseq);
    GT_SAIN_SHOWTIMER(sainseq->roundtable == NULL
                      ? "induce S suffixes" : "fast induce S suffixes");
    gt_sain_induceStypesuffixes1(sainseq,(GtSsainindextype *) suftab,
                                 nonspecialentries);
    if (sainseq->roundtable == NULL)
    {
      GT_SAIN_SHOWTIMER("simple moverSstar2front");
      gt_sain_moveSstar2front(countSstartype,(GtSsainindextype *) suftab,
                              nonspecialentries);
      GT_SAIN_SHOWTIMER("simple assignSstarlength");
      gt_sain_assignSstarlength(sainseq,suftab + countSstartype);
      GT_SAIN_SHOWTIMER("simple assignSstarnames");
      numberofnames = gt_sain_assignSstarnames(sainseq,countSstartype,suftab);
    } else
    {
      GT_SAIN_SHOWTIMER("fast moveSstar2front");
      numberofnames = gt_sain_fast_moveSstar2front(sainseq->totallength,
                                                   countSstartype,
                                                   (GtSsainindextype *) suftab,
                                                   nonspecialentries);
      if (!sainseq->roundtablepoints2suftab)
      {
        gt_free(sainseq->roundtable);
        sainseq->roundtable = NULL;
      }
      GT_SAIN_SHOWTIMER("fast assignSstarnames");
      gt_sain_fast_assignSstarnames(sainseq->totallength,countSstartype,
                                    suftab,numberofnames,nonspecialentries);
    }
    gt_assert(numberofnames <= countSstartype);
    if (numberofnames < countSstartype)
    {
    /* Now the name sequence is in the range from
       countSstartype .. 2 * countSstartype - 1 */
      GtUsainindextype *subseq = suftab + countSstartype;
      GtSainseq *sainseq_rec;

      GT_SAIN_SHOWTIMER("movenames2front");
      gt_sain_setundefined(true,suftab,0,countSstartype-1);
      gt_sain_movenames2front(suftab,countSstartype,sainseq->totallength);
      if (level == 0)
      {
        firstusable = (GtUsainindextype) GT_MULT2(countSstartype);
      }
      sainseq_rec = gt_sainseq_new_from_array(subseq,countSstartype,
                                              numberofnames,
                                              suftab,
                                              firstusable,
                                              suftabentries);
      gt_sain_rec_sortsuffixes(level+1,
                               sainseq_rec,
                               suftab,
                               firstusable,
                               countSstartype,
                               suftabentries,
                               intermediatecheck,
                               finalcheck,
                               logger,
                               timer);
      gt_sainseq_delete(sainseq_rec);
      GT_SAIN_SHOWTIMER("expandorder2original");
      gt_sain_expandorder2original(sainseq,countSstartype,suftab);
    } else
    {
      if (sainseq->seqtype == GT_SAIN_INTSEQ)
      {
        GtUword charidx;
        gt_assert(sainseq->sstarfirstcharcount == NULL);
        sainseq->sstarfirstcharcount = sainseq->bucketfillptr;
        for (charidx = 0; charidx < sainseq->numofchars; charidx++)
        {
          sainseq->sstarfirstcharcount[charidx] = 0;
          sainseq->bucketsize[charidx] = 0;
        }
        gt_sain_determineSstarfirstchardist(sainseq);
      }
    }
  }
  if (intermediatecheck && countSstartype > 0)
  {
    gt_sain_checkorder(sainseq,suftab,0,countSstartype-1);
  }
  GT_SAIN_SHOWTIMER("insert sorted Sstar suffixes");
  if (countSstartype > 0)
  {
    gt_sain_insertsortedSstarsuffixes (sainseq, suftab, countSstartype - 1,
                                       nonspecialentries);
  }
  gt_sain_startbuckets(sainseq);
  GT_SAIN_SHOWTIMER("final induce L suffixes");
  gt_sain_induceLtypesuffixes2(sainseq,(GtSsainindextype *) suftab,
                               nonspecialentries);
  gt_sain_endbuckets(sainseq);
  GT_SAIN_SHOWTIMER("final induce S suffixes");
  gt_sain_induceStypesuffixes2(sainseq,(GtSsainindextype *) suftab,
                               nonspecialentries);
  if (intermediatecheck && nonspecialentries > 0)
  {
    gt_sain_checkorder(sainseq,suftab,0,nonspecialentries-1);
  }
  if (sainseq->seqtype != GT_SAIN_INTSEQ)
  {
    if (sainseq->seqtype == GT_SAIN_BARE_ENCSEQ ||
        sainseq->seqtype == GT_SAIN_ENCSEQ)
    {
      GT_SAIN_SHOWTIMER("fill tail suffixes");
      gt_sain_filltailsuffixes(suftab + nonspecialentries,
                               sainseq,
                               sainseq->readmode);
    } else
    {
      suftab[sainseq->totallength] = (GtUsainindextype) sainseq->totallength;
    }
  }
  if (finalcheck && sainseq->seqtype != GT_SAIN_INTSEQ)
  {
    GT_SAIN_SHOWTIMER("check suffix order");
    switch (sainseq->seqtype)
    {
      case GT_SAIN_BARE_ENCSEQ:
        gt_suftab_lightweightcheck(sain_accesschar_bare_encseq,
                                   sain_charcount_bare_encseq,
                                   (void *) sainseq->bare_encseq,
                                   sainseq->readmode,
                                   sainseq->totallength,
                                   (unsigned int) sainseq->numofchars,
                                   suftab,
                                   sizeof *suftab,
                                   NULL);
        break;
      case GT_SAIN_ENCSEQ:
        gt_suftab_lightweightcheck(sain_accesschar_encseq,
                                   sain_charcount_encseq,
                                   (void *) sainseq->seq.encseq,
                                   sainseq->readmode,
                                   sainseq->totallength,
                                   (unsigned int) sainseq->numofchars,
                                   suftab,
                                   sizeof *suftab,
                                   NULL);
        break;
      case GT_SAIN_PLAINSEQ:
        gt_assert(sainseq->readmode == GT_READMODE_FORWARD);
        gt_suftab_lightweightcheck(sain_accesschar_plainseq,
                                   sain_charcount_plainseq,
                                   (void *) sainseq,
                                   sainseq->readmode,
                                   sainseq->totallength,
                                   (unsigned int) sainseq->numofchars,
                                   suftab,
                                   sizeof *suftab,
                                   NULL);
        break;
      default:
        gt_assert(false);
    }
  }
}

GtUsainindextype *gt_sain_encseq_sortsuffixes(const GtEncseq *encseq,
                                              GtReadmode readmode,
                                              bool intermediatecheck,
                                              bool finalcheck,
                                              GtLogger *logger,
                                              GtTimer *timer)
{
  GtUword nonspecialentries, suftabentries, totallength;
  GtUsainindextype *suftab;
  GtSainseq *sainseq;

  totallength = gt_encseq_total_length(encseq);
  nonspecialentries = totallength - gt_encseq_specialcharacters(encseq);
  suftabentries = totallength+1;
  suftab = (GtUsainindextype *) gt_calloc((size_t) suftabentries,
                                          sizeof *suftab);
  sainseq = gt_sainseq_new_from_encseq(encseq,readmode);
  gt_sain_rec_sortsuffixes(0,
                           sainseq,
                           suftab,
                           0,
                           nonspecialentries,
                           suftabentries,
                           intermediatecheck,
                           finalcheck,
                           logger,
                           timer);
#ifdef GT_SAIN_WITHCOUNTS
  printf("countcharaccess="GT_WU" (%.2f)\n",gt_sain_countcharaccess,
          (double) gt_sain_countcharaccess/sainseq->totallength);
#endif
  gt_sainseq_delete(sainseq);
  return suftab;
}

GtUsainindextype *gt_sain_bare_encseq_sortsuffixes(
                                      const GtBareEncseq *bare_encseq,
                                      GtReadmode readmode,
                                      bool intermediatecheck,
                                      bool finalcheck,
                                      GtLogger *logger,
                                      GtTimer *timer)
{
  GtUword nonspecialentries, suftabentries, totallength;
  GtUsainindextype *suftab;
  GtSainseq *sainseq;

  totallength = gt_bare_encseq_total_length(bare_encseq);
  nonspecialentries = totallength -
                      gt_bare_encseq_specialcharacters(bare_encseq);
  suftabentries = totallength+1;
  suftab = (GtUsainindextype *) gt_calloc((size_t) suftabentries,
                                          sizeof *suftab);
  sainseq = gt_sainseq_new_from_bare_encseq(bare_encseq,readmode);
  gt_sain_rec_sortsuffixes(0,
                           sainseq,
                           suftab,
                           0,
                           nonspecialentries,
                           suftabentries,
                           intermediatecheck,
                           finalcheck,
                           logger,
                           timer);
#ifdef GT_SAIN_WITHCOUNTS
  printf("countcharaccess="GT_WU" (%.2f)\n",gt_sain_countcharaccess,
          (double) gt_sain_countcharaccess/sainseq->totallength);
#endif
  gt_sainseq_delete(sainseq);
  return suftab;
}

GtUsainindextype *gt_sain_plain_sortsuffixes(const GtUchar *plainseq,
                                             GtUword len,
                                             bool intermediatecheck,
                                             bool finalcheck,
                                             GtLogger *logger,
                                             GtTimer *timer)
{
  GtUword suftabentries;
  GtUsainindextype *suftab;
  GtSainseq *sainseq;

  suftabentries = len+1;
  suftab = (GtUsainindextype *) gt_calloc((size_t) suftabentries,
                                          sizeof *suftab);
  sainseq = gt_sainseq_new_from_plainseq(plainseq,len);
  gt_sain_rec_sortsuffixes(0,
                           sainseq,
                           suftab,
                           0,
                           sainseq->totallength,
                           suftabentries,
                           intermediatecheck,
                           finalcheck,
                           logger,
                           timer);
#ifdef GT_SAIN_WITHCOUNTS
  printf("countcharaccess="GT_WU" (%.2f)\n",gt_sain_countcharaccess,
          (double) gt_sain_countcharaccess/sainseq->totallength);
#endif
  gt_sainseq_delete(sainseq);
  return suftab;
}

int gt_sain_checkmaxsequencelength(GtUword len,bool forencseq,GtError *err)
{
  GtUword maxsequencelength;

  if (forencseq)
  {
    maxsequencelength = (GtUword) (~GT_FIRSTBIT) - 1 - GT_COMPAREOFFSET;
  } else
  {
    maxsequencelength = (GtUword) (~GT_FIRSTBIT) - 1;
  }
  if (len > maxsequencelength)
  {
    gt_error_set(err,"sequence of size "GT_WU" is too long: sain algorithm "
                     "can only compute sequence of length up to "GT_WU"",
                     len,maxsequencelength);
    return -1;
  }
  return 0;
}

struct GtSainSufLcpIterator
{
  GtUsainindextype *suftab;
  unsigned int *plcptab;
  GtBareEncseq *bare_encseq;
  GtReadmode readmode;
  GtUword currentsuftabindex;
};

GtSainSufLcpIterator *gt_sain_suf_lcp_iterator_new(bool withlcp,
                                                   GtUchar *sequence,
                                                   GtUword len,
                                                   GtReadmode readmode,
                                                   GtUword numofchars,
                                                   GtError *err)
{
  GtSainSufLcpIterator *suflcpiterator;

  if (gt_sain_checkmaxsequencelength(len,false,err) != 0)
  {
    return NULL;
  }
  suflcpiterator = gt_malloc(sizeof *suflcpiterator);
  suflcpiterator->suftab = NULL;
  suflcpiterator->plcptab = NULL;
  suflcpiterator->currentsuftabindex = 0;
  suflcpiterator->readmode = readmode;
  suflcpiterator->bare_encseq = gt_bare_encseq_new(sequence,len,numofchars);
  gt_assert(suflcpiterator->bare_encseq != NULL);
  if (readmode != GT_READMODE_FORWARD)
  {
    bare_encseq_convert(suflcpiterator->bare_encseq,
                        GT_ISDIRREVERSE(readmode) ? false : true,
                        GT_ISDIRCOMPLEMENT(readmode) ? false : true);
  }
  suflcpiterator->suftab
    = gt_sain_bare_encseq_sortsuffixes(suflcpiterator->bare_encseq,
                                       readmode,
                                       false,
                                       false,
                                       NULL,
                                       NULL);
  if (withlcp)
  {
    GtUword maxlcp = 0,
            totallength
              = gt_bare_encseq_total_length(suflcpiterator->bare_encseq),
            partwidth
              = totallength -
                gt_bare_encseq_specialcharacters(suflcpiterator->bare_encseq);
    suflcpiterator->plcptab = gt_plain_lcp_phialgorithm(true,
                                                        &maxlcp,
                                                        sequence,
                                                        true,
                                                        partwidth,
                                                        totallength,
                                                        suflcpiterator->suftab);
  }
  return suflcpiterator;
}

GtReadmode gt_sain_suf_lcp_iterator_readmode(const GtSainSufLcpIterator
                                             *suflcpiterator)
{
  return suflcpiterator->readmode;
}

void gt_sain_suf_lcp_iterator_delete(GtSainSufLcpIterator *suflcpiterator)
{
  if (suflcpiterator != NULL)
  {
    gt_bare_encseq_delete(suflcpiterator->bare_encseq);
    gt_free(suflcpiterator->suftab);
    gt_free(suflcpiterator->plcptab);
    gt_free(suflcpiterator);
  }
}

GtUword gt_sain_suf_lcp_iterator_nonspecials(const GtSainSufLcpIterator
                                                   *suflcpiterator)
{
  GtUword totallength;

  gt_assert(suflcpiterator != NULL);
  totallength = gt_bare_encseq_total_length(suflcpiterator->bare_encseq);
  return totallength -
         gt_bare_encseq_specialcharacters(suflcpiterator->bare_encseq);
}

const GtBareEncseq *gt_sain_suf_lcp_iterator_bare_encseq(
         const GtSainSufLcpIterator *suflcpiterator)
{
  return suflcpiterator->bare_encseq;
}

GtUword gt_sain_suf_lcp_iterator_next(GtUword *lcpvalue,
                                      GtSainSufLcpIterator *suflcpiterator)
{
  GtUsainindextype previoussuffix, currentsuffix;
  gt_assert(suflcpiterator != NULL);

  previoussuffix = suflcpiterator->suftab[suflcpiterator->currentsuftabindex++];
  currentsuffix = suflcpiterator->suftab[suflcpiterator->currentsuftabindex];
  *lcpvalue = (GtUword) suflcpiterator->plcptab[currentsuffix];
  return (GtUword) previoussuffix;
}
