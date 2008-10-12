/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/fa.h"
#include "core/unused_api.h"
#include "core/seqiterator.h"
#include "core/chardef.h"
#include "revcompl.h"
#include "alphadef.h"
#include "format64.h"
#include "encseq-def.h"
#include "divmodmul.h"
#include "opensfxfile.h"
#include "tyr-search.h"
#include "spacedef.h"

struct Tyrcountinfo
{
  void *mappedmctfileptr;
  const GtStr *indexfilename;
  Uchar *smallcounts;
  Largecount *largecounts;
  unsigned long numoflargecounts;
};

struct Tyrindex
{
  void *mappedfileptr;
  const GtStr *indexfilename;
  unsigned int alphasize;
  size_t numofmers;
  unsigned long mersize,
                merbytes;
  Uchar *mertable,
        *lastmer;
};

typedef struct
{
  Uchar *leftmer,
        *rightmer;
} Merbounds;

unsigned long decodesingleinteger(const Uchar *start)
{
  unsigned long idx, value;

  value = (unsigned long) start[0];
  for (idx=1UL; idx < (unsigned long) sizeof (unsigned long); idx++)
  {
    value |= (((unsigned long) start[idx]) << MULT8(idx));
  }
  return value;
}

Tyrindex *tyrindex_new(const GtStr *tyrindexname,GtError *err)
{
  bool haserr = false;
  size_t numofbytes, rest;
  Tyrindex *tyrindex;

  ALLOCASSIGNSPACE(tyrindex,NULL,Tyrindex,1);
  tyrindex->indexfilename = tyrindexname;
  tyrindex->mappedfileptr = genericmaponlytable(tyrindexname,
                                                MERSUFFIX,&numofbytes,err);
  if (tyrindex->mappedfileptr == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    tyrindex->mertable = (Uchar *) tyrindex->mappedfileptr;
    rest = sizeof (unsigned long) * EXTRAINTEGERS;
    if (rest > numofbytes)
    {
      gt_error_set(err,"index must contain at least %lu bytes",
                        (unsigned long) rest);
      haserr = true;
    }
  }
  if (!haserr)
  {
    assert(tyrindex->mertable != NULL);
    tyrindex->mersize
      = decodesingleinteger(tyrindex->mertable + numofbytes - rest);
    tyrindex->alphasize
      = (unsigned int) decodesingleinteger(tyrindex->mertable +
                                           numofbytes - rest +
                                           sizeof (unsigned long));
    tyrindex->merbytes = MERBYTES(tyrindex->mersize);
    if ((numofbytes - rest) % tyrindex->merbytes != 0)
    {
      gt_error_set(err,"size of index is %lu which is not a multiple of %lu",
                   (unsigned long) (numofbytes - rest),
                   tyrindex->merbytes);
      haserr = true;
    }
  }
  if (!haserr)
  {
    tyrindex->numofmers = (numofbytes - rest) / tyrindex->merbytes;
    assert(tyrindex->mertable != NULL);
    if (tyrindex->numofmers == 0)
    {
      tyrindex->lastmer = tyrindex->mertable - 1;
    } else
    {
      tyrindex->lastmer
        = tyrindex->mertable + numofbytes - rest - tyrindex->merbytes;
    }
  }
  if (haserr)
  {
    if (tyrindex->mappedfileptr != NULL)
    {
      gt_fa_xmunmap(tyrindex->mappedfileptr);
      tyrindex->mappedfileptr = NULL;
    }
    FREESPACE(tyrindex);
  }
  return haserr ? NULL : tyrindex;
}

void tyrindex_show(const Tyrindex *tyrindex)
{
  printf("# indexfilename = %s\n",gt_str_get(tyrindex->indexfilename));
  printf("# alphasize = %u\n",tyrindex->alphasize);
  printf("# mersize = %lu\n",tyrindex->mersize);
  printf("# numofmers = %lu\n",(unsigned long) tyrindex->numofmers);
  printf("# merbytes = %lu\n",tyrindex->merbytes);
}

void tyrindex_delete(Tyrindex **tyrindexptr)
{
  Tyrindex *tyrindex = *tyrindexptr;

  gt_fa_xmunmap(tyrindex->mappedfileptr);
  tyrindex->mappedfileptr = NULL;
  FREESPACE(tyrindex);
  *tyrindexptr = NULL;
}

Tyrcountinfo *tyrcountinfo_new(size_t numofmers,
                               const GtStr *tyrindexname,
                               GtError *err)
{
  size_t numofbytes;
  void *tmp;
  bool haserr = false;
  Tyrcountinfo *tyrcountinfo;

  ALLOCASSIGNSPACE(tyrcountinfo,NULL,Tyrcountinfo,1);
  tyrcountinfo->indexfilename = tyrindexname;
  tyrcountinfo->mappedmctfileptr
    = genericmaponlytable(tyrindexname,COUNTSSUFFIX,&numofbytes,err);
  if (tyrcountinfo->mappedmctfileptr == NULL)
  {
    tyrcountinfo->smallcounts = NULL;
    haserr = true;
  } else
  {
    tyrcountinfo->smallcounts
      = (Uchar *) tyrcountinfo->mappedmctfileptr;
    tmp = &tyrcountinfo->smallcounts[numofmers];
    tyrcountinfo->largecounts = (Largecount *) tmp;
    if (numofbytes < numofmers)
    {
      gt_error_set(err,"size of file \"%s.%s\" is smaller than minimum size "
                       "%lu",gt_str_get(tyrindexname),COUNTSSUFFIX,
                       (unsigned long) numofmers);
      haserr = true;
    }
  }
  if (!haserr && (numofbytes - numofmers) % sizeof (Largecount) != 0)
  {
    gt_error_set(err,"(numofbytes - numofmers) = %lu must be a multiple of %lu",
           (unsigned long) (numofbytes - numofmers),
           (unsigned long) sizeof (Largecount));
    haserr = true;
  }
  if (!haserr)
  {
    tyrcountinfo->numoflargecounts
      = (unsigned long) (numofbytes - numofmers)/
        (unsigned long) sizeof (Largecount);
  }
  if (haserr)
  {
    if (tyrcountinfo->mappedmctfileptr != NULL)
    {
      gt_fa_xmunmap(tyrcountinfo->mappedmctfileptr);
      tyrcountinfo->mappedmctfileptr = NULL;
    }
    FREESPACE(tyrcountinfo);
  }
  return haserr ? NULL : tyrcountinfo;
}

void tyrcountinfo_delete(Tyrcountinfo **tyrcountinfoptr)
{
  Tyrcountinfo *tyrcountinfo = *tyrcountinfoptr;

  gt_fa_xmunmap(tyrcountinfo->mappedmctfileptr);
  tyrcountinfo->mappedmctfileptr = NULL;
  FREESPACE(tyrcountinfo);
  *tyrcountinfoptr = NULL;
}

static int mymemcmp(unsigned long *offset,const Uchar *s1,const Uchar *s2,
                    unsigned long len)
{
  unsigned long idx;

  for (idx=*offset; idx<len; idx++)
  {
    if (s1[idx] < s2[idx])
    {
      *offset = idx;
      return -1;
    }
    if (s1[idx] > s2[idx])
    {
      *offset = idx;
      return 1;
    }
  }
  return 0;
}

static /*@null@*/ const Uchar *binmersearch(unsigned long merbytes,
                                            unsigned long offset,
                                            const Uchar *key,
                                            const Uchar *left,
                                            const Uchar *right)
{
  const Uchar *leftptr = left,
              *midptr,
              *rightptr = right;
  int cmpval;
  unsigned long leftlength = offset, rightlength = offset, len;

  while (leftptr <= rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/MULT2(merbytes);
    midptr = leftptr + merbytes * len;
    cmpval = mymemcmp(&offset,midptr,key,merbytes);
    if (cmpval < 0)
    {
      leftptr = midptr + merbytes;
      leftlength = offset;
      if (offset > rightlength)
      {
        offset = rightlength;
      }
    } else
    {
      if (cmpval > 0)
      {
        rightptr = midptr - merbytes;
        rightlength = offset;
        if (offset > leftlength)
        {
          offset = leftlength;
        }
      } else
      {
        return midptr;
      }
    }
  }
  return NULL;
}

static void checktyrindex(const Tyrindex *tyrindex)
{
  Uchar *mercodeptr;
  const Uchar *result;
  unsigned long position, previousposition = 0;

  for (mercodeptr = tyrindex->mertable;
       mercodeptr <= tyrindex->lastmer;
       mercodeptr += tyrindex->merbytes)
  {
    result = binmersearch(tyrindex->merbytes,
                          0,
                          mercodeptr,
                          tyrindex->mertable,
                          tyrindex->lastmer);
    assert(result != NULL);
    if ((result - tyrindex->mertable) % tyrindex->merbytes != 0)
    {
      fprintf(stderr,"result is not multiple of %lu\n",tyrindex->merbytes);
      exit(EXIT_FAILURE);
    }
    position = (unsigned long) (result-tyrindex->mertable)/
                               tyrindex->merbytes;
    if (position > 0)
    {
      if (previousposition + 1 != position)
      {
        fprintf(stderr,"position %lu is not increasing\n",position);
        exit(EXIT_FAILURE);
      }
    }
    previousposition = position;
  }
}

static unsigned long containsspecialbytestring(const Uchar *seq,
                                               unsigned long len)
{
  const Uchar *sptr;

  for (sptr=seq; sptr < seq + len; sptr++)
  {
    if (ISSPECIAL(*sptr))
    {
      return (unsigned long) (sptr - seq);
    }
  }
  return len;
}

typedef struct
{
  Uchar *bytecode,  /* buffer for encoded word to be searched */
        *rcbuf;
  unsigned int showmode,
               searchstrand;
  Alphabet *dnaalpha;
} Tyrsearchinfo;

static void tyrsearchinfo_init(Tyrsearchinfo *tyrsearchinfo,
                               const Tyrindex *tyrindex,
                               unsigned int showmode,
                               unsigned int searchstrand)
{
  tyrsearchinfo->showmode = showmode;
  tyrsearchinfo->searchstrand = searchstrand;
  tyrsearchinfo->dnaalpha = assigninputalphabet(true,false,NULL,NULL,NULL);
  ALLOCASSIGNSPACE(tyrsearchinfo->bytecode,NULL,Uchar,
                   tyrindex->merbytes);
  ALLOCASSIGNSPACE(tyrsearchinfo->rcbuf,NULL,Uchar,tyrindex->mersize);
}

void tyrsearchinfo_delete(Tyrsearchinfo *tyrsearchinfo)
{
  freeAlphabet(&tyrsearchinfo->dnaalpha);
  FREESPACE(tyrsearchinfo->bytecode);
  FREESPACE(tyrsearchinfo->rcbuf);
}

/*@null@*/ const Uchar *searchsinglemer(const Uchar *qptr,
                                        const Tyrindex *tyrindex,
                                        const Tyrsearchinfo *tyrsearchinfo)
{
  const Uchar *result;

  plainseq2bytecode(tyrsearchinfo->bytecode,qptr,tyrindex->mersize);
  result = binmersearch(tyrindex->merbytes,
                        0,
                        tyrsearchinfo->bytecode,
                        tyrindex->mertable,
                        tyrindex->lastmer);
  return result;
}

#define ADDTABULATOR\
        if (firstitem)\
        {\
          firstitem = false;\
        } else\
        {\
          (void) putchar('\t');\
        }

static void mermatchoutput(const Tyrindex *tyrindex,
                           const Tyrsearchinfo *tyrsearchinfo,
                           GT_UNUSED const Uchar *result,
                           const Uchar *query,
                           const Uchar *qptr,
                           uint64_t unitnum,
                           bool forward)
{
  bool firstitem = true;
  unsigned long queryposition;

  queryposition = (unsigned long) (qptr-query);
  if (tyrsearchinfo->showmode & SHOWQSEQNUM)
  {
    printf(Formatuint64_t,PRINTuint64_tcast(unitnum));
    firstitem = false;
  }
  if (tyrsearchinfo->showmode & SHOWQPOS)
  {
    ADDTABULATOR;
    printf("%c%lu",forward ? '+' : '-',queryposition);
  }
  if (tyrsearchinfo->showmode & SHOWCOUNTS)
  {
    /*
    unsigned long mernumber
      = (unsigned long) (result - tyrindex->mertable)/
                        tyrindex->merbytes;
    ADDTABULATOR;
    if (showmercounts(&tyrsearchinfo->mctinfo,mernumber) != 0)
    {
      return (Sint) -1;
    }
    */
  }
  if (tyrsearchinfo->showmode & SHOWSEQUENCE)
  {
    ADDTABULATOR;
    fprintfsymbolstring(stdout,tyrsearchinfo->dnaalpha,qptr,
                        tyrindex->mersize);
  }
  if (tyrsearchinfo->showmode & (SHOWSEQUENCE | SHOWQPOS | SHOWCOUNTS))
  {
    (void) putchar('\n');
  }
}

static void singleseqtyrsearch(const Tyrindex *tyrindex,
                               const Tyrsearchinfo *tyrsearchinfo,
                               uint64_t unitnum,
                               const Uchar *query,
                               unsigned long querylen,
                               GT_UNUSED const char *desc)
{
  const Uchar *qptr, *result;
  unsigned long skipvalue;

  if (tyrindex->mersize > querylen)
  {
    return;
  }
  qptr = query;
  while (qptr <= query + querylen - tyrindex->mersize)
  {
    skipvalue = containsspecialbytestring(qptr,tyrindex->mersize);
    if (skipvalue == tyrindex->mersize)
    {
      if (tyrsearchinfo->searchstrand & STRAND_FORWARD)
      {
        result = searchsinglemer(qptr,tyrindex,tyrsearchinfo);
        if (result != NULL)
        {
          mermatchoutput(tyrindex,
                         tyrsearchinfo,
                         result,
                         query,
                         qptr,
                         unitnum,
                         true);
        }
      }
      if (tyrsearchinfo->searchstrand & STRAND_REVERSE)
      {
        assert(tyrsearchinfo->rcbuf != NULL);
        copy_reversecomplement(tyrsearchinfo->rcbuf,qptr,
                               tyrindex->mersize);
        result = searchsinglemer(tyrsearchinfo->rcbuf,tyrindex,
                                 tyrsearchinfo);
        if (result != NULL)
        {
          mermatchoutput(tyrindex,
                         tyrsearchinfo,
                         result,
                         query,
                         qptr,
                         unitnum,
                         false);
        }
      }
      qptr++;
    } else
    {
      qptr += (skipvalue+1);
    }
  }
}

int tyrsearch(const GtStr *tyrindexname,
              const GtStrArray *queryfilenames,
              unsigned int showmode,
              unsigned int searchstrand,
              bool verbose,
              bool performtest,
              GtError *err)
{
  Tyrindex *tyrindex;
  Tyrcountinfo *tyrcountinfo = NULL;
  bool haserr = false;
  GtSeqIterator *seqit;

  tyrindex = tyrindex_new(tyrindexname,err);
  if (tyrindex == NULL)
  {
    haserr = true;
  } else
  {
    if (verbose)
    {
      tyrindex_show(tyrindex);
    }
    if (performtest)
    {
      checktyrindex(tyrindex);
    }
  }
  if (!haserr)
  {
    assert(tyrindex != NULL);
    if ((showmode & SHOWCOUNTS) && tyrindex->numofmers > 0)
    {
      tyrcountinfo = tyrcountinfo_new(tyrindex->numofmers,tyrindexname,err);
      if (tyrcountinfo == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Tyrsearchinfo tyrsearchinfo;

    assert(tyrindex != NULL);
    tyrsearchinfo_init(&tyrsearchinfo,tyrindex,showmode,searchstrand);
    seqit = gt_seqiterator_new(queryfilenames,
                               getsymbolmapAlphabet(tyrsearchinfo.dnaalpha),
                               true);
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = gt_seqiterator_next(seqit,
                                   &query,
                                   &querylen,
                                   &desc,
                                   err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      singleseqtyrsearch(tyrindex,
                              &tyrsearchinfo,
                              unitnum,
                              query,
                              querylen,
                              desc);
      gt_free(desc);
    }
    gt_seqiterator_delete(seqit);
    tyrsearchinfo_delete(&tyrsearchinfo);
  }
  if (tyrcountinfo != NULL)
  {
    tyrcountinfo_delete(&tyrcountinfo);
  }
  if (tyrindex != NULL)
  {
    tyrindex_delete(&tyrindex);
  }
  return haserr ? -1 : 0;
}
