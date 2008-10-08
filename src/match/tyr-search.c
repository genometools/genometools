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
#include "match/alphadef.h"
#include "divmodmul.h"
#include "opensfxfile.h"
#include "tyr-search.h"
#include "spacedef.h"

struct Tallymercountinfo
{
  void *mappedmctfileptr;
  const GtStr *indexfilename;
  Uchar *smallcounts;
  Largecount *largecounts;
  unsigned long numoflargecounts;
};

struct Tallymerindex
{
  void *mappedfileptr;
  const GtStr *indexfilename;
  unsigned int alphasize;
  unsigned long numofmers,
                mersize,
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

Tallymerindex *tallymerindex_new(const GtStr *tallymerindexname,GtError *err)
{
  bool haserr = false;
  size_t numofbytes, rest;
  Tallymerindex *tallymerindex;

  ALLOCASSIGNSPACE(tallymerindex,NULL,Tallymerindex,1);
  tallymerindex->indexfilename = tallymerindexname;
  tallymerindex->mappedfileptr = genericmaponlytable(tallymerindexname,
                                                 MERSUFFIX,&numofbytes,err);
  if (tallymerindex->mappedfileptr == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    tallymerindex->mertable = (Uchar *) tallymerindex->mappedfileptr;
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
    assert(tallymerindex->mertable != NULL);
    tallymerindex->mersize
      = decodesingleinteger(tallymerindex->mertable + numofbytes - rest);
    tallymerindex->alphasize
      = (unsigned int) decodesingleinteger(tallymerindex->mertable +
                                           numofbytes - rest +
                                           sizeof (unsigned long));
    tallymerindex->merbytes = MERBYTES(tallymerindex->mersize);
    if ((numofbytes - rest) % tallymerindex->merbytes != 0)
    {
      gt_error_set(err,"size of index is %lu which is not a multiple of %lu",
                   (unsigned long) (numofbytes - rest),
                   tallymerindex->merbytes);
      haserr = true;
    }
  }
  if (!haserr)
  {
    tallymerindex->numofmers
      = (unsigned long) (numofbytes - rest) / tallymerindex->merbytes;
    assert(tallymerindex->mertable != NULL);
    if (tallymerindex->numofmers == 0)
    {
      tallymerindex->lastmer = tallymerindex->mertable - 1;
    } else
    {
      tallymerindex->lastmer
        = tallymerindex->mertable + numofbytes - rest - tallymerindex->merbytes;
    }
  }
  if (haserr)
  {
    if (tallymerindex->mappedfileptr != NULL)
    {
      gt_fa_xmunmap(tallymerindex->mappedfileptr);
      tallymerindex->mappedfileptr = NULL;
    }
    FREESPACE(tallymerindex);
  }
  return haserr ? NULL : tallymerindex;
}

void tallymerindex_show(const Tallymerindex *tallymerindex)
{
  printf("# indexfilename = %s\n",gt_str_get(tallymerindex->indexfilename));
  printf("# alphasize = %u\n",tallymerindex->alphasize);
  printf("# mersize = %lu\n",tallymerindex->mersize);
  printf("# numofmers = %lu\n",tallymerindex->numofmers);
  printf("# merbytes = %lu\n",tallymerindex->merbytes);
}

void tallymerindex_delete(Tallymerindex **tallymerindexptr)
{
  Tallymerindex *tallymerindex = *tallymerindexptr;

  gt_fa_xmunmap(tallymerindex->mappedfileptr);
  tallymerindex->mappedfileptr = NULL;
  FREESPACE(tallymerindex);
  *tallymerindexptr = NULL;
}

Tallymercountinfo *tallymercountinfo_new(size_t numofmers,
                                         const GtStr *tallymerindexname,
                                         GtError *err)
{
  size_t numofbytes;
  void *tmp;
  bool haserr = false;
  Tallymercountinfo *tallymercountinfo;

  ALLOCASSIGNSPACE(tallymercountinfo,NULL,Tallymercountinfo,1);
  tallymercountinfo->indexfilename = tallymerindexname;
  tallymercountinfo->mappedmctfileptr
    = genericmaponlytable(tallymerindexname,COUNTSSUFFIX,&numofbytes,err);
  if (tallymercountinfo->mappedmctfileptr == NULL)
  {
    tallymercountinfo->smallcounts = NULL;
    haserr = true;
  } else
  {
    tallymercountinfo->smallcounts
      = (Uchar *) tallymercountinfo->mappedmctfileptr;
    tmp = &tallymercountinfo->smallcounts[numofmers];
    tallymercountinfo->largecounts = (Largecount *) tmp;
    if (numofbytes < numofmers)
    {
      gt_error_set(err,"size of file \"%s.%s\" is smaller than minimum size "
                       "%lu",gt_str_get(tallymerindexname),COUNTSSUFFIX,
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
    tallymercountinfo->numoflargecounts
      = (unsigned long) (numofbytes - numofmers)/
        (unsigned long) sizeof (Largecount);
  }
  if (haserr)
  {
    if (tallymercountinfo->mappedmctfileptr != NULL)
    {
      gt_fa_xmunmap(tallymercountinfo->mappedmctfileptr);
      tallymercountinfo->mappedmctfileptr = NULL;
    }
    FREESPACE(tallymercountinfo);
  }
  return haserr ? NULL : tallymercountinfo;
}

void tallymercountinfo_delete(Tallymercountinfo **tallymercountinfoptr)
{
  Tallymercountinfo *tallymercountinfo = *tallymercountinfoptr;

  gt_fa_xmunmap(tallymercountinfo->mappedmctfileptr);
  tallymercountinfo->mappedmctfileptr = NULL;
  FREESPACE(tallymercountinfo);
  *tallymercountinfoptr = NULL;
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

static void checktallymerindex(const Tallymerindex *tallymerindex)
{
  Uchar *mercodeptr;
  const Uchar *result;
  unsigned long position, previousposition = 0;

  for (mercodeptr = tallymerindex->mertable;
       mercodeptr <= tallymerindex->lastmer;
       mercodeptr += tallymerindex->merbytes)
  {
    result = binmersearch(tallymerindex->merbytes,
                          0,
                          mercodeptr,
                          tallymerindex->mertable,
                          tallymerindex->lastmer);
    assert(result != NULL);
    if ((result - tallymerindex->mertable) % tallymerindex->merbytes != 0)
    {
      fprintf(stderr,"result is not multiple of %lu\n",tallymerindex->merbytes);
      exit(EXIT_FAILURE);
    }
    position = (unsigned long) (result-tallymerindex->mertable)/
                               tallymerindex->merbytes;
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

static void singleseqtallymersearch(unsigned long mersize,
                                    GT_UNUSED uint64_t unitnum,
                                    const Uchar *query,
                                    unsigned long querylen,
                                    GT_UNUSED const char *desc)
{
  const Uchar *qptr;

  for (qptr = query; qptr <= query + querylen - mersize; qptr++)
  {
  }
}

int tallymersearch(const GtStr *tallymerindexname,
                   const GtStrArray *queryfilenames,
                   GT_UNUSED unsigned int showmode,
                   GT_UNUSED unsigned int strand,
                   bool verbose,
                   bool performtest,
                   GtError *err)
{
  Tallymerindex *tallymerindex;
  bool haserr = false;
  GtSeqIterator *seqit;
  Alphabet *dnaalpha = NULL;

  tallymerindex = tallymerindex_new(tallymerindexname,err);
  if (tallymerindex == NULL)
  {
    haserr = true;
  } else
  {
    if (verbose)
    {
      tallymerindex_show(tallymerindex);
    }
    if (performtest)
    {
      checktallymerindex(tallymerindex);
    }
    dnaalpha = assigninputalphabet(true,false,NULL,NULL,err);
    if (dnaalpha == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;

    assert(tallymerindex != NULL);
    seqit = gt_seqiterator_new(queryfilenames,
                               getsymbolmapAlphabet(dnaalpha),
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
      singleseqtallymersearch(tallymerindex->mersize,
                              unitnum,
                              query,
                              querylen,
                              desc);
      FREESPACE(desc);
    }
    gt_seqiterator_delete(seqit);
  }
  if (tallymerindex != NULL)
  {
    tallymerindex_delete(&tallymerindex);
  }
  if (dnaalpha != NULL)
  {
    freeAlphabet(&dnaalpha);
  }
  return haserr ? -1 : 0;
}
