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
} Vmersearchinfo;

static void vmersearchinfo_init(Vmersearchinfo *vmersearchinfo,
                                const Tallymerindex *tallymerindex,
                                unsigned int showmode,
                                unsigned int searchstrand)
{
  vmersearchinfo->showmode = showmode;
  vmersearchinfo->searchstrand = searchstrand;
  vmersearchinfo->dnaalpha = assigninputalphabet(true,false,NULL,NULL,NULL);
  ALLOCASSIGNSPACE(vmersearchinfo->bytecode,NULL,Uchar,tallymerindex->merbytes);
  ALLOCASSIGNSPACE(vmersearchinfo->rcbuf,NULL,Uchar,tallymerindex->mersize);
}

void vmersearchinfo_delete(Vmersearchinfo *vmersearchinfo)
{
  freeAlphabet(&vmersearchinfo->dnaalpha);
  FREESPACE(vmersearchinfo->bytecode);
  FREESPACE(vmersearchinfo->rcbuf);
}

/*@null@*/ const Uchar *searchsinglemer(const Uchar *qptr,
                                        const Tallymerindex *tallymerindex,
                                        const Vmersearchinfo *vmersearchinfo)
{
  const Uchar *result;

  plainseq2bytecode(vmersearchinfo->bytecode,qptr,tallymerindex->mersize);
  result = binmersearch(tallymerindex->merbytes,
                        0,
                        vmersearchinfo->bytecode,
                        tallymerindex->mertable,
                        tallymerindex->lastmer);
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

static void mermatchoutput(const Tallymerindex *tallymerindex,
                           const Vmersearchinfo *vmersearchinfo,
                           GT_UNUSED const Uchar *result,
                           const Uchar *query,
                           const Uchar *qptr,
                           uint64_t unitnum,
                           bool forward)
{
  bool firstitem = true;
  unsigned long queryposition;

  queryposition = (unsigned long) (qptr-query);
  if (vmersearchinfo->showmode & SHOWQSEQNUM)
  {
    printf(Formatuint64_t,PRINTuint64_tcast(unitnum));
    firstitem = false;
  }
  if (vmersearchinfo->showmode & SHOWQPOS)
  {
    ADDTABULATOR;
    printf("%c%lu",forward ? '+' : '-',queryposition);
  }
  if (vmersearchinfo->showmode & SHOWCOUNTS)
  {
    /*
    unsigned long mernumber
      = (unsigned long) (result - tallymerindex->mertable)/
                        tallymerindex->merbytes;
    ADDTABULATOR;
    if (showmercounts(&vmersearchinfo->mctinfo,mernumber) != 0)
    {
      return (Sint) -1;
    }
    */
  }
  if (vmersearchinfo->showmode & SHOWSEQUENCE)
  {
    ADDTABULATOR;
    fprintfsymbolstring(stdout,vmersearchinfo->dnaalpha,qptr,
                        tallymerindex->mersize);
  }
  if (vmersearchinfo->showmode & (SHOWSEQUENCE | SHOWQPOS | SHOWCOUNTS))
  {
    printf("\n");
  }
}

static void singleseqtallymersearch(const Tallymerindex *tallymerindex,
                                    const Vmersearchinfo *vmersearchinfo,
                                    uint64_t unitnum,
                                    const Uchar *query,
                                    unsigned long querylen,
                                    GT_UNUSED const char *desc)
{
  const Uchar *qptr, *result;
  unsigned long skipvalue;

  if (tallymerindex->mersize > querylen)
  {
    return;
  }
  qptr = query;
  while (qptr <= query + querylen - tallymerindex->mersize)
  {
    skipvalue = containsspecialbytestring(qptr,tallymerindex->mersize);
    if (skipvalue == tallymerindex->mersize)
    {
      if (vmersearchinfo->searchstrand & STRAND_FORWARD)
      {
        result = searchsinglemer(qptr,tallymerindex,vmersearchinfo);
        if (result != NULL)
        {
          mermatchoutput(tallymerindex,
                         vmersearchinfo,
                         result,
                         query,
                         qptr,
                         unitnum,
                         true);
        }
      }
      if (vmersearchinfo->searchstrand & STRAND_REVERSE)
      {
        assert(vmersearchinfo->rcbuf != NULL);
        copy_reversecomplement(vmersearchinfo->rcbuf,qptr,
                               tallymerindex->mersize);
        result = searchsinglemer(vmersearchinfo->rcbuf,tallymerindex,
                                 vmersearchinfo);
        if (result != NULL)
        {
          mermatchoutput(tallymerindex,
                         vmersearchinfo,
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

int tallymersearch(const GtStr *tallymerindexname,
                   const GtStrArray *queryfilenames,
                   unsigned int showmode,
                   unsigned int searchstrand,
                   bool verbose,
                   bool performtest,
                   GtError *err)
{
  Tallymerindex *tallymerindex;
  bool haserr = false;
  GtSeqIterator *seqit;

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
  }
  if (!haserr)
  {
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Vmersearchinfo vmersearchinfo;

    assert(tallymerindex != NULL);
    vmersearchinfo_init(&vmersearchinfo,tallymerindex,showmode,searchstrand);
    seqit = gt_seqiterator_new(queryfilenames,
                               getsymbolmapAlphabet(vmersearchinfo.dnaalpha),
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
      singleseqtallymersearch(tallymerindex,
                              &vmersearchinfo,
                              unitnum,
                              query,
                              querylen,
                              desc);
      gt_free(desc);
    }
    gt_seqiterator_delete(seqit);
    vmersearchinfo_delete(&vmersearchinfo);
  }
  if (tallymerindex != NULL)
  {
   tallymerindex_delete(&tallymerindex);
  }
  return haserr ? -1 : 0;
}
