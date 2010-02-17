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

#include "core/alphabet.h"
#include "core/fa.h"
#include "core/unused_api.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/chardef.h"
#include "revcompl.h"
#include "format64.h"
#include "encseq-def.h"
#include "tyr-map.h"
#include "tyr-search.h"
#include "tyr-show.h"
#include "tyr-mersplit.h"
#include "spacedef.h"

typedef struct
{
  GtUchar *bytecode,  /* buffer for encoded word to be searched */
        *rcbuf;
  const GtUchar *mertable, *lastmer;
  unsigned long mersize;
  unsigned int showmode,
               searchstrand;
  GtAlphabet *dnaalpha;
} Tyrsearchinfo;

static void tyrsearchinfo_init(Tyrsearchinfo *tyrsearchinfo,
                               const Tyrindex *tyrindex,
                               unsigned int showmode,
                               unsigned int searchstrand)
{
  unsigned long merbytes;

  merbytes = tyrindex_merbytes(tyrindex);
  tyrsearchinfo->mersize = tyrindex_mersize(tyrindex);
  tyrsearchinfo->mertable = tyrindex_mertable(tyrindex);
  tyrsearchinfo->lastmer = tyrindex_lastmer(tyrindex);
  tyrsearchinfo->showmode = showmode;
  tyrsearchinfo->searchstrand = searchstrand;
  tyrsearchinfo->dnaalpha = gt_alphabet_new(true,false,NULL,NULL,NULL);
  ALLOCASSIGNSPACE(tyrsearchinfo->bytecode,NULL,GtUchar,merbytes);
  ALLOCASSIGNSPACE(tyrsearchinfo->rcbuf,NULL,GtUchar,tyrsearchinfo->mersize);
}

void tyrsearchinfo_delete(Tyrsearchinfo *tyrsearchinfo)
{
  gt_alphabet_delete(tyrsearchinfo->dnaalpha);
  FREESPACE(tyrsearchinfo->bytecode);
  FREESPACE(tyrsearchinfo->rcbuf);
}

/*@null@*/ const GtUchar *searchsinglemer(const GtUchar *qptr,
                                        const Tyrindex *tyrindex,
                                        const Tyrsearchinfo *tyrsearchinfo,
                                        const Tyrbckinfo *tyrbckinfo)
{
  const GtUchar *result;

  plainseq2bytecode(tyrsearchinfo->bytecode,qptr,tyrsearchinfo->mersize);
  if (tyrbckinfo == NULL)
  {
    result = tyrindex_binmersearch(tyrindex,0,tyrsearchinfo->bytecode,
                                   tyrsearchinfo->mertable,
                                   tyrsearchinfo->lastmer);
  } else
  {
    result = searchinbuckets(tyrindex,tyrbckinfo,tyrsearchinfo->bytecode);
  }
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
                           const Tyrcountinfo *tyrcountinfo,
                           const Tyrsearchinfo *tyrsearchinfo,
                           const GtUchar *result,
                           const GtUchar *query,
                           const GtUchar *qptr,
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
    unsigned long mernumber = tyrindex_ptr2number(tyrindex,result);
    ADDTABULATOR;
    printf("%lu",tyrcountinfo_get(tyrcountinfo,mernumber));
  }
  if (tyrsearchinfo->showmode & SHOWSEQUENCE)
  {
    ADDTABULATOR;
    gt_alphabet_printf_symbolstring(tyrsearchinfo->dnaalpha,qptr,
                                    tyrsearchinfo->mersize);
  }
  if (tyrsearchinfo->showmode & (SHOWSEQUENCE | SHOWQPOS | SHOWCOUNTS))
  {
    (void) putchar('\n');
  }
}

static void singleseqtyrsearch(const Tyrindex *tyrindex,
                               const Tyrcountinfo *tyrcountinfo,
                               const Tyrsearchinfo *tyrsearchinfo,
                               const Tyrbckinfo *tyrbckinfo,
                               uint64_t unitnum,
                               const GtUchar *query,
                               unsigned long querylen,
                               GT_UNUSED const char *desc)
{
  const GtUchar *qptr, *result;
  unsigned long offset, skipvalue;

  if (tyrsearchinfo->mersize > querylen)
  {
    return;
  }
  qptr = query;
  offset = 0;
  while (qptr <= query + querylen - tyrsearchinfo->mersize)
  {
    skipvalue = containsspecialbytestring(qptr,offset,tyrsearchinfo->mersize);
    if (skipvalue == tyrsearchinfo->mersize)
    {
      offset = tyrsearchinfo->mersize-1;
      if (tyrsearchinfo->searchstrand & STRAND_FORWARD)
      {
        result = searchsinglemer(qptr,tyrindex,tyrsearchinfo,tyrbckinfo);
        if (result != NULL)
        {
          mermatchoutput(tyrindex,
                         tyrcountinfo,
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
        gt_assert(tyrsearchinfo->rcbuf != NULL);
        copy_reversecomplement(tyrsearchinfo->rcbuf,qptr,
                               tyrsearchinfo->mersize);
        result = searchsinglemer(tyrsearchinfo->rcbuf,tyrindex,
                                 tyrsearchinfo,tyrbckinfo);
        if (result != NULL)
        {
          mermatchoutput(tyrindex,
                         tyrcountinfo,
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
      offset = 0;
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
  Tyrbckinfo *tyrbckinfo = NULL;
  bool haserr = false;

  gt_error_check(err);
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
      tyrindex_check(tyrindex);
    }
  }
  if (!haserr)
  {
    gt_assert(tyrindex != NULL);
    if ((showmode & SHOWCOUNTS) && !tyrindex_isempty(tyrindex))
    {
      tyrcountinfo = tyrcountinfo_new(tyrindex,tyrindexname,err);
      if (tyrcountinfo == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    gt_assert(tyrindex != NULL);
    if (!tyrindex_isempty(tyrindex))
    {
      tyrbckinfo = tyrbckinfo_new(tyrindexname,tyrindex_alphasize(tyrindex),
                                  err);
      if (tyrbckinfo == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    const GtUchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Tyrsearchinfo tyrsearchinfo;
    GtSeqIterator *seqit;

    gt_assert(tyrindex != NULL);
    tyrsearchinfo_init(&tyrsearchinfo,tyrindex,showmode,searchstrand);
    seqit = gt_seqiterator_sequence_buffer_new(queryfilenames, err);
    if (!seqit)
      haserr = true;
    if (!haserr)
    {
      gt_seqiterator_set_symbolmap(seqit,
                                 gt_alphabet_symbolmap(tyrsearchinfo.dnaalpha));
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
                           tyrcountinfo,
                           &tyrsearchinfo,
                           tyrbckinfo,
                           unitnum,
                           query,
                           querylen,
                           desc);
        gt_free(desc);
      }
      gt_seqiterator_delete(seqit);
    }
    tyrsearchinfo_delete(&tyrsearchinfo);
  }
  if (tyrbckinfo != NULL)
  {
    tyrbckinfo_delete(&tyrbckinfo);
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
