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
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/chardef.h"
#include "core/format64.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "revcompl.h"
#include "tyr-map.h"
#include "tyr-search.h"
#include "tyr-show.h"
#include "tyr-mersplit.h"

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

static void gt_tyrsearchinfo_init(Tyrsearchinfo *tyrsearchinfo,
                               const Tyrindex *tyrindex,
                               unsigned int showmode,
                               unsigned int searchstrand)
{
  unsigned long merbytes;

  merbytes = gt_tyrindex_merbytes(tyrindex);
  tyrsearchinfo->mersize = gt_tyrindex_mersize(tyrindex);
  tyrsearchinfo->mertable = gt_tyrindex_mertable(tyrindex);
  tyrsearchinfo->lastmer = gt_tyrindex_lastmer(tyrindex);
  tyrsearchinfo->showmode = showmode;
  tyrsearchinfo->searchstrand = searchstrand;
  tyrsearchinfo->dnaalpha = gt_alphabet_new_dna();
  tyrsearchinfo->bytecode = gt_malloc(sizeof *tyrsearchinfo->bytecode
                                      * merbytes);
  tyrsearchinfo->rcbuf = gt_malloc(sizeof *tyrsearchinfo->rcbuf
                                   * tyrsearchinfo->mersize);
}

static void gt_tyrsearchinfo_delete(Tyrsearchinfo *tyrsearchinfo)
{
  if (tyrsearchinfo != NULL)
  {
    gt_alphabet_delete(tyrsearchinfo->dnaalpha);
    gt_free(tyrsearchinfo->bytecode);
    gt_free(tyrsearchinfo->rcbuf);
  }
}

/*@null@*/ const GtUchar *gt_searchsinglemer(const GtUchar *qptr,
                                        const Tyrindex *tyrindex,
                                        const Tyrsearchinfo *tyrsearchinfo,
                                        const Tyrbckinfo *tyrbckinfo)
{
  const GtUchar *result;

  gt_encseq_plainseq2bytecode(tyrsearchinfo->bytecode,qptr,
                                       tyrsearchinfo->mersize);
  if (tyrbckinfo == NULL)
  {
    result = gt_tyrindex_binmersearch(tyrindex,0,tyrsearchinfo->bytecode,
                                   tyrsearchinfo->mertable,
                                   tyrsearchinfo->lastmer);
  } else
  {
    result = gt_searchinbuckets(tyrindex,tyrbckinfo,tyrsearchinfo->bytecode);
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
    unsigned long mernumber = gt_tyrindex_ptr2number(tyrindex,result);
    ADDTABULATOR;
    printf("%lu",gt_tyrcountinfo_get(tyrcountinfo,mernumber));
  }
  if (tyrsearchinfo->showmode & SHOWSEQUENCE)
  {
    ADDTABULATOR;
    gt_alphabet_decode_seq_to_fp(tyrsearchinfo->dnaalpha,
                                 stdout,
                                 qptr,
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
        result = gt_searchsinglemer(qptr,tyrindex,tyrsearchinfo,tyrbckinfo);
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
        gt_copy_reversecomplement(tyrsearchinfo->rcbuf,qptr,
                               tyrsearchinfo->mersize);
        result = gt_searchsinglemer(tyrsearchinfo->rcbuf,tyrindex,
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

int gt_tyrsearch(const char *tyrindexname,
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
  tyrindex = gt_tyrindex_new(tyrindexname,err);
  if (tyrindex == NULL)
  {
    haserr = true;
  } else
  {
    if (verbose)
    {
      gt_tyrindex_show(tyrindex);
    }
    if (performtest)
    {
      gt_tyrindex_check(tyrindex);
    }
  }
  if (!haserr)
  {
    gt_assert(tyrindex != NULL);
    if ((showmode & SHOWCOUNTS) && !gt_tyrindex_isempty(tyrindex))
    {
      tyrcountinfo = gt_tyrcountinfo_new(tyrindex,tyrindexname,err);
      if (tyrcountinfo == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    gt_assert(tyrindex != NULL);
    if (!gt_tyrindex_isempty(tyrindex))
    {
      tyrbckinfo = gt_tyrbckinfo_new(tyrindexname,
                                     gt_tyrindex_alphasize(tyrindex),
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
    gt_tyrsearchinfo_init(&tyrsearchinfo,tyrindex,showmode,searchstrand);
    seqit = gt_seq_iterator_sequence_buffer_new(queryfilenames, err);
    if (!seqit)
      haserr = true;
    if (!haserr)
    {
      gt_seq_iterator_set_symbolmap(seqit,
                                 gt_alphabet_symbolmap(tyrsearchinfo.dnaalpha));
      for (unitnum = 0; /* Nothing */; unitnum++)
      {
        retval = gt_seq_iterator_next(seqit,
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
      }
      gt_seq_iterator_delete(seqit);
    }
    gt_tyrsearchinfo_delete(&tyrsearchinfo);
  }
  if (tyrbckinfo != NULL)
  {
    gt_tyrbckinfo_delete(&tyrbckinfo);
  }
  if (tyrcountinfo != NULL)
  {
    gt_tyrcountinfo_delete(&tyrcountinfo);
  }
  if (tyrindex != NULL)
  {
    gt_tyrindex_delete(&tyrindex);
  }
  return haserr ? -1 : 0;
}
