/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include "core/assert_api.h"
#include "core/defined-types.h"
#include "core/logger_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/readmode.h"
#include "core/encseq.h"
#include "sfx-suffixgetset.h"

struct GtSuffixsortspace
{
  bool currentexport;
  unsigned int clonenumber; /* 0 for no clone, > 0 otherwise */
  Definedunsignedlong longestidx;
  GtSuffixsortspace_exportptr exportptr;
  GtUword maxindex,
          maxvalue,
          partoffset,
          bucketleftidx,
          widthrelative2bucketleftidx,
          *ulongtab;      /* clone */
  uint32_t *uinttab;      /* clone */
};

static bool gt_decide_to_use_uint(bool useuint,GtUword maxvalue)
{
  if (useuint && maxvalue <= (GtUword) UINT_MAX)
  {
    return true;
  }
  return false;
}

uint64_t gt_suffixsortspace_requiredspace(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint)
{
  uint64_t requiredspace = (uint64_t) sizeof (GtSuffixsortspace);

  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    gt_assert(maxvalue <= (GtUword) UINT_MAX);
    requiredspace += (uint64_t) numofentries * (uint64_t) sizeof (uint32_t);
  } else
  {
    requiredspace += (uint64_t) numofentries *
                     (uint64_t) sizeof (GtUword);
  }
  return requiredspace;
}

static void gt_suffixsortspace_overflow_abort(GT_UNUSED const char *f,
                                              GT_UNUSED int l,
                                              void *data)
{
  fprintf(stderr, "error: overflow detected while calculating size of "
                  "suffix sorting space: "GT_WU" * "GT_WU" bytes is too large "
                  "for " "the current platform, please recompile GenomeTools "
                  "with support for a larger address space to prevent this "
                  "(e.g. 64 bit instead of 32 bit) or use the `-parts' "
                  "option.\n",
                  (GtUword) sizeof (GtUword),
                  *(GtUword*) data);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

static GtSuffixsortspace *gt_suffixsortspace_new_generic(GtUword numofentries,
                                                  GtUword maxvalue,
                                                  bool useuint,
                                                  void *tab2clone,
                                                  GT_UNUSED GtLogger *logger)
{
  GtSuffixsortspace *suffixsortspace;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
  suffixsortspace->maxindex = numofentries-1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
  suffixsortspace->exportptr.ulongtabsectionptr = NULL;
  suffixsortspace->exportptr.uinttabsectionptr = NULL;
  suffixsortspace->currentexport = false;
  suffixsortspace->partoffset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->widthrelative2bucketleftidx = 0;
#if defined (_LP64) || defined (_WIN64)
  gt_logger_log(logger,"%s suffix_sort_space: suftab uses %dbit values: "
                       "maxvalue="GT_WU",numofentries="GT_WU,
                       tab2clone == NULL ? "create" : "clone",
                       gt_decide_to_use_uint(useuint,maxvalue) ? 32 : 64,
                       maxvalue,numofentries);
#endif
  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    if (tab2clone == NULL)
    {
      GtUword sufspacesize
        = gt_safe_mult_ulong_check((GtUword)
                                   sizeof (*suffixsortspace->uinttab),
                                   numofentries,
                                   gt_suffixsortspace_overflow_abort,
                                   &numofentries);
      suffixsortspace->uinttab = gt_malloc((size_t) sufspacesize);
      suffixsortspace->clonenumber = 0;
    } else
    {
      suffixsortspace->uinttab = (uint32_t *) tab2clone;
      suffixsortspace->clonenumber = UINT_MAX; /* undefined */
    }
    suffixsortspace->ulongtab = NULL;
  } else
  {
    if (tab2clone == NULL)
    {
      GtUword sufspacesize
        = gt_safe_mult_ulong_check((GtUword)
                                   sizeof (*suffixsortspace->ulongtab),
                                   numofentries,
                                   gt_suffixsortspace_overflow_abort,
                                   &numofentries);
      suffixsortspace->ulongtab = gt_malloc((size_t) sufspacesize);
      suffixsortspace->clonenumber = 0;
    } else
    {
      suffixsortspace->ulongtab = (GtUword *) tab2clone;
      suffixsortspace->clonenumber = UINT_MAX; /* undefined */
    }
    suffixsortspace->uinttab = NULL;
  }
  return suffixsortspace;
}

GtSuffixsortspace *gt_suffixsortspace_new(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint,
                                          GtLogger *logger)
{
  return gt_suffixsortspace_new_generic(numofentries,
                                        maxvalue,
                                        useuint,
                                        NULL,
                                        logger);
}

GtSuffixsortspace *gt_suffixsortspace_clone(GtSuffixsortspace *sssp,
                                            unsigned int clonenumber,
                                            GtLogger *logger)
{
  bool useuint;
  void *tab2clone;
  GtSuffixsortspace *cloned_sssp;

  gt_assert(sssp != NULL);
  useuint = sssp->uinttab != NULL ? true : false;
  tab2clone = useuint ? (void *) sssp->uinttab
                      : (void *) sssp->ulongtab;
  gt_assert(tab2clone != NULL);
  cloned_sssp = gt_suffixsortspace_new_generic(sssp->maxindex + 1,
                                               sssp->maxvalue,
                                               useuint,
                                               tab2clone,
                                               logger);
  gt_assert(clonenumber > 0);
  cloned_sssp->clonenumber = clonenumber;
  cloned_sssp->partoffset = sssp->partoffset;
  return cloned_sssp;
}

void gt_suffixsortspace_delete(GtSuffixsortspace *suffixsortspace,
                               bool checklongestdefined)
{
  if (suffixsortspace != NULL)
  {
    if (checklongestdefined && !suffixsortspace->longestidx.defined)
    {
      fprintf(stderr,"%s, l. %d: longest is not defined\n",__FILE__,__LINE__);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (suffixsortspace->clonenumber == 0)
    {
      gt_free(suffixsortspace->uinttab);
      gt_free(suffixsortspace->ulongtab);
    }
    gt_free(suffixsortspace);
  }
}

void gt_suffixsortspace_nooffsets(GT_UNUSED const GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL && sssp->partoffset == 0 && sssp->bucketleftidx == 0);
}

GtUword gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,GtUword idx)
{
  gt_assert(sssp != NULL);
  if (idx > sssp->maxindex)
  {
    printf("idx = " GT_WU " > " GT_WU " = maxindex\n",idx,sssp->maxindex);
  }
  gt_assert(idx <= sssp->maxindex);
  if (sssp->ulongtab != NULL)
  {
    return sssp->ulongtab[idx];
  }
  gt_assert(sssp->uinttab != NULL);
  return (GtUword) sssp->uinttab[idx];
}

void gt_suffixsortspace_updatelongest(GtSuffixsortspace *sssp,GtUword idx)
{
  gt_assert(sssp != NULL);

  sssp->longestidx.defined = true;
  sssp->longestidx.valueunsignedlong = sssp->bucketleftidx + idx;
}

static void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                         GtUword idx,
                                         GtUword value)
{
  gt_assert(sssp != NULL && idx <= sssp->maxindex && value <= sssp->maxvalue);
  if (value == 0)
  {
    sssp->longestidx.defined = true;
    sssp->longestidx.valueunsignedlong = sssp->partoffset + idx;
  }
  if (sssp->ulongtab != NULL)
  {
    sssp->ulongtab[idx] = value;
  } else
  {
    gt_assert (sssp->uinttab != NULL && value <= (GtUword) UINT_MAX);
    sssp->uinttab[idx] = (uint32_t) value;
  }
}

void gt_suffixsortspace_showrange(const GtSuffixsortspace *sssp,
                                  GtUword subbucketleft,
                                  GtUword width)
{
  GtUword idx;

  gt_assert(sssp != NULL);
  printf(GT_WU","GT_WU"=",sssp->bucketleftidx+subbucketleft-sssp->partoffset,
                    sssp->bucketleftidx+subbucketleft+width-1-sssp->partoffset);
  for (idx=sssp->bucketleftidx+subbucketleft-sssp->partoffset;
       idx<sssp->bucketleftidx+subbucketleft+width-sssp->partoffset;
       idx++)
  {
    printf(" "GT_WU, gt_suffixsortspace_getdirect(sssp,idx));
  }
}

GtSuffixsortspace_exportptr *gt_suffixsortspace_exportptr(
                                  GtSuffixsortspace *sssp,
                                  GtUword subbucketleft)
{
  gt_assert(sssp != NULL);
  if (sssp->ulongtab != NULL)
  {
    sssp->exportptr.ulongtabsectionptr = sssp->ulongtab + sssp->bucketleftidx
                                                        + subbucketleft
                                                        - sssp->partoffset;
    sssp->exportptr.uinttabsectionptr = NULL;
  } else
  {
    gt_assert(sssp->uinttab != NULL);
    sssp->exportptr.uinttabsectionptr = sssp->uinttab + sssp->bucketleftidx
                                                      + subbucketleft
                                                      - sssp->partoffset;
    sssp->exportptr.ulongtabsectionptr = NULL;
  }
  sssp->currentexport = true;
  return &sssp->exportptr;
}

void gt_suffixsortspace_export_done (GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL);
  sssp->exportptr.ulongtabsectionptr = NULL;
  sssp->exportptr.uinttabsectionptr = NULL;
  sssp->currentexport = false;
}

GtUword gt_suffixsortspace_get (const GtSuffixsortspace *sssp,
                                GtUword subbucketleft,
                                GtUword idx)
{
  gt_assert(sssp != NULL);
  return gt_suffixsortspace_getdirect(sssp, sssp->bucketleftidx
                                              + subbucketleft
                                              + idx
                                              - sssp->partoffset);
}

const GtUword *gt_suffixsortspace_getptr_ulong (const GtSuffixsortspace *sssp,
                                                GtUword subbucketleft)
{
  gt_assert(sssp != NULL);
  if (sssp->ulongtab != NULL)
  {
    return sssp->ulongtab + sssp->bucketleftidx + subbucketleft
                          - sssp->partoffset;
  } else
  {
    return NULL;
  }
}

const uint32_t *gt_suffixsortspace_getptr_uint32 (const GtSuffixsortspace *sssp,
                                                  GtUword subbucketleft)
{
  gt_assert(sssp != NULL);
  if (sssp->uinttab != NULL)
  {
    return sssp->uinttab + sssp->bucketleftidx + subbucketleft
                         - sssp->partoffset;
  } else
  {
    return NULL;
  }
}

void gt_suffixsortspace_set (GtSuffixsortspace *sssp,
                             GtUword subbucketleft,
                             GtUword idx,
                             GtUword value)
{
  GtUword updateindex;

  gt_assert(sssp != NULL);
  updateindex = sssp->bucketleftidx + subbucketleft + idx - sssp->partoffset;
  gt_assert (sssp->widthrelative2bucketleftidx  == 0 ||
             updateindex <
             sssp->bucketleftidx + sssp->widthrelative2bucketleftidx);
  gt_suffixsortspace_setdirect(sssp, updateindex,value);
}

void gt_suffixsortspace_init_seqstartpos(GtSuffixsortspace *sssp,
                                         const GtEncseq *encseq)
{
  GtUword idx, numofsequences = gt_encseq_num_of_sequences(encseq);

  for (idx = 0; idx < numofsequences; idx++)
  {
    gt_suffixsortspace_setdirect(sssp, idx, gt_encseq_seqstartpos(encseq, idx));
  }
}

void gt_suffixsortspace_init_identity(GtSuffixsortspace *sssp,
                                      GtUword numofsuffixes)
{
  GtUword idx;

  for (idx=0; idx<numofsuffixes; idx++)
  {
    gt_suffixsortspace_setdirect(sssp,idx,idx);
  }
}

GtUword gt_suffixsortspace_bucketleftidx_get (const GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL);
  return sssp->bucketleftidx;
}

void gt_suffixsortspace_bucketrange_set(GtSuffixsortspace *sssp,
                                        GtUword bucketleftidx,
                                        GtUword widthrelative2bucketleftidx)
{
  gt_assert(sssp != NULL && (sssp->bucketleftidx == bucketleftidx ||
                             !sssp->currentexport));
  sssp->bucketleftidx = bucketleftidx;
  sssp->widthrelative2bucketleftidx = widthrelative2bucketleftidx;
}

void gt_suffixsortspace_bucketrange_reset(GtSuffixsortspace *sssp)
{
  sssp->bucketleftidx = 0;
  sssp->widthrelative2bucketleftidx = 0;
}

void gt_suffixsortspace_partoffset_set (GtSuffixsortspace *sssp,
                                        GtUword partoffset)
{
  gt_assert(sssp != NULL && (sssp->partoffset == partoffset ||
                             !sssp->currentexport));
  sssp->partoffset = partoffset;
}

void gt_suffixsortspace_sortspace_delete (GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL);
  gt_free(sssp->ulongtab);
  sssp->ulongtab = NULL;
  gt_free(sssp->uinttab);
  sssp->uinttab = NULL;
}

void gt_suffixsortspace_delete_cloned(GtSuffixsortspace **sssp_tab,
                                      unsigned int parts)
{
  bool found = false;
  unsigned int p;

  gt_assert(sssp_tab != NULL && parts > 1U && sssp_tab[0]->clonenumber == 0);
  for (p = 1U; p < parts; p++)
  {
    GtSuffixsortspace *sssp = sssp_tab[p];

    gt_assert(sssp->clonenumber == p);
    if (!found && sssp->longestidx.defined)
    {
      GtUword zerotest;
      gt_assert(sssp->longestidx.valueunsignedlong >= sssp->partoffset);
      zerotest = sssp->longestidx.valueunsignedlong - sssp->partoffset;
      if ((sssp->ulongtab != NULL && sssp->ulongtab[zerotest] == 0) ||
          (sssp->uinttab != NULL && sssp->uinttab[zerotest] == 0))
      {
        sssp_tab[0]->longestidx.defined = true;
        sssp_tab[0]->longestidx.valueunsignedlong
          = sssp->longestidx.valueunsignedlong;
        found = true;
      }
    }
    gt_suffixsortspace_delete(sssp_tab[p],false);
    sssp_tab[p] = NULL;
  }
}

const GtUword *gt_suffixsortspace_ulong_get (const GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL && sssp->ulongtab != NULL);
  return sssp->ulongtab;
}

GtUword gt_suffixsortspace_longest(const GtSuffixsortspace *sssp)
{
  gt_assert(sssp != NULL && sssp->longestidx.defined);
  return sssp->longestidx.valueunsignedlong;
}

void gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                 const GtSuffixsortspace *sssp,
                                 GtUword numberofsuffixes)
{
  gt_assert(sssp != NULL);
  if (sssp->ulongtab != NULL)
  {
    gt_xfwrite((void *) sssp->ulongtab,sizeof (*sssp->ulongtab),
               (size_t) numberofsuffixes, outfpsuftab);
  } else
  {
    gt_assert(sssp->uinttab != NULL);
    gt_xfwrite((void *) sssp->uinttab,sizeof (*sssp->uinttab),
               (size_t) numberofsuffixes, outfpsuftab);
  }
}

void gt_suffixsortspace_compressed_to_file (const GtSuffixsortspace *sssp,
                                            GtBitbuffer *bb,
                                            GtUword numberofsuffixes)
{
  gt_assert(sssp != NULL);
  if (sssp->ulongtab != NULL)
  {
    gt_bitbuffer_next_ulongtab(bb,sssp->ulongtab,numberofsuffixes);
  } else
  {
    gt_assert (sssp->uinttab != NULL);
    gt_bitbuffer_next_uint32tab(bb,sssp->uinttab,numberofsuffixes);
  }
}

static GtUword gt_suffixsortspace_insertfullspecialrange(
                                                  GtSuffixsortspace *sssp,
                                                  GtReadmode readmode,
                                                  GtUword nextfree,
                                                  GtUword totallength,
                                                  GtUword leftpos,
                                                  GtUword rightpos)
{
  GtUword pos;

  gt_assert(leftpos < rightpos);
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (GT_ISDIRREVERSE(readmode))
    {
      gt_assert(pos < totallength);
      gt_suffixsortspace_setdirect(sssp,nextfree,
                                   GT_REVERSEPOS(totallength,pos));
      nextfree++;
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      gt_suffixsortspace_setdirect(sssp, nextfree, pos);
      nextfree++;
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
  return nextfree;
}

struct GtSSSPbuf
{
  GtRange overhang;
  GtUword size, nextfree;
};

GtSSSPbuf *gt_SSSPbuf_new(GtUword size)
{
  GtSSSPbuf *sssp_buf = gt_malloc(sizeof *sssp_buf);

  sssp_buf->size = size;
  sssp_buf->nextfree = 0;
  sssp_buf->overhang.start = sssp_buf->overhang.end = 0;
  return sssp_buf;
}

GtUword gt_SSSPbuf_filled(const GtSSSPbuf *sssp_buf)
{
  gt_assert(sssp_buf != NULL);
  return sssp_buf->nextfree;
}

void gt_SSSPbuf_delete(GtSSSPbuf *sssp_buf)
{
  if (sssp_buf != NULL)
  {
    gt_free(sssp_buf);
  }
}

static void gt_SSSPbuf_insertfullspecialrange(GtSuffixsortspace *sssp,
                                              GtReadmode readmode,
                                              GtSSSPbuf *sssp_buf,
                                              GtUword totallength,
                                              GtUword leftpos,
                                              GtUword rightpos)
{
  gt_assert(sssp_buf != NULL);
  sssp_buf->nextfree
    = gt_suffixsortspace_insertfullspecialrange(sssp,
                                                readmode,
                                                sssp_buf->nextfree,
                                                totallength,
                                                leftpos,
                                                rightpos);
}

bool gt_SSSPbuf_fillspecialnextpage(GtSuffixsortspace *sssp,
                                    GtReadmode readmode,
                                    GtSpecialrangeiterator *sri,
                                    GtUword totallength,
                                    GtSSSPbuf *sssp_buf)
{
  bool exhausted = false;

  gt_assert(sssp_buf != NULL);
  sssp_buf->nextfree = 0;
  while (true)
  {
    if (sssp_buf->overhang.start < sssp_buf->overhang.end)
    {
      GtUword width = sssp_buf->overhang.end - sssp_buf->overhang.start;
      if (sssp_buf->nextfree + width > sssp_buf->size)
      {
        /* does not fit into the buffer, so only output a part */
        GtUword rest = sssp_buf->nextfree + width - sssp_buf->size;
        gt_assert(rest > 0);
        if (GT_ISDIRREVERSE(readmode))
        {
          gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                            sssp_buf->overhang.start + rest,
                                            sssp_buf->overhang.end);
          sssp_buf->overhang.end = sssp_buf->overhang.start + rest;
        } else
        {
          gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                            sssp_buf->overhang.start,
                                            sssp_buf->overhang.end - rest);
          sssp_buf->overhang.start = sssp_buf->overhang.end - rest;
        }
        break;
      }
      if (sssp_buf->nextfree + width == sssp_buf->size)
      { /* overhang fits into the buffer and buffer is full */
        gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                          sssp_buf->overhang.start,
                                          sssp_buf->overhang.end);
        sssp_buf->overhang.start = sssp_buf->overhang.end = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                        sssp_buf->overhang.start,
                                        sssp_buf->overhang.end);
      sssp_buf->overhang.start = sssp_buf->overhang.end = 0;
    } else
    {
      GtRange range;

      if (sri != NULL && gt_specialrangeiterator_next(sri,&range))
      {
        GtUword width = range.end - range.start;
        gt_assert(width > 0);
        if (sssp_buf->nextfree + width > sssp_buf->size)
        { /* does not fit into the buffer, so only output a part */
          GtUword rest = sssp_buf->nextfree + width - sssp_buf->size;
          if (GT_ISDIRREVERSE(readmode))
          {
            gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,
                                              totallength,
                                              range.start + rest,
                                              range.end);
            sssp_buf->overhang.start = range.start;
            sssp_buf->overhang.end = range.start + rest;
          } else
          {
            gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,
                                              totallength,
                                              range.start,
                                              range.end - rest);
            sssp_buf->overhang.start = range.end - rest;
            sssp_buf->overhang.end = range.end;
          }
          break;
        }
        if (sssp_buf->nextfree + width == sssp_buf->size)
        { /* overhang fits into the buffer and buffer is full */
          gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                            range.start,range.end);
          sssp_buf->overhang.start = sssp_buf->overhang.end = 0;
          break;
        }
        /* overhang fits into the buffer and buffer is not full */
        gt_SSSPbuf_insertfullspecialrange(sssp,readmode,sssp_buf,totallength,
                                          range.start,range.end);
        sssp_buf->overhang.start = sssp_buf->overhang.end = 0;
      } else
      {
        if (sssp_buf->nextfree < sssp_buf->size)
        {
          /* add last suffix for start position totallength */
          gt_suffixsortspace_setdirect(sssp,sssp_buf->nextfree,totallength);
          sssp_buf->nextfree++;
          exhausted = true;
        }
        break;
      }
    }
  }
  gt_assert(sssp_buf->nextfree > 0);
  return exhausted;
}
