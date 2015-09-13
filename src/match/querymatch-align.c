/*
  Copyright (c) 2007-2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "extended/alignment.h"
#include "extended/linearalign.h"
#undef SKDEBUG
#ifdef SKDEBUG
#include "match/greedyedist.h"
#endif
#include "querymatch-align.h"
#include "seed-extend.h"

static void eoplist_reverse_order(uint8_t *start,uint8_t *end)
{
  if (start < end)
  {
    uint8_t *fwd, *bwd;

    for (fwd = start, bwd = end; fwd < bwd; fwd++, bwd--)
    {
      uint8_t tmp = *fwd;
      *fwd = *bwd;
      *bwd = tmp;
    }
  }
}

static void converteoplist2alignment(GtAlignment *alignment,
                                     const GtArrayuint8_t *eoplist)
{
  uint8_t *ptr;

  for (ptr = eoplist->spaceuint8_t;
       ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t; ptr++)
  {
    if (*ptr == FT_EOPCODE_DELETION)
    {
      gt_alignment_add_deletion(alignment);
    } else
    {
      if (*ptr == FT_EOPCODE_INSERTION)
      {
        gt_alignment_add_insertion(alignment);
      } else
      {
        gt_alignment_add_replacement_multi(alignment,(GtUword) ((*ptr) + 1));
      }
    }
  }
}

void gt_querymatch_showeoplist(const GtArrayuint8_t *eoplist)
{
  uint8_t *ptr;

  gt_assert(eoplist != NULL);
  if (eoplist->nextfreeuint8_t == 0)
  {
    printf("[]\n");
    return;
  }
  for (ptr = eoplist->spaceuint8_t;
       ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t;
       ptr++)
  {
    if (*ptr == FT_EOP_DELETION)
    {
      printf("D");
    } else
    {
      if (*ptr == FT_EOP_INSERTION)
      {
        printf("I");
      } else
      {
        printf("R %d",(int) ((*ptr) + 1));
      }
    }
    if (ptr < eoplist->spaceuint8_t + eoplist->nextfreeuint8_t - 1)
    {
      printf(",");
    } else
    {
      printf("\n");
    }
  }
}

struct GtQuerymatchoutoptions
{
  GtUword totallength,
          alignmentwidth, /* if > 0, then with alignment of this width */
          useqbuffer_size,
          vseqbuffer_size;
  GtArrayuint8_t eoplist;
  Fronttrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  GtUchar *useqbuffer, *vseqbuffer, *alignment_show_buffer;
  const GtUchar *characters;
  GtAlignment *alignment;
  GtEncseqReader *esr_for_align_show;
  GtUchar wildcardshow;
  GtSeqpaircoordinates correction_info;
};

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(
                                GtUword alignmentwidth,
                                GtUword errorpercentage,
                                GtUword maxalignedlendifference,
                                GtUword history,
                                GtUword perc_mat_history,
                                GtExtendCharAccess extend_char_access,
                                GtUword sensitivity)
{
  GtQuerymatchoutoptions *querymatchoutoptions
    = gt_malloc(sizeof *querymatchoutoptions);

  querymatchoutoptions->alignmentwidth = alignmentwidth;
  querymatchoutoptions->front_trace = NULL;
  querymatchoutoptions->ggemi = NULL;
  querymatchoutoptions->useqbuffer = NULL;
  querymatchoutoptions->useqbuffer_size = 0;
  querymatchoutoptions->vseqbuffer = NULL;
  querymatchoutoptions->vseqbuffer_size = 0;
  GT_INITARRAY(&querymatchoutoptions->eoplist,uint8_t);
  querymatchoutoptions->totallength = GT_UWORD_MAX; /* not yet known */
  querymatchoutoptions->esr_for_align_show = NULL;
  querymatchoutoptions->characters = NULL;
  querymatchoutoptions->alignment = NULL;
  if (alignmentwidth > 0)
  {
    querymatchoutoptions->alignment_show_buffer
      = gt_alignment_buffer_new(alignmentwidth);
  } else
  {
    querymatchoutoptions->alignment_show_buffer = NULL;
  }
  if (errorpercentage > 0)
  {
    if (alignmentwidth > 0)
    {
      querymatchoutoptions->front_trace = front_trace_new();
      querymatchoutoptions->alignment = gt_alignment_new();
    }
    querymatchoutoptions->ggemi
      = gt_greedy_extend_matchinfo_new(errorpercentage,
                                       maxalignedlendifference,
                                       history, /* default value */
                                       perc_mat_history,
                                       0,/* userdefinedleastlength not used */
                                       extend_char_access,
                                       sensitivity);
  }
  return querymatchoutoptions;
}

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    front_trace_delete(querymatchoutoptions->front_trace);
    gt_greedy_extend_matchinfo_delete(querymatchoutoptions->ggemi);
    gt_free(querymatchoutoptions->useqbuffer);
    if (querymatchoutoptions->vseqbuffer_size > 0)
    {
      gt_free(querymatchoutoptions->vseqbuffer);
    }
    gt_alignment_delete(querymatchoutoptions->alignment);
    GT_FREEARRAY(&querymatchoutoptions->eoplist,uint8_t);
    gt_alignment_buffer_delete(querymatchoutoptions->alignment_show_buffer);
    gt_encseq_reader_delete(querymatchoutoptions->esr_for_align_show);
    gt_free(querymatchoutoptions);
  }
}

static bool seededmatch2eoplist(GtQuerymatchoutoptions *querymatchoutoptions,
                                const GtEncseq *encseq,
                                const GtUchar *query,
                                GtUword query_totallength,
                                GtUword dbstart,
                                GtUword dblen,
                                GtUword querystartabsolute,
                                GtUword querylen,
                                GtUword seedpos1,
                                GtUword seedpos2,
                                GtUword seedlen,
                                bool greedyextension)
{
  GtUword ulen, vlen, ustart, vstart;
  bool alignment_succeeded = true;
  Polished_point right_best_polished_point = {0,0,0};
  Polished_point left_best_polished_point = {0,0,0};
  GtSeqpaircoordinates *coords;

  gt_assert(querymatchoutoptions != NULL);
  querymatchoutoptions->eoplist.nextfreeuint8_t = 0;
  if (querymatchoutoptions->totallength == GT_UWORD_MAX)
  {
    querymatchoutoptions->totallength = gt_encseq_total_length(encseq);
  }
  ustart = seedpos1 + seedlen;
  vstart = seedpos2 + seedlen;
  gt_assert(dbstart + dblen >= ustart);
  ulen = dbstart + dblen - ustart;
  gt_assert(querystartabsolute + querylen >= vstart);
  vlen = querystartabsolute + querylen - vstart;
  if (ulen > 0 && vlen > 0)
  {
    if (gt_align_front_prune_edist(true,
                                   &right_best_polished_point,
                                   querymatchoutoptions->front_trace,
                                   encseq,
                                   query,
                                   query_totallength,
                                   querymatchoutoptions->ggemi,
                                   greedyextension,
                                   ustart,
                                   ulen,
                                   vstart,
                                   vlen) == ulen + vlen + 1)
    {
      alignment_succeeded = false;
    } else
    {
      if (querymatchoutoptions->front_trace != NULL)
      {
        front_trace2eoplist(&querymatchoutoptions->eoplist,
                            querymatchoutoptions->front_trace,
                            &right_best_polished_point,
                            ulen,
                            vlen);
      }
    }
    if (querymatchoutoptions->front_trace != NULL)
    {
      gt_assert(querymatchoutoptions->alignmentwidth > 0);
      front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
    }
  }
  if (alignment_succeeded)
  {
    if (querymatchoutoptions->alignmentwidth > 0)
    {
      front_trace_multireplacement(&querymatchoutoptions->eoplist,seedlen);
    }
    if (seedpos1 > dbstart && seedpos2 > querystartabsolute)
    {
      ustart = GT_REVERSEPOS(querymatchoutoptions->totallength,seedpos1 - 1);
      ulen = seedpos1 - dbstart;
      if (query == NULL)
      {
        vstart = GT_REVERSEPOS(querymatchoutoptions->totallength,seedpos2 - 1);
      } else
      {
        vstart = seedpos2 - 1;
      }
      vlen = seedpos2 - querystartabsolute;
      if (gt_align_front_prune_edist(false,
                                     &left_best_polished_point,
                                     querymatchoutoptions->front_trace,
                                     encseq,
                                     query,
                                     query_totallength,
                                     querymatchoutoptions->ggemi,
                                     greedyextension,
                                     ustart,
                                     ulen,
                                     vstart,
                                     vlen) == ulen + vlen + 1)
      {
        alignment_succeeded = false;
      } else
      {
        if (querymatchoutoptions->front_trace != NULL)
        {
          GtUword eoplistlen = querymatchoutoptions->eoplist.nextfreeuint8_t;

          front_trace2eoplist(&querymatchoutoptions->eoplist,
                              querymatchoutoptions->front_trace,
                              &left_best_polished_point,
                              ulen,
                              vlen);
          eoplist_reverse_order(querymatchoutoptions->eoplist.spaceuint8_t +
                                 eoplistlen,
                                querymatchoutoptions->eoplist.spaceuint8_t +
                                 querymatchoutoptions->eoplist.nextfreeuint8_t -
                                 1);
        }
      }
      if (querymatchoutoptions->front_trace != NULL)
      {
        front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
      }
    }
  }
  coords = &querymatchoutoptions->correction_info;
  if (alignment_succeeded)
  {
    GtUword leftcolumn, rightcolumn;

    gt_assert(dbstart <= seedpos1 - left_best_polished_point.row);
    coords->uoffset = seedpos1 - left_best_polished_point.row - dbstart;
    coords->ulen = seedlen + left_best_polished_point.row +
                   right_best_polished_point.row;
    leftcolumn = left_best_polished_point.alignedlen -
                 left_best_polished_point.row;
    rightcolumn = right_best_polished_point.alignedlen -
                  right_best_polished_point.row;
    gt_assert(seedpos2 >= leftcolumn &&
              querystartabsolute <= seedpos2 - leftcolumn);
    coords->voffset = seedpos2 - leftcolumn - querystartabsolute;
    coords->vlen = seedlen + leftcolumn + rightcolumn;
    coords->sumdist = left_best_polished_point.distance +
                      right_best_polished_point.distance;
  } else
  {
    coords->uoffset = coords->voffset = 0;
    coords->ulen = dblen;
    coords->vlen = querylen;
    coords->sumdist = 0;
  }
  return alignment_succeeded;
}

#ifdef SKDEBUG
static void check_correct_edist(const GtUchar *useq,
                                GtUword ulen,
                                const GtUchar *vseq,
                                GtUword vlen,
                                GtUword edist)
{
  GtUword realedist;
  GtFrontResource *ftres = gt_frontresource_new(edist);
  GtSeqabstract *useq_abstract = gt_seqabstract_new_gtuchar(useq, ulen, 0);
  GtSeqabstract *vseq_abstract = gt_seqabstract_new_gtuchar(vseq, vlen, 0);

  realedist = greedyunitedist(ftres,useq_abstract,vseq_abstract);
  gt_assert(realedist <= edist);
  gt_seqabstract_delete(useq_abstract);
  gt_seqabstract_delete(vseq_abstract);
  gt_frontresource_delete(ftres);
}
#endif

bool gt_querymatchoutoptions_alignment_prepare(GtQuerymatchoutoptions
                                                *querymatchoutoptions,
                                               const GtEncseq *encseq,
                                               const GtUchar *query,
                                               GtUword query_totallength,
                                               GtUword dbstart,
                                               GtUword dblen,
                                               GtUword querystartabsolute,
                                               GtUword querylen,
                                               GtUword edist,
                                               GtUword seedpos1,
                                               GtUword seedpos2,
                                               GtUword seedlen,
                                               bool greedyextension)
{
  bool seededalignment = false;

  gt_assert(querymatchoutoptions != NULL);
  if (querymatchoutoptions->alignmentwidth > 0 &&
      querymatchoutoptions->characters == NULL)
  {
    querymatchoutoptions->characters
      = gt_encseq_alphabetcharacters(encseq);
    querymatchoutoptions->wildcardshow
      = gt_alphabet_wildcard_show(gt_encseq_alphabet(encseq));
  }
  if (querymatchoutoptions->esr_for_align_show == NULL)
  {
    querymatchoutoptions->esr_for_align_show
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (dblen > querymatchoutoptions->useqbuffer_size)
  {
    querymatchoutoptions->useqbuffer
      = gt_realloc(querymatchoutoptions->useqbuffer,
                   sizeof *querymatchoutoptions->useqbuffer * dblen);
    querymatchoutoptions->useqbuffer_size = dblen;
  }
  gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->esr_for_align_show,
                            encseq,
                            querymatchoutoptions->useqbuffer,
                            dbstart,
                            dbstart + dblen - 1);
  if (edist > 0)
  {
    if (query == NULL)
    {
      if (querylen > querymatchoutoptions->vseqbuffer_size)
      {
        querymatchoutoptions->vseqbuffer
          = gt_realloc(querymatchoutoptions->vseqbuffer,
                       sizeof *querymatchoutoptions->vseqbuffer * querylen);
        querymatchoutoptions->vseqbuffer_size = querylen;
      }
      gt_encseq_extract_encoded_with_reader(
                              querymatchoutoptions->esr_for_align_show,
                              encseq,
                              querymatchoutoptions->vseqbuffer,
                              querystartabsolute,
                              querystartabsolute + querylen - 1);
    } else
    {
      querymatchoutoptions->vseqbuffer
        = (GtUchar *) (query + querystartabsolute);
    }
    if (seededmatch2eoplist(querymatchoutoptions,
                            encseq,
                            query,
                            query_totallength,
                            dbstart,
                            dblen,
                            querystartabsolute,
                            querylen,
                            seedpos1,
                            seedpos2,
                            seedlen,
                            greedyextension))
    {
      if (querymatchoutoptions->alignmentwidth > 0)
      {
        gt_alignment_set_seqs(querymatchoutoptions->alignment,
                              querymatchoutoptions->useqbuffer +
                                  querymatchoutoptions->correction_info.uoffset,
                              querymatchoutoptions->correction_info.ulen,
                              querymatchoutoptions->vseqbuffer +
                                  querymatchoutoptions->correction_info.voffset,
                              querymatchoutoptions->correction_info.vlen);
        converteoplist2alignment(querymatchoutoptions->alignment,
                                 &querymatchoutoptions->eoplist);
      }
      seededalignment = true;
    } else
    {
      if (querymatchoutoptions->alignmentwidth > 0)
      {
#ifndef NDEBUG
        GtUword linedist;
#endif
        gt_assert(!greedyextension);
        gt_alignment_set_seqs(querymatchoutoptions->alignment,
                              querymatchoutoptions->useqbuffer,
                              dblen,
                              querymatchoutoptions->vseqbuffer,
                              querylen);
#ifndef NDEBUG
        linedist =
#else
        (void)
#endif
                   gt_computelinearspace(querymatchoutoptions->alignment,
                                         querymatchoutoptions->useqbuffer,
                                         0,
                                         dblen,
                                         querymatchoutoptions->vseqbuffer,
                                         0,
                                         querylen,
                                         0,
                                         1,
                                         1);
        gt_assert(linedist <= edist);
        /*printf("linedist = " GT_WU " <= " GT_WU "= edist\n",linedist,edist);*/
      }
    }
#ifdef SKDEBUG
    check_correct_edist(querymatchoutoptions->useqbuffer,
                        dblen,
                        querymatchoutoptions->vseqbuffer,
                        querylen,
                        edist);
#endif
  }
  return seededalignment;
}

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions,
                                            GtUword edist,
                                            GtUword len)
{
  if (querymatchoutoptions != NULL && querymatchoutoptions->alignmentwidth > 0)
  {
    if (edist > 0)
    {
      gt_alignment_show_generic(querymatchoutoptions->alignment_show_buffer,
                                querymatchoutoptions->alignment,
                                stdout,
                                (unsigned int)
                                querymatchoutoptions->alignmentwidth,
                                querymatchoutoptions->characters,
                                querymatchoutoptions->wildcardshow);
      gt_alignment_reset(querymatchoutoptions->alignment);
    } else
    {
      gt_alignment_exact_show(querymatchoutoptions->alignment_show_buffer,
                              querymatchoutoptions->useqbuffer,
                              len,
                              stdout,
                              querymatchoutoptions->alignmentwidth,
                              querymatchoutoptions->characters);
    }
  }
}

const GtSeqpaircoordinates *gt_querymatchoutoptions_correction_get(
              const GtQuerymatchoutoptions *querymatchoutoptions)
{
  return &querymatchoutoptions->correction_info;
}
