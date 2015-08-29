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
  Fronttrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  GtUword seedpos1, seedpos2, seedlen, totallength,
          alignmentwidth, /* if > 0, then with alignment of this width */
          useqbuffer_size,
          vseqbuffer_size;
  GtUchar *useqbuffer, *vseqbuffer;
  GtAlignment *alignment;
  GtArrayuint8_t eoplist;
  GtUchar *alignment_show_buffer;
  GtEncseqReader *esr_for_align_show;
  const GtUchar *characters;
  GtUchar wildcardshow;
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
  querymatchoutoptions->alignment_show_buffer = NULL;
  querymatchoutoptions->front_trace = NULL;
  querymatchoutoptions->ggemi = NULL;
  querymatchoutoptions->useqbuffer = NULL;
  querymatchoutoptions->useqbuffer_size = 0;
  querymatchoutoptions->vseqbuffer = NULL;
  querymatchoutoptions->vseqbuffer_size = 0;
  querymatchoutoptions->alignment = NULL;
  GT_INITARRAY(&querymatchoutoptions->eoplist,uint8_t);
  querymatchoutoptions->totallength = GT_UWORD_MAX; /* not yet known */
  querymatchoutoptions->esr_for_align_show = NULL;
  querymatchoutoptions->characters = NULL;
  gt_assert(alignmentwidth > 0);
  if (alignmentwidth > 0)
  {
    querymatchoutoptions->alignment_show_buffer
      = gt_alignment_buffer_new(alignmentwidth);
    if (errorpercentage > 0)
    {
      querymatchoutoptions->front_trace = front_trace_new();
      querymatchoutoptions->ggemi
        = gt_greedy_extend_matchinfo_new(errorpercentage,
                                         maxalignedlendifference,
                                         history, /* default value */
                                         perc_mat_history,
                                         0,/* userdefinedleastlength not used */
                                         extend_char_access,
                                         sensitivity);
    }
    querymatchoutoptions->alignment = gt_alignment_new();
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
    gt_free(querymatchoutoptions->vseqbuffer);
    gt_alignment_delete(querymatchoutoptions->alignment);
    GT_FREEARRAY(&querymatchoutoptions->eoplist,uint8_t);
    gt_alignment_buffer_delete(querymatchoutoptions->alignment_show_buffer);
    gt_encseq_reader_delete(querymatchoutoptions->esr_for_align_show);
    gt_free(querymatchoutoptions);
  }
}

typedef struct
{
  GtUword uoffset, voffset, ulen, vlen;
} GtSeqpaircoordinates;

static bool seededmatch2eoplist(GtSeqpaircoordinates *coords,
                                GtUword dbstart,
                                GtUword dblen,
                                GtUword querystartabsolute,
                                GtUword querylen,
                                GtQuerymatchoutoptions *querymatchoutoptions,
                                const GtEncseq *encseq)
{
  GtUword ulen, vlen, rightdistance = 0, ustart, vstart;
  bool alignment_succeeded = true;
  Polished_point right_best_polished_point = {0,0,0};
  Polished_point left_best_polished_point = {0,0,0};

  gt_assert(querymatchoutoptions != NULL);
  querymatchoutoptions->eoplist.nextfreeuint8_t = 0;
  if (querymatchoutoptions->totallength == GT_UWORD_MAX)
  {
    querymatchoutoptions->totallength = gt_encseq_total_length(encseq);
  }
  ustart = querymatchoutoptions->seedpos1 + querymatchoutoptions->seedlen;
  vstart = querymatchoutoptions->seedpos2 + querymatchoutoptions->seedlen;
  gt_assert(dbstart + dblen >= ustart);
  ulen = dbstart + dblen - ustart;
  gt_assert(querystartabsolute + querylen >= vstart);
  vlen = querystartabsolute + querylen - vstart;
  if (ulen > 0 && vlen > 0)
  {
    rightdistance = align_front_prune_edist(true,
                                            &right_best_polished_point,
                                            querymatchoutoptions->front_trace,
                                            encseq,
                                            querymatchoutoptions->ggemi,
                                            ustart,
                                            ulen,
                                            vstart,
                                            vlen);
    if (rightdistance == ulen + vlen + 1)
    {
      alignment_succeeded = false;
    } else
    {
      front_trace2eoplist(&querymatchoutoptions->eoplist,
                          querymatchoutoptions->front_trace,
                          &right_best_polished_point,
                          ulen,vlen);
    }
    front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
  }
  if (alignment_succeeded)
  {
    GtUword leftdistance;

    front_trace_multireplacement(&querymatchoutoptions->eoplist,
                                 querymatchoutoptions->seedlen);
    if (querymatchoutoptions->seedpos1 > dbstart &&
        querymatchoutoptions->seedpos2 > querystartabsolute)
    {
      ustart = GT_REVERSEPOS(querymatchoutoptions->totallength,
                             querymatchoutoptions->seedpos1 - 1);
      ulen = querymatchoutoptions->seedpos1 - dbstart;
      vstart = GT_REVERSEPOS(querymatchoutoptions->totallength,
                             querymatchoutoptions->seedpos2-1);
      vlen = querymatchoutoptions->seedpos2 - querystartabsolute;
      leftdistance = align_front_prune_edist(false,
                                             &left_best_polished_point,
                                             querymatchoutoptions->front_trace,
                                             encseq,
                                             querymatchoutoptions->ggemi,
                                             ustart,
                                             ulen,
                                             vstart,
                                             vlen);
      if (leftdistance == ulen + vlen + 1)
      {
        alignment_succeeded = false;
      } else
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
      front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
    }
  }
  if (alignment_succeeded)
  {
    GtUword leftcolumn, rightcolumn;

    gt_assert(dbstart <= querymatchoutoptions->seedpos1 -
                         left_best_polished_point.row);
    coords->uoffset = querymatchoutoptions->seedpos1 -
                      left_best_polished_point.row - dbstart;
    coords->ulen = querymatchoutoptions->seedlen +
                   left_best_polished_point.row +
                   right_best_polished_point.row;
    leftcolumn = left_best_polished_point.alignedlen -
                 left_best_polished_point.row;
    rightcolumn = right_best_polished_point.alignedlen -
                  right_best_polished_point.row;
    gt_assert(querymatchoutoptions->seedpos2 >= leftcolumn &&
              querystartabsolute <= querymatchoutoptions->seedpos2 -
                                    leftcolumn);
    coords->voffset = querymatchoutoptions->seedpos2 -
                      leftcolumn - querystartabsolute;
    coords->vlen = querymatchoutoptions->seedlen +
                   leftcolumn + rightcolumn;
  } else
  {
    coords->uoffset = coords->voffset = 0;
    coords->ulen = dblen;
    coords->vlen = querylen;
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

void gt_querymatch_alignment_prepare(GtQuerymatchoutoptions
                                     *querymatchoutoptions,
                                     const GtEncseq *encseq,
                                     GtUword dbstart,
                                     GtUword dblen,
                                     GtUword querystartabsolute,
                                     GtUword querylen,
                                     GtUword edist,
                                     GT_UNUSED bool greedyextension)
{
  gt_assert(querymatchoutoptions != NULL);
  if (dblen > querymatchoutoptions->useqbuffer_size)
  {
    querymatchoutoptions->useqbuffer
      = gt_realloc(querymatchoutoptions->useqbuffer,
                   sizeof *querymatchoutoptions->useqbuffer * dblen);
    querymatchoutoptions->useqbuffer_size = dblen;
  }
  if (querylen > querymatchoutoptions->vseqbuffer_size)
  {
    querymatchoutoptions->vseqbuffer
      = gt_realloc(querymatchoutoptions->vseqbuffer,
                   sizeof *querymatchoutoptions->vseqbuffer * querylen);
    querymatchoutoptions->vseqbuffer_size = querylen;
  }
  if (querymatchoutoptions->esr_for_align_show == NULL)
  {
    querymatchoutoptions->esr_for_align_show
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (querymatchoutoptions->characters == NULL)
  {
    querymatchoutoptions->characters
      = gt_encseq_alphabetcharacters(encseq);
    querymatchoutoptions->wildcardshow
      = gt_alphabet_wildcard_show(gt_encseq_alphabet(encseq));
  }
  gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->esr_for_align_show,
                            encseq,
                            querymatchoutoptions->useqbuffer,
                            dbstart,
                            dbstart + dblen - 1);
  gt_encseq_extract_encoded_with_reader(
                            querymatchoutoptions->esr_for_align_show,
                            encseq,
                            querymatchoutoptions->vseqbuffer,
                            querystartabsolute,
                            querystartabsolute + querylen -1);
  if (edist > 0)
  {
    GtSeqpaircoordinates coords;

    if (seededmatch2eoplist(&coords,
                            dbstart,
                            dblen,
                            querystartabsolute,
                            querylen,
                            querymatchoutoptions,
                            encseq))
    {
      gt_alignment_set_seqs(querymatchoutoptions->alignment,
                            querymatchoutoptions->useqbuffer +
                                coords.uoffset,
                            coords.ulen,
                            querymatchoutoptions->vseqbuffer +
                                coords.voffset,
                            coords.vlen);
      converteoplist2alignment(querymatchoutoptions->alignment,
                               &querymatchoutoptions->eoplist);
    } else
    {
      GtUword linedist;

      gt_assert(!greedyextension);
      gt_alignment_set_seqs(querymatchoutoptions->alignment,
                            querymatchoutoptions->useqbuffer,
                            dblen,
                            querymatchoutoptions->vseqbuffer,
                            querylen);
      linedist = gt_computelinearspace(querymatchoutoptions->alignment,
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
#ifdef SKDEBUG
    check_correct_edist(querymatchoutoptions->useqbuffer,
                        dblen,
                        querymatchoutoptions->vseqbuffer,
                        querylen,
                        edist);
#endif
  }
}

void gt_querymatchoutoptions_set_seed(
                            GtQuerymatchoutoptions *querymatchoutoptions,
                            GtUword pos1,GtUword pos2,GtUword len)
{
  gt_assert(querymatchoutoptions != NULL);
  querymatchoutoptions->seedpos1 = pos1;
  querymatchoutoptions->seedpos2 = pos2;
  querymatchoutoptions->seedlen = len;
}

void gt_querymatchoutoptions_alignment_show(const GtQuerymatchoutoptions
                                              *querymatchoutoptions)
{
  gt_assert(querymatchoutoptions != NULL);
  gt_alignment_show_generic(querymatchoutoptions->alignment_show_buffer,
                            querymatchoutoptions->alignment,
                            stdout,
                            (unsigned int)
                            querymatchoutoptions->alignmentwidth,
                            querymatchoutoptions->characters,
                            querymatchoutoptions->wildcardshow);
  gt_alignment_reset(querymatchoutoptions->alignment);
}

void gt_querymatchoutoptions_exact_alignment_show(const GtQuerymatchoutoptions
                                                     *querymatchoutoptions,
                                                  GtUword len)
{
  gt_assert(querymatchoutoptions != NULL);
  gt_alignment_exact_show(querymatchoutoptions->alignment_show_buffer,
                          querymatchoutoptions->useqbuffer,
                          len,
                          stdout,
                          querymatchoutoptions->alignmentwidth,
                          querymatchoutoptions->characters);
}
