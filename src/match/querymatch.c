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
#include "core/readmode.h"
#include "core/format64.h"
#include "core/minmax.h"
#include "extended/alignment.h"
#include "extended/linearalign.h"
#undef SKDEBUG
#ifdef SKDEBUG
#include "match/greedyedist.h"
#endif
#include "querymatch.h"
#include "seed-extend.h"
#include "ft-front-generation.h"

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
          useqbuffer_size, vseqbuffer_size;
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
  querymatchoutoptions->totallength = GT_UWORD_MAX;
  querymatchoutoptions->esr_for_align_show = NULL;
  querymatchoutoptions->characters = NULL;
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

static void gt_querymatch_alignment_prepare(GtQuerymatchoutoptions
                                             *querymatchoutoptions,
                                            const GtEncseq *encseq,
                                            GtUword dbstart,
                                            GtUword dblen,
                                            GtUword querystartabsolute,
                                            GtUword querylen,
                                            GtUword edist,
                                            GT_UNUSED bool greedyextension)
{
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

void gt_querymatch_set_seed(GtQuerymatchoutoptions *querymatchoutoptions,
                            GtUword pos1,GtUword pos2,GtUword len)
{
  querymatchoutoptions->seedpos1 = pos1;
  querymatchoutoptions->seedpos2 = pos2;
  querymatchoutoptions->seedlen = len;
}

struct GtQuerymatch
{
   GtUword
      dblen, /* length of match in dbsequence */
      querylen, /* same as dblen for exact matches */
      dbstart, /* absolute start position of match in database seq */
      querystart, /* start of match in query, relative to start of query */
      edist, /* 0 for exact match */
      dbseqnum, querystart_fwdstrand, dbstart_relative;
   GtWord score; /* 0 for exact match */
   uint64_t queryseqnum; /* ordinal number of match in query */
   double similarity;
   GtReadmode readmode; /* readmode by which reference sequence was accessed */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        query_as_reversecopy, /* matched the reverse copy of the query */
        greedyextension;
};

GtQuerymatch *gt_querymatch_new(void)
{
  return gt_malloc(sizeof (GtQuerymatch));
}

GtUword gt_querymatch_dbseqnum(const GtEncseq *encseq,
                               const GtQuerymatch *querymatch)
{
  return gt_encseq_seqnum(encseq,querymatch->dbstart);
}

void gt_querymatch_fill(GtQuerymatch *querymatch,
                        const GtEncseq *encseq,
                        GtUword dblen,
                        GtUword dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength)
{
  GtUword dbseqstartpos;

  querymatch->dblen = dblen;
  querymatch->dbstart = dbstart;
  querymatch->readmode = readmode;
  querymatch->query_as_reversecopy = query_as_reversecopy;
  querymatch->score = score;
  querymatch->edist = edist;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->greedyextension = false;
  gt_assert(encseq != NULL);
  if (gt_encseq_has_multiseq_support(encseq))
  {
    querymatch->dbseqnum = gt_querymatch_dbseqnum(encseq,querymatch);
    dbseqstartpos = gt_encseq_seqstartpos(encseq, querymatch->dbseqnum);
  } else
  {
    querymatch->dbseqnum = 0;
    dbseqstartpos = 0;
  }
  gt_assert((int) querymatch->readmode < 4);
  if (querymatch->readmode == GT_READMODE_REVERSE ||
      querymatch->readmode == GT_READMODE_REVCOMPL)
  {
    gt_assert(querymatch->querystart + querymatch->querylen <=
              query_totallength);
    querymatch->querystart_fwdstrand
      = query_totallength - querymatch->querystart - querymatch->querylen;
  } else
  {
    querymatch->querystart_fwdstrand = querymatch->querystart;
  }
  gt_assert(querymatch->dbstart >= dbseqstartpos);
  querymatch->dbstart_relative = querymatch->dbstart - dbseqstartpos;
  if (querymatch->edist == 0)
  {
    querymatch->similarity = 100.0;
  } else
  {
    querymatch->similarity
      = 100.0 - gt_querymatch_error_rate(querymatch->edist,
                                          querymatch->dblen +
                                          querymatch->querylen);
  }
}

void gt_querymatch_delete(GtQuerymatch *querymatch)
{
  if (querymatch != NULL)
  {
    gt_free(querymatch);
  }
}

#ifdef VERIFY
static void verifyexactmatch(const GtEncseq *encseq,
                             GtUword len,
                             GtUword pos1,
                             uint64_t seqnum2,
                             GtUword pos2,
                             GtReadmode readmode)
{
  if (readmode == GT_READMODE_REVERSE)
  {
    GtUword offset, seqstartpos, totallength = gt_encseq_total_length(encseq);
    GtUchar cc1, cc2;

    seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
    pos2 += seqstartpos;
    for (offset = 0; offset < len; offset++)
    {
      gt_assert(pos1 + len - 1 < totallength);
      gt_assert(pos2 + len - 1 < totallength);
      cc1 = gt_encseq_get_encoded_char(encseq,pos1+offset,GT_READMODE_FORWARD);
      cc2 = gt_encseq_get_encoded_char(encseq,pos2+len-1-offset,
                                       GT_READMODE_FORWARD);
      gt_assert(cc1 == cc2 && ISNOTSPECIAL(cc1));
    }
    if (pos1 + len < totallength)
    {
      cc1 = gt_encseq_get_encoded_char(encseq,pos1+len,GT_READMODE_FORWARD);
    } else
    {
      cc1 = SEPARATOR;
    }
    if (pos2 > 0)
    {
      cc2 = gt_encseq_get_encoded_char(encseq,pos2-1,GT_READMODE_FORWARD);
    } else
    {
      cc2 = SEPARATOR;
    }
    gt_assert(cc1 != cc2 || ISSPECIAL(cc1));
  }
}
#endif

int gt_querymatch_output(void *info,
                         const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         GT_UNUSED GtUword query_totallength,
                         GT_UNUSED GtError *err)
{
  if (!querymatch->selfmatch ||
      (uint64_t) querymatch->dbseqnum != querymatch->queryseqnum ||
      querymatch->dbstart_relative <= querymatch->querystart_fwdstrand)
  {
    GtQuerymatchoutoptions *querymatchoutoptions
      = (GtQuerymatchoutoptions *) info;
    const char *outflag = "FRCP";

#ifdef VERIFY
    verifyexactmatch(encseq,
                     querymatch->len,
                     querymatch->dbstart,
                     querymatch->queryseqnum,
                     querystart_fwdstrand,
                     querymatch->readmode);
#endif
    if (querymatchoutoptions != NULL &&
        querymatchoutoptions->alignmentwidth > 0)
    {
      if (querymatch->selfmatch)
      {
        GtUword querystartabsolute
          = gt_encseq_seqstartpos(encseq,querymatch->queryseqnum) +
            querymatch->querystart_fwdstrand;
        gt_querymatch_alignment_prepare(querymatchoutoptions,
                                        encseq,
                                        querymatch->dbstart,
                                        querymatch->dblen,
                                        querystartabsolute,
                                        querymatch->querylen,
                                        querymatch->edist,
                                        querymatch->greedyextension);
      } else
      {
        gt_assert(false); /* case not implemented yet */
      }
    }
    printf(GT_WU " " GT_WU " " GT_WU " %c " GT_WU " " Formatuint64_t
                   " " GT_WU,
           querymatch->dblen,
           querymatch->dbseqnum,
           querymatch->dbstart_relative,
           outflag[querymatch->readmode],
           querymatch->querylen,
           PRINTuint64_tcast(querymatch->queryseqnum),
           querymatch->querystart_fwdstrand);
    if (querymatch->score > 0)
    {
      printf(" " GT_WD " " GT_WU " %.2f",
             querymatch->score,querymatch->edist,querymatch->similarity);
    }
    printf("\n");
    if (querymatchoutoptions != NULL &&
        querymatchoutoptions->alignmentwidth > 0)
    {
      if (querymatch->selfmatch)
      {
        if (querymatch->edist > 0)
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
          gt_assert(querymatch->dblen == querymatch->querylen);
          gt_alignment_exact_show(querymatchoutoptions->alignment_show_buffer,
                                  querymatchoutoptions->useqbuffer,
                                  querymatch->dblen,
                                  stdout,
                                  querymatchoutoptions->alignmentwidth,
                                  querymatchoutoptions->characters);
        }
      } else
      {
        gt_assert(false); /* case not implemented yet */
      }
    }
  }
  return 0;
}

int gt_querymatch_fill_and_output(
                        GtUword dblen,
                        GtUword dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtQuerymatchoutoptions *querymatchoutoptions,
                        const GtEncseq *encseq,
                        const GtUchar *query,
                        GtUword query_totallength,
                        bool greedyextension,
                        GtError *err)
{
  GtQuerymatch querymatch;

  gt_querymatch_fill(&querymatch,
                     encseq,
                     dblen,
                     dbstart,
                     readmode,
                     query_as_reversecopy,
                     score,
                     edist,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength);
  querymatch.greedyextension = greedyextension;
  return gt_querymatch_output(querymatchoutoptions,
                       encseq,
                       &querymatch,
                       query,
                       query_totallength,
                       err);
}

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch)
{
  return querymatch->querylen;
}

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch)
{
  return querymatch->dbstart;
}

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch)
{
  return querymatch->querystart;
}

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch)
{
  return querymatch->queryseqnum;
}

bool gt_querymatch_queryreverse(const GtQuerymatch *querymatch)
{
  return querymatch->query_as_reversecopy;
}

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
}
