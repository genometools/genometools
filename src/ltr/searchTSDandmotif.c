/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#include "core/arraydef.h"
#include "core/error.h"
#include "core/log.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/encseq.h"
#include "match/spacedef.h"
#include "match/esa-mmsearch.h"

#include "ltrharvest-opt.h"
#include "repeattypes.h"
#include "repeats.h"
#include "searchTSDandmotif.h"

/*
 The following function searches for TSDs and/or a specified palindromic
 motif at the 5'-border of left LTR and 3'-border of right LTR. Thereby,
 all maximal repeats from the vicinity are processed one after another
 to find the TSD with the minimum deviation with regard to the boundary
 position from the x-drop alignment. If also a motif is searched,
 a simple motif check at the boundaries of the TSDs is performed.
 */

static void searchforbestTSDandormotifatborders(SubRepeatInfo *info,
                                                LTRharvestoptions *lo,
                                                const GtEncseq *encseq,
                                                LTRboundaries *boundaries,
                                                unsigned int
                                                *motifmismatchesleftLTR,
                                                unsigned int
                                                *motifmismatchesrightLTR)
{
  unsigned long i;
  unsigned long motifpos1,
         motifpos2,
         back, forward,
         oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary1 = 0,
         difffromoldboundary2 = 0;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR,
               hitcounter = 0;

  if (info->repeats.nextfreeRepeat > 0 )
  {
    boundaries->tsd = true;
  }
  boundaries->motif_near_tsd = false;

  for (i = 0; i < info->repeats.nextfreeRepeat; i++)
  {
    /* motifpos1 is the first position after the left repeat */
    motifpos1 = info->repeats.spaceRepeat[i].pos1 +
                           info->repeats.spaceRepeat[i].len;
    /* motifpos2 is two positions before the right repeat */
    motifpos2 = info->repeats.spaceRepeat[i].pos1
                  + info->repeats.spaceRepeat[i].offset - 2;

    for (back = 0;
         back < info->repeats.spaceRepeat[i].len - info->lmin + 1;
         back++)
    {
      for (forward = 0;
           forward < info->repeats.spaceRepeat[i].len -
                     info->lmin + 1 - back;
           forward++)
      {
        tmp_motifmismatchesleftLTR = tmp_motifmismatchesrightLTR = 0;
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                            motifpos1 - back, GT_READMODE_FORWARD)
            != lo->motif.firstleft)
        {
          tmp_motifmismatchesleftLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos1 + 1 - back,
                                              GT_READMODE_FORWARD)
            != lo->motif.secondleft)
        {
          tmp_motifmismatchesleftLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos2 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif.firstright)
        {
          tmp_motifmismatchesrightLTR++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ encseq,
                                              motifpos2 + 1 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif.secondright)
        {
          tmp_motifmismatchesrightLTR++;
        }

        if (tmp_motifmismatchesleftLTR <= lo->motif.allowedmismatches
            &&
            tmp_motifmismatchesrightLTR <= lo->motif.allowedmismatches)
        {
          unsigned long tsd_len;
          tsd_len = info->repeats.spaceRepeat[i].len - back - forward;

          /* TSD length not too big */
          if (tsd_len <= (unsigned long)info->lmax)
          {
            if ( !boundaries->motif_near_tsd )
            {
              unsigned long max, min;

              /* save number of mismatches */
              *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
              *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

              /* adjust boundaries */
              boundaries->motif_near_tsd = true;
              boundaries->leftLTR_5  = motifpos1 - back;
              boundaries->rightLTR_3 = motifpos2 + 1 + forward;

              /* store TSD length */
              boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;

              max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
              min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
              difffromoldboundary1 = max - min;

              max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
              min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
              difffromoldboundary2 = max - min;

              hitcounter++;
            }
            else
            {
              unsigned long max, min, difffromnewboundary1,
                   difffromnewboundary2;

              /* test if hit is nearer to old boundaries than previous hit */
              max = MAX(oldleftLTR_5, (motifpos1 - back));
              min = MIN(oldleftLTR_5, (motifpos1 - back));
              difffromnewboundary1 = max - min;
              max = MAX(oldrightLTR_3, (motifpos2 + 1 + forward));
              min = MIN(oldrightLTR_3, (motifpos2 + 1 + forward));
              difffromnewboundary2 = max - min;

              if (difffromnewboundary1 + difffromnewboundary2 <
                  difffromoldboundary1 + difffromoldboundary2)
              {
                /* save number of mismatches */
                *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
                *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

                /* adjust boundaries */
                boundaries->leftLTR_5  = motifpos1 - back;
                boundaries->rightLTR_3 = motifpos2 + 1 + forward;

                /* store TSD length */
                boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;

                difffromoldboundary1 = difffromnewboundary1;
                difffromoldboundary2 = difffromnewboundary2;
                hitcounter++;
              }
            }
          }
        }
      }
    }
  }
}

/*
 The following function searches only for a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static void searchformotifonlyborders(LTRharvestoptions *lo,
    LTRboundaries *boundaries,
    const GtEncseq *encseq,
    unsigned long startleftLTR,
    unsigned long endleftLTR,
    unsigned long startrightLTR,
    unsigned long endrightLTR,
    unsigned int *motifmismatchesleftLTR,
    unsigned int *motifmismatchesrightLTR
    )
{
  bool motif1 = false,
       motif2 = false;
  unsigned int tmp_motifmismatchesleftLTR,
         tmp_motifmismatchesrightLTR,
         motifmismatches_frombestmatch = 0;
  unsigned long idx,
         oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary = 0;

  /**** search for left motif around leftLTR_5 ****/

  for (idx = startleftLTR; idx < endleftLTR; idx++)
  {
    tmp_motifmismatchesleftLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
        != lo->motif.firstleft)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx+1,
                                          GT_READMODE_FORWARD) !=
        lo->motif.secondleft)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR)
                                <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if ( !motif1 )
       {
         unsigned long max, min;

         motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
         boundaries->leftLTR_5 = idx;
         motif1 = true;
         max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
         min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
         difffromoldboundary = max - min;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_5, idx);
         minval = MIN(oldleftLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
           boundaries->leftLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  for (idx = startrightLTR + 1; idx <= endrightLTR; idx++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD) !=
                       lo->motif.secondright)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx-1,
                                          GT_READMODE_FORWARD) !=
                       lo->motif.firstright)
    {
      tmp_motifmismatchesrightLTR++;
    }

    if (tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if ( !motif2 )
       {
         unsigned long max, min;

         motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
         boundaries->rightLTR_3 = idx;
         motif2 = true;
         max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
         min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
         difffromoldboundary = max - min;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_3, idx);
         minval = MIN(oldrightLTR_3, idx);
         difffromnewboundary = maxval - minval;
         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
           boundaries->rightLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  if (idx > endrightLTR && (!motif2))
  {
    gt_log_log("no right motif found.\n");
  }
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;

  if (motif1 && motif2)
  {
    boundaries->motif_near_tsd = true;
  }
  else
  {
    boundaries->motif_near_tsd = false;
  }
}

/*
 The following function searches for a specified palindromic motif at the
 3'-border of left LTR and the 5'-border of right LTR.
 */

static void searchformotifonlyinside(LTRharvestoptions *lo,
                                     LTRboundaries *boundaries,
                                     const GtEncseq *encseq,
                                     unsigned int *motifmismatchesleftLTR,
                                     unsigned int *motifmismatchesrightLTR)
{
  bool motif1 = false,
       motif2 = false;
  unsigned long startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         oldleftLTR_3  = boundaries->leftLTR_3,
         oldrightLTR_5 = boundaries->rightLTR_5,
         difffromoldboundary = 0,
         idx;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR,
               motifmismatches_frombestmatch = 0;

  /** vicinity of 3'-border of left LTR **/
  /* do not align over 5'border of left LTR,
     in case of need decrease alignment length */
  if ( (startleftLTR = boundaries->leftLTR_3 -
         lo->vicinityforcorrectboundaries) <
      boundaries->leftLTR_5 + 2)
  {
    startleftLTR = boundaries->leftLTR_5 + 2;
  }
  /* do not align over 5'-border of right LTR */
  if ( (endleftLTR = boundaries->leftLTR_3 +
       lo->vicinityforcorrectboundaries) >
      boundaries->rightLTR_5 - 1)
  {
    endleftLTR = boundaries->rightLTR_5 - 1;
  }
  /** vicinity of 5'-border of right LTR **/
  /* do not align over 3'-border of left LTR */
  if ( (startrightLTR = boundaries->rightLTR_5 -
         lo->vicinityforcorrectboundaries)
       < boundaries->leftLTR_3 + 1)
  {
    startrightLTR = boundaries->leftLTR_3 + 1;
  }
  /* do not align over 3'border of right LTR */
  if ( (endrightLTR = boundaries->rightLTR_5 +
       lo->vicinityforcorrectboundaries) >
      boundaries->rightLTR_3 - 2)
  {
    endrightLTR = boundaries->rightLTR_3 - 2;
  }

  /**** search for right motif around leftLTR_3 ****/

  for (idx = startleftLTR + 1; idx <= endleftLTR; idx++)
  {
    tmp_motifmismatchesleftLTR = (unsigned int)0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
                       != lo->motif.secondright)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq,
                                          idx-1,
                                          GT_READMODE_FORWARD) !=
                       lo->motif.firstright)
    {
      tmp_motifmismatchesleftLTR++;
    }
    if (tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR)
                                <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if ( !motif1 )
       {
         unsigned long maxval, minval;

         motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
         boundaries->leftLTR_3 = idx;
         motif1 = true;
         maxval = MAX(oldleftLTR_3, boundaries->leftLTR_3);
         minval = MIN(oldleftLTR_3, boundaries->leftLTR_3);
         difffromoldboundary = maxval - minval;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_3, idx);
         minval = MIN(oldleftLTR_3, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
           boundaries->leftLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  /**** search for left motif around rightLTR_5 ****/

  for (idx = startrightLTR ; idx < endrightLTR; idx++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx,
                                          GT_READMODE_FORWARD)
                       != lo->motif.firstleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (gt_encseq_get_encoded_char(/* XXX */ encseq, idx+1,
                                          GT_READMODE_FORWARD)
                       != lo->motif.secondleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if (tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if ( !motif2 )
       {
         unsigned long maxval, minval;

         motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
         boundaries->rightLTR_5 = idx;
         motif2 = true;
         maxval = MAX(oldrightLTR_5, boundaries->rightLTR_5);
         minval = MIN(oldrightLTR_5, boundaries->rightLTR_5);
         difffromoldboundary = maxval - minval;
       }
       /* next hit */
       else
       {
         unsigned long maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_5, idx);
         minval = MIN(oldrightLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if ( difffromnewboundary < difffromoldboundary )
         {
           motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
           boundaries->rightLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;
  if (motif1 && motif2)
  {
    boundaries->motif_far_tsd = true;
  }
  else
  {
    boundaries->motif_far_tsd = false;
  }
}

/*
 The following function searches for TSDs and/or a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static int searchforTSDandorMotifoutside(
  LTRharvestoptions *lo,
  LTRboundaries *boundaries,
  const GtEncseq *encseq,
  unsigned int *motifmismatchesleftLTR,
  unsigned int *motifmismatchesrightLTR,
  GtError *err)
{
  unsigned long startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         leftlen,
         rightlen,
         sequenceendpos,
         seqstartpos,
         seqlength;
  unsigned long contignumber = boundaries->contignumber;
  SubRepeatInfo subrepeatinfo;
  bool haserr = false;

  gt_error_check(err);

  /* check border cases */

  /* vicinity of 5'-border of left LTR */
  seqstartpos = gt_encseq_seqstartpos(encseq, contignumber);
  seqlength = gt_encseq_seqlength(encseq, contignumber);
  if (contignumber == 0)
  {
    /* do not align over left sequence boundary,
       in case of need decrease alignment length */
    if ( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries)
    {
      startleftLTR = seqstartpos;
    } else
    {
      startleftLTR = boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries;
    }
  }
  else
  {
    /* do not align over left separator symbol
       in case of need decrease alignment length */
    if ( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries )
    {
      startleftLTR = seqstartpos;
    }
    else
    {
      if ( ((startleftLTR =
              boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries) <
                seqstartpos)
            &&
          (boundaries->leftLTR_5 >= seqstartpos)
        )
      {
        startleftLTR = seqstartpos;
      }
    }
  }
  /* do not align over 3'-border of left LTR */
  if ( (endleftLTR = boundaries->leftLTR_5 + lo->vicinityforcorrectboundaries)
       > boundaries->leftLTR_3 - 2 /* -2 because of possible motif */
    )
  {
    endleftLTR = boundaries->leftLTR_3 - 2;
  }
  leftlen = endleftLTR - startleftLTR + 1;

  /* vicinity of 3'-border of right LTR
     do not align over 5'border of right LTR */
  if ( (startrightLTR =
         boundaries->rightLTR_3 - lo->vicinityforcorrectboundaries) <
      boundaries->rightLTR_5 + 2  /* +2 because of possible motif */
    )
  {
    startrightLTR = boundaries->rightLTR_5 + 2;
  }
  sequenceendpos = seqstartpos + seqlength - 1;
  /* do not align into next sequence in case of need decrease alignment
     length */
  endrightLTR = boundaries->rightLTR_3 + lo->vicinityforcorrectboundaries;
  if (endrightLTR > sequenceendpos &&
      boundaries->rightLTR_3 <= sequenceendpos)
  {
    endrightLTR = sequenceendpos;
  }
  rightlen = endrightLTR - startrightLTR + 1;

  /* now, search for correct boundaries */

  /* search for TSDs and/or motif */
  if (lo->minlengthTSD > 1U)
  {
    GtUchar *dbseq, *query;
    ALLOCASSIGNSPACE(dbseq,NULL,GtUchar,leftlen);
    ALLOCASSIGNSPACE(query,NULL,GtUchar,rightlen);

    gt_encseq_extract_substring(encseq,dbseq,startleftLTR,endleftLTR);
    gt_encseq_extract_substring(encseq,query,startrightLTR,
                                         endrightLTR);
    GT_INITARRAY(&subrepeatinfo.repeats, Repeat);
    subrepeatinfo.lmin = (unsigned long) lo->minlengthTSD;
    subrepeatinfo.lmax = (unsigned long) lo->maxlengthTSD;
    gt_assert(startleftLTR < startrightLTR);
    subrepeatinfo.offset1 = startleftLTR;
    subrepeatinfo.offset2 = startrightLTR;

    if (gt_sarrquerysubstringmatch(dbseq,
                                leftlen,
                                query,
                                (unsigned long) rightlen,
                                lo->minlengthTSD,
                                gt_encseq_alphabet(encseq),
                                gt_subsimpleexactselfmatchstore,
                                &subrepeatinfo,
                                NULL,
                                err) != 0)
    {
       haserr = true;
    }

    FREESPACE(dbseq);
    FREESPACE(query);

    if (!haserr)
    {
      searchforbestTSDandormotifatborders(&subrepeatinfo,
                                          lo,
                                          encseq,
                                          boundaries,
                                          motifmismatchesleftLTR,
                                          motifmismatchesrightLTR);
    }
    GT_FREEARRAY (&subrepeatinfo.repeats, Repeat);
  } else /* no search for TSDs, search for motif only */
  {
    searchformotifonlyborders(lo,
                              boundaries,
                              encseq,
                              startleftLTR,
                              endleftLTR,
                              startrightLTR,
                              endrightLTR,
                              motifmismatchesleftLTR,
                              motifmismatchesrightLTR);
  }
  return haserr ? -1 : 0;
}

/*
 The following function searches for TSD and/or a specified palindromic motif
 at the borders of left LTR and the right LTR, respectively.
 */
int gt_findcorrectboundaries(LTRharvestoptions *lo,
                          LTRboundaries *boundaries,
                          const GtEncseq *encseq,
                          GtError *err)
{
  unsigned int motifmismatchesleftLTR = 0,
               motifmismatchesrightLTR = 0;

  gt_error_check(err);
  gt_log_log("searching for correct boundaries in vicinity...\n");
  /* first: 5'-border of left LTR and 3'-border of right LTR */

  if (searchforTSDandorMotifoutside(lo,
                                    boundaries,
                                    encseq,
                                    &motifmismatchesleftLTR,
                                    &motifmismatchesrightLTR,
                                    err) != 0 )
  {
    return -1;
  }

  /* second: 3'-border of left LTR and 5'-border of right LTR */
  if ( lo->motif.allowedmismatches < (unsigned int)4 )
  {
    gt_log_log("second: searching for motif only around 3'border of left LTR "
               "and 5'-border of right LTR...\n");
    searchformotifonlyinside(lo,boundaries,encseq,
                             &motifmismatchesleftLTR,
                             &motifmismatchesrightLTR);
  }
  return 0;
}
