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

#include "core/arraydef.h"
#include "core/error.h"
#include "core/minmax.h"
#include "core/mathsupport.h"
#include "core/encodedsequence.h"
#include "match/spacedef.h"

#include "searchforLTRs.h"
#include "repeattypes.h"
#include "repeats.h"
#include "myxdrop.h"
#include "searchTSDandmotif.h"

#include "match/greedyedist.pr"

/*
 The following function checks, if the remaining candidate pairs still
 fulfill the length and distance constraints of LTRs.
 */

static bool checklengthanddistanceconstraints(LTRboundaries *boundaries,
                                             RepeatInfo *repeatinfo)
{
  unsigned long ulen, vlen, dist_between_LTRs;

  ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1;
  vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
  dist_between_LTRs = boundaries->rightLTR_5 - boundaries->leftLTR_5;
  if (ulen > (unsigned long)repeatinfo->lmax
        || vlen > (unsigned long)repeatinfo->lmax
        || ulen < (unsigned long)repeatinfo->lmin
        || vlen < (unsigned long)repeatinfo->lmin
        || dist_between_LTRs > (unsigned long)repeatinfo->dmax
        || dist_between_LTRs < (unsigned long)repeatinfo->dmin
        || boundaries->leftLTR_3 >= boundaries->rightLTR_5)
  {
    boundaries->lengthdistconstraint = false;
    boundaries->similarity = 0.0;
    return false;
  } else
  {
    boundaries->lengthdistconstraint = true;
  }
  return true;
}

static void adjustboundariesfromXdropextension(Myxdropbest xdropbest_left,
                                               Myxdropbest xdropbest_right,
                                               unsigned long seed1_startpos,
                                               unsigned long seed2_startpos,
                                               unsigned long seed1_endpos,
                                               unsigned long seed2_endpos,
                                               LTRboundaries *boundaries)
{
  /* left alignment */
  boundaries->leftLTR_5  = seed1_startpos - xdropbest_left.ivalue;
  boundaries->rightLTR_5 = seed2_startpos - xdropbest_left.jvalue;

  /* right alignment */
  boundaries->leftLTR_3  = seed1_endpos + xdropbest_right.ivalue;
  boundaries->rightLTR_3 = seed2_endpos + xdropbest_right.jvalue;

  /*
  printf("Adjusted boundaries from Xdrop-alignment-extension:\n");
  printf("boundaries->leftLTR_5  = " FormatSeqpos "\n",
            PRINTSeqposcast( boundaries->leftLTR_5 - offset));
  printf("boundaries->leftLTR_3  = " FormatSeqpos "\n",
            PRINTSeqposcast( boundaries->leftLTR_3 - offset));
  printf("len leftLTR = " FormatSeqpos "\n",
            PRINTSeqposcast( (boundaries->leftLTR_3 - boundaries->leftLTR_5
                              + 1)));
  printf("boundaries->rightLTR_5 = " FormatSeqpos "\n",
            PRINTSeqposcast( boundaries->rightLTR_5 - offset));
  printf("boundaries->rightLTR_3 = " FormatSeqpos "\n",
            PRINTSeqposcast( boundaries->rightLTR_3 - offset));
  printf("len rightLTR = " FormatSeqpos "\n",
            PRINTSeqposcast( (boundaries->rightLTR_3 - boundaries->rightLTR_5
                              + 1)));
 */
}

/*
 The following function applies the filter algorithms one after another
 to all candidate pairs.
*/
int gt_searchforLTRs(LTRharvestoptions *lo,
                  GtArrayLTRboundaries *arrayLTRboundaries,
                  const GtEncodedsequence *encseq,
                  GtError *err)
{
  unsigned long repeatcounter;
  GtArrayMyfrontvalue fronts;
  Myxdropbest xdropbest_left;
  Myxdropbest xdropbest_right;
  unsigned long alilen = 0,
         totallength,
         ulen,
         vlen,
         maxulen = 0,
         maxvlen = 0;
  GtUchar *useq = NULL,
        *vseq = NULL;
  unsigned long edist;
  Repeat *repeatptr;
  LTRboundaries *boundaries;
  bool haserr = false;

  gt_error_check(err);

  /*
  printf("xdropbelowscore = %d\n", lo->xdropbelowscore);
  printf("scores:\n");
  printf("mat  = %d\n", lo->arbitscores.mat);
  printf("mis  = %d\n", lo->arbitscores.mis);
  printf("ins  = %d\n", lo->arbitscores.ins);
  printf("del  = %d\n", lo->arbitscores.del);
  */

  for (repeatcounter = 0; repeatcounter < lo->repeatinfo.repeats.nextfreeRepeat;
       repeatcounter++)
  {
    /* printf("\n\nAlignments of repeat nr. = %u :\n", repeatcounter+1); */

    repeatptr = &(lo->repeatinfo.repeats.spaceRepeat[repeatcounter]);
    alilen = ((unsigned long)lo->repeatinfo.lmax) - repeatptr->len;

    /**** left (reverse) xdrop alignment ****/
    GT_INITARRAY (&fronts, Myfrontvalue);
    if (alilen <= repeatptr->pos1)
    {
      gt_evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) alilen,
                               (int) alilen,
                               (Xdropscore)lo->xdropbelowscore);
    }
    else /* do not align over left sequence boundary */
    {
      gt_evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) repeatptr->pos1,
                               (int) (repeatptr->pos1 + repeatptr->offset),
                               (Xdropscore)lo->xdropbelowscore);
    }
    GT_FREEARRAY (&fronts, Myfrontvalue);

    /**** right xdrop alignment ****/
    GT_INITARRAY (&fronts, Myfrontvalue);
    totallength = gt_encodedsequence_total_length(encseq);
    if (alilen <= totallength - (repeatptr->pos1 + repeatptr->offset +
                                repeatptr->len) )
    {
      gt_evalxdroparbitscoresright (&lo->arbitscores,
                                 &xdropbest_right,
                                 &fronts,
                                 encseq,
                                 encseq,
                                 repeatptr->pos1 + repeatptr->len,
                                 repeatptr->pos1 + repeatptr->offset +
                                 repeatptr->len,
                                 (int) alilen,
                                 (int) alilen,
                                 lo->xdropbelowscore);
    }
    else /* do not align over right sequence boundary */
    {
      gt_evalxdroparbitscoresright(&lo->arbitscores,
                                &xdropbest_right,
                                &fronts,
                                encseq,
                                encseq,
                                repeatptr->pos1 + repeatptr->len,
                                repeatptr->pos1 + repeatptr->offset +
                                repeatptr->len,
                                (int) (totallength -
                                (repeatptr->pos1 + repeatptr->len)),
                                (int) (totallength -
                                (repeatptr->pos1 + repeatptr->offset +
                                 repeatptr->len)),
                                lo->xdropbelowscore);
    }
    GT_FREEARRAY (&fronts, Myfrontvalue);

    GT_GETNEXTFREEINARRAY(boundaries,arrayLTRboundaries,LTRboundaries,5);

    boundaries->contignumber = repeatptr->contignumber;
    boundaries->leftLTR_5 = (unsigned long) 0;
    boundaries->leftLTR_3 = (unsigned long) 0;
    boundaries->rightLTR_5 = (unsigned long) 0;
    boundaries->rightLTR_3 = (unsigned long) 0;
    boundaries->lenleftTSD = (unsigned long) 0;
    boundaries->lenrightTSD = (unsigned long) 0;
    boundaries->tsd = false;
    boundaries->motif_near_tsd = false;
    boundaries->motif_far_tsd = false;
    boundaries->skipped = false;
    boundaries->similarity = 0.0;

    /* test
    printf("contig number: %lu\n",
               boundaries->contignumber);
    printf("offset to contig startpos: " FormatSeqpos "\n",
               PRINTSeqposcast(offset));

    printf("Boundaries from vmatmaxoutdynamic:\n");

    printf("boundaries->leftLTR_5 abs. = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1));
    printf("boundaries->rightLTR_5 abs. = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset));

    printf("boundaries->leftLTR_5  = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 - offset));
    printf("boundaries->leftLTR_3  = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 + repeatptr->len - 1
                        - offset));
    printf("boundaries->rightLTR_5 = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset
                        - offset));
    printf("boundaries->rightLTR_3 = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset +
                        repeatptr->len - 1 - offset));
    */

    /* store new boundaries-positions in boundaries */
    adjustboundariesfromXdropextension(
                   xdropbest_left,
                   xdropbest_right,
                   repeatptr->pos1,  /*seed1 startpos*/
                   repeatptr->pos1 + repeatptr->offset, /*seed2 startpos*/
                   repeatptr->pos1 + repeatptr->len - 1, /*seed1 endpos*/
                   repeatptr->pos1 + repeatptr->offset + repeatptr->len - 1,
                                                         /*seed2 endpos*/
                   boundaries);

    /* if search for motif and/or TSD */
    if ( lo->motif.allowedmismatches < 4U || lo->minlengthTSD > 1U)
    {
      if ( gt_findcorrectboundaries(lo, boundaries, encseq, err) != 0 )
      {
        haserr = true;
        break;
      }

      /* if search for TSDs and (not) motif */
      if ( boundaries->tsd &&
          (lo->motif.allowedmismatches >= (unsigned int)4 ||
          (boundaries->motif_near_tsd && boundaries->motif_far_tsd)) )
      {
        /* predicted as full LTR-pair, keep it */
      }
      else
      {
        /* if search for motif only (and not TSD) */
        if ( lo->minlengthTSD <= 1U &&
            boundaries->motif_near_tsd &&
            boundaries->motif_far_tsd )
        {
          /* predicted as full LTR-pair, keep it */
        }
        else
        {
          /* delete this LTR-pair candidate */
          arrayLTRboundaries->nextfreeLTRboundaries--;
          continue;
        }
      }
    }

    /* check length and distance constraints again */
    if (!checklengthanddistanceconstraints(boundaries, &lo->repeatinfo))
    {
      /* delete this LTR-pair candidate */
      arrayLTRboundaries->nextfreeLTRboundaries--;
      continue;
    }

    /* check similarity from candidate pair
       copy LTR sequences for greedyunitedist function */
    ulen = boundaries->leftLTR_3 - boundaries->leftLTR_5 + 1;
    vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
    if (ulen > maxulen)
    {
      maxulen = ulen;
      ALLOCASSIGNSPACE(useq, useq, GtUchar, maxulen);
    }
    if (vlen > maxvlen)
    {
      maxvlen = vlen;
      ALLOCASSIGNSPACE(vseq, vseq, GtUchar, maxvlen);
    }

    gt_encodedsequence_extract_substring(encseq, useq,
                                         boundaries->leftLTR_5,
                                         boundaries->leftLTR_3);
    gt_encodedsequence_extract_substring(encseq, vseq,
                                         boundaries->rightLTR_5,
                                         boundaries->rightLTR_3);
    edist = greedyunitedist(useq,(unsigned long) ulen, /*Implement for encseq */
                            vseq,(unsigned long) vlen);

    /* determine similarity */
    boundaries->similarity = 100.0 *
                             (1 - (((double) edist)/(MAX(ulen,vlen))));

    if (gt_double_smaller_double(boundaries->similarity,
                                 lo->similaritythreshold))
    /* if ( boundaries->similarity < lo->similaritythreshold ) */
    {
      /* delete this LTR-pair candidate */
      arrayLTRboundaries->nextfreeLTRboundaries--;
    }

  }
  FREESPACE(useq);
  FREESPACE(vseq);

  return haserr ? -1 : 0;
}
