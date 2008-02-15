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

#include "libgtcore/arraydef.h"
#include "libgtcore/error.h"
#include "libgtcore/minmax.h"
#include "libgtcore/unused.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/spacedef.h"
#include "libgtmatch/intcode-def.h"
#include "searchforLTRs.h"
#include "repeattypes.h"
#include "repeats.h"
#include "myxdrop.h"
#include "searchTSDandmotif.h"

#include "libgtmatch/greedyedist.pr"
#include "libgtmatch/pos2seqnum.pr"

/*
 The following function checks, if the remaining candidate pairs still
 fulfill the length and distance constraints of LTRs.
 */

static int checklengthanddistanceconstraints(LTRboundaries *boundaries,
                                             RepeatInfo *repeatinfo)
{
  Seqpos ulen, vlen, dist_between_LTRs;

  ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1;
  vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
  dist_between_LTRs = boundaries->rightLTR_5 - boundaries->leftLTR_5;
  if (ulen > (Seqpos)repeatinfo->lmax || vlen > (Seqpos)repeatinfo->lmax ||
     ulen < (Seqpos)repeatinfo->lmin || vlen < (Seqpos)repeatinfo->lmin ||
     dist_between_LTRs > (Seqpos)repeatinfo->dmax ||
     dist_between_LTRs < (Seqpos)repeatinfo->dmin ||
     boundaries->leftLTR_3 >= boundaries->rightLTR_5)
  {
    boundaries->lengthdistconstraint = false;
    boundaries->similarity = 0.0;
    return -1;
  }
  else
  {
    boundaries->lengthdistconstraint = true;
  }

  return 0;
}

static void adjustboundariesfromXdropextension(Myxdropbest xdropbest_left,
                                               Myxdropbest xdropbest_right,
                                               Seqpos seed1_startpos,
                                               Seqpos seed2_startpos,
                                               Seqpos seed1_endpos,
                                               Seqpos seed2_endpos,
                                               LTRboundaries *boundaries,
                                               UNUSED Seqpos offset)
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

#define INITBOUNDARIES(B)\
        (B)->contignumber = repeatptr->contignumber;\
        (B)->leftLTR_5 = (Seqpos)0;\
        (B)->leftLTR_3 = (Seqpos)0;\
        (B)->rightLTR_5 = (Seqpos)0;\
        (B)->rightLTR_3 = (Seqpos)0;\
        (B)->lenleftTSD = (Seqpos)0;\
        (B)->lenrightTSD = (Seqpos)0;\
        (B)->tsd = false;\
        (B)->motif_near_tsd = false;\
        (B)->motif_far_tsd = false;\
        (B)->skipped = false;\
        (B)->similarity = 0.0;

/*
 The following function applies the filter algorithms one after another
 to all candidate pairs.
*/
int searchforLTRs(Sequentialsuffixarrayreader *ssar, LTRharvestoptions *lo,
                  const Seqpos *markpos, Error *err)
{
  unsigned long repeatcounter;
  ArrayMyfrontvalue fronts;
  Myxdropbest xdropbest_left;
  Myxdropbest xdropbest_right;
  Seqpos alilen = 0,
         totallength,
         ulen,
         vlen;
  Uchar *useq = NULL,
        *vseq = NULL;
  unsigned long edist;
  Repeat *repeatptr;
  LTRboundaries *boundaries;
  Seqpos offset = 0;
  const Encodedsequence *encseq =
          encseqSequentialsuffixarrayreader(ssar);

  error_check(err);

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
    alilen = ((Seqpos)lo->repeatinfo.lmax) - repeatptr->len;

    /**** left (reverse) xdrop alignment ****/
    INITARRAY (&fronts, Myfrontvalue);
    if (alilen <= repeatptr->pos1)
    {
      evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) alilen,
                               (int) alilen,
                               (Xdropscore)lo->xdropbelowscore,
                               err);
    }
    else /* do not align over left sequence boundary */
    {
      evalxdroparbitscoresleft(&lo->arbitscores,
                               &xdropbest_left,
                               &fronts,
                               encseq,
                               encseq,
                               repeatptr->pos1,
                               repeatptr->pos1 + repeatptr->offset,
                               (int) repeatptr->pos1,
                               (int) (repeatptr->pos1 + repeatptr->offset),
                               (Xdropscore)lo->xdropbelowscore,
                               err);
    }
    FREEARRAY (&fronts, Myfrontvalue);

    /**** right xdrop alignment ****/
    INITARRAY (&fronts, Myfrontvalue);
    totallength = getencseqtotallength(encseq);
    if (alilen <= totallength - (repeatptr->pos1 + repeatptr->offset +
                                repeatptr->len) )
    {
      evalxdroparbitscoresright (&lo->arbitscores,
                                 &xdropbest_right,
                                 &fronts,
                                 encseq,
                                 encseq,
                                 repeatptr->pos1 + repeatptr->len,
                                 repeatptr->pos1 + repeatptr->offset +
                                 repeatptr->len,
                                 (int) alilen,
                                 (int) alilen,
                                 lo->xdropbelowscore,
                                 err);
    }
    else /* do not align over right sequence boundary */
    {
      evalxdroparbitscoresright(&lo->arbitscores,
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
                                lo->xdropbelowscore,
                                err);
    }
    FREEARRAY (&fronts, Myfrontvalue);

    GETNEXTFREEINARRAY(boundaries,
                       &lo->arrayLTRboundaries,
                       LTRboundaries,
                       5);
    INITBOUNDARIES(boundaries);

    if ( boundaries->contignumber == 0)
    {
      offset = 0;
    }
    else
    {
      offset = markpos[boundaries->contignumber-1];
    }

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
                   boundaries,
                   offset);

    /* if search for motif and/or TSD */
    if ( lo->motif.allowedmismatches < 4U || lo->minlengthTSD > 1U)
    {
      if ( findcorrectboundaries(lo, boundaries, ssar,
                                markpos, err) != 0 )
      {
        return -1;
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
          lo->arrayLTRboundaries.nextfreeLTRboundaries--;
          continue;
        }
      }
    }

    /* check length and distance constraints again */
    if (checklengthanddistanceconstraints(boundaries, &lo->repeatinfo)) {
      /* delete this LTR-pair candidate */
      lo->arrayLTRboundaries.nextfreeLTRboundaries--;
      continue;
    }

    /* check similarity from candidate pair
       copy LTR sequences for greedyunitedist function */
    ulen = boundaries->leftLTR_3 - boundaries->leftLTR_5 + 1;
    vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
    ALLOCASSIGNSPACE(useq, NULL, Uchar, ulen);
    ALLOCASSIGNSPACE(vseq, NULL, Uchar, vlen);

    encseqextract(useq,encseq,boundaries->leftLTR_5,boundaries->leftLTR_3);
    encseqextract(vseq,encseq,boundaries->rightLTR_5,boundaries->rightLTR_3);
    edist = greedyunitedist(useq,(unsigned long) ulen, /*Implement for encseq */
                            vseq,(unsigned long) vlen);

    /* determine similarity */
    boundaries->similarity = 100.0 *
    (1 - (((double) edist)/(MAX(ulen,vlen))));

    if ( boundaries->similarity < lo->similaritythreshold )
    {
      /* delete this LTR-pair candidate */
      lo->arrayLTRboundaries.nextfreeLTRboundaries--;
    }

    FREESPACE(useq);
    FREESPACE(vseq);
  }

  return 0;
}
