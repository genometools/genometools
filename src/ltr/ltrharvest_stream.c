/*
  Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/array_api.h"
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/log.h"
#include "core/mathsupport.h"
#include "core/md5_seqid.h"
#include "core/minmax.h"
#include "core/multithread_api.h"
#include "core/str_api.h"
#include "core/thread_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "match/esa-seqread.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/greedyedist.h"
#include "match/xdrop.h"
#include "ltr/ltrharvest_stream.h"

#define GT_LTRHARVEST_NAME "LTRharvest"

typedef struct
{
  GtUword pos1,         /* first position of maximal repeat (seed) */
          offset,       /* second position = pos1 + offset */
          len,          /* length of maximal repeat  */
          contignumber; /* number of contig for this repeat */
} Repeat;

GT_DECLAREARRAYSTRUCT(Repeat);

/* The datatype RepeatInfo stores all maximal repeats (seeds) and */
/* information about the length and distance constraints. */
typedef struct
{
  GtArrayRepeat repeats; /* array of maximal repeats (seeds) */
  GtUword       lmin,    /* minimum allowed length of a LTR */
                lmax,    /* maximum allowed length of a LTR */
                dmin,    /* minimum distance between LTRs */
                dmax;    /* maximum distance between LTRs,
                            this values is determined as the minimum
                            of the corresponding option value and the
                            maximum length of any sequence in the encseq */
  GtRange       ltrsearchseqrange; /* if start and end are 0, then no range */
} RepeatInfo;

/* The datatype SubRepeatInfo stores information about the maximal repeats */
/* for the TSD detection. */
typedef struct
{
  GtArrayRepeat repeats; /* array of maximal repeats for TSDs */
  GtUword tsd_lmin,    /* minimal length of TSD */
                tsd_lmax,    /* maximal length of TSD */
                offset1, /* offset1 for absolute position 1 in sequence */
                offset2; /* offset2 for absolute position 2 in sequence */
                         /* pos1 < pos2 */
} SubRepeatInfo;

typedef struct
{
  GtUword contignumber,   /* ordinal number of sequence in encseq */
                leftLTR_5,      /* 5' boundary of left LTR */
                leftLTR_3,      /* 3' boundary of left LTR */
                rightLTR_5,     /* 5' boundary of right LTR */
                rightLTR_3,     /* 3' boundary of right LTR */
                lenleftTSD,
                lenrightTSD;
  bool          tsd,            /* If true, then TSDs exist. */
                motif_near_tsd, /* If true, then motif near the TSD exists. */
                motif_far_tsd,  /* If true, then motif at the inner borders
                                   of LTRs exist. */
                lengthdistconstraint, /* If true, length and distance
                                         constraints are satisfied */
                skipped; /* if skipped then because of an overlap
                            with a higher similarity prediction or
                            because of "overlap=no" option */
  double similarity;     /* similarity value of LTRs */
} LTRboundaries;

GT_DECLAREARRAYSTRUCT(LTRboundaries);

typedef enum {
  GT_LTRHARVEST_STREAM_STATE_START,
  GT_LTRHARVEST_STREAM_STATE_REGIONS,
  GT_LTRHARVEST_STREAM_STATE_COMMENTS,
  GT_LTRHARVEST_STREAM_STATE_FEATURES
} GtLTRharvestStreamState;

struct GtLTRharvestStream
{
  const GtNodeStream parent_instance;
  RepeatInfo repeatinfo;
  GtStr *str_indexname;
  GtUword minseedlength;
  double similaritythreshold;
  int xdropbelowscore;
  GtXdropArbitraryscores arbitscores;
  unsigned int allowedmismatches;
  GtLTRFourCharMotif *motif;
  GtUword offset;
  unsigned int minlengthTSD,
               maxlengthTSD;
  GtUword numofboundaries;
  GtUword vicinityforcorrectboundaries;
  GtUword prevseqnum;
  const LTRboundaries **bdptrtab;
  GtArrayLTRboundaries arrayLTRboundaries;
  const GtEncseq *encseq;
  Sequentialsuffixarrayreader *ssar;
  bool verbosemode,
       nooverlaps,
       bestoverlaps,
       output_md5,
       output_seqids;
  GtUword cur_elem_index;
  GtLTRharvestStreamState state;
};

#define gt_ltrharvest_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrharvest_stream_class(), GS);

static int bdcompare(const void *a, const void *b)
{
  const LTRboundaries *bda, *bdb;

  bda = (const LTRboundaries *) a;
  bdb = (const LTRboundaries *) b;
  if (bda == bdb) return 0;
  if (bda->contignumber < bdb->contignumber)
  {
    return -1;
  }
  if (bda->contignumber > bdb->contignumber)
  {
    return 1;
  }
  gt_assert(bda->contignumber == bdb->contignumber);
  if (bda->leftLTR_5 < bdb->leftLTR_5)
  {
    return -1;
  }
  if (bda->leftLTR_5 > bdb->leftLTR_5)
  {
    return 1;
  }
  gt_assert(bda->leftLTR_5 == bdb->leftLTR_5);
  if (bda->rightLTR_3 < bdb->rightLTR_3)
  {
    return -1;
  }
  if (bda->rightLTR_3 > bdb->rightLTR_3)
  {
    return 1;
  }
  return 0;
}

static int gt_simpleexactselfmatchstore(void *info,
                                        const GtGenericEncseq *genericencseq,
                                        GtUword len, GtUword pos1,
                                        GtUword pos2,
                                        GT_UNUSED GtError *err)
{
  GtUword distance;
  RepeatInfo *repeatinfo = (RepeatInfo *) info;

  gt_error_check(err);
  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (repeatinfo->ltrsearchseqrange.start > 0 ||
      repeatinfo->ltrsearchseqrange.end > 0)
  {
    if (pos1 < repeatinfo->ltrsearchseqrange.start  ||
        pos2 + len - 1 > repeatinfo->ltrsearchseqrange.end)
    {
      return 0;
    }
  }
  distance = pos2 - pos1;
  /*test maximal length of candidate pair and distance constraints*/
  if (len <= repeatinfo->lmax && repeatinfo->dmin <= distance
                              && distance <= repeatinfo->dmax)
  {
    GtUword seqnum1, seqnum2;

    gt_assert(genericencseq->hasencseq);
    seqnum1 = gt_encseq_seqnum(genericencseq->seqptr.encseq,pos1);
    seqnum2 = gt_encseq_seqnum(genericencseq->seqptr.encseq,pos2);
    if (seqnum1 == seqnum2)
    {
      Repeat *nextfreerepeatptr;

      GT_GETNEXTFREEINARRAY(nextfreerepeatptr,&repeatinfo->repeats,Repeat,32);
      nextfreerepeatptr->pos1 = pos1;
      nextfreerepeatptr->offset = distance;
      nextfreerepeatptr->len = len;
      nextfreerepeatptr->contignumber = seqnum1;
    }
  }
  return 0;
}

static int gt_subsimpleexactselfmatchstore(void *info,
                                            GT_UNUSED const GtEncseq *encseq,
                                            const GtQuerymatch *querymatch,
                                            GT_UNUSED const GtUchar *query,
                                            GT_UNUSED GtUword
                                              query_totallength,
                                            GT_UNUSED GtError *err)
{
  Repeat *nextfreerepeatptr;
  SubRepeatInfo *sri = (SubRepeatInfo *) info;

  GT_GETNEXTFREEINARRAY (nextfreerepeatptr, &sri->repeats, Repeat, 10);
  nextfreerepeatptr->pos1 = sri->offset1 + gt_querymatch_dbstart(querymatch);
  nextfreerepeatptr->offset = sri->offset2 +
                              gt_querymatch_querystart(querymatch) -
                              (sri->offset1 +
                               gt_querymatch_dbstart(querymatch));
  nextfreerepeatptr->len = gt_querymatch_querylen(querymatch);
  return 0;
}

static const LTRboundaries **compactboundaries(GtUword *numofboundaries,
                                               const GtArrayLTRboundaries *ltr)
{
  GtUword countboundaries = 0, nextfill = 0;
  const LTRboundaries *bd, **bdptrtab;

  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      countboundaries++;
    }
  }
  bdptrtab = gt_malloc(sizeof (LTRboundaries *) * countboundaries);
  nextfill = 0;
  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      bdptrtab[nextfill++] = bd;
    }
  }
  *numofboundaries = countboundaries;
  return bdptrtab;
}

typedef struct
{
  unsigned int left, right;
} LTRMotifmismatches;

/* The following function searches for TSDs and/or a specified palindromic
   motif at the 5' border of the left LTR and 3' border of the right LTR.
   Thereby, all maximal repeats from the vicinity are processed one after
   another to find the TSD with the minimum deviation with regard to the
   boundary position from the x-drop alignment. If also a motif is searched,
   a simple motif check at the boundaries of the TSDs is performed. */
static void searchforbestTSDandormotifatborders(const SubRepeatInfo
                                                                 *subrepeatinfo,
                                                const GtLTRharvestStream *lo,
                                                LTRboundaries *boundaries,
                                                LTRMotifmismatches *mismatches)
{
  GtUword motifpos1,
                motifpos2,
                back, forward,
                oldleftLTR_5  = boundaries->leftLTR_5,
                oldrightLTR_3 = boundaries->rightLTR_3,
                difffromoldboundary1 = 0,
                difffromoldboundary2 = 0;
  LTRMotifmismatches tmp_mm;
  unsigned int hitcounter = 0;
  Repeat *rep;

  if (subrepeatinfo->repeats.nextfreeRepeat > 0)
  {
    boundaries->tsd = true;
  }
  boundaries->motif_near_tsd = false;

  for (rep = subrepeatinfo->repeats.spaceRepeat;
       rep < subrepeatinfo->repeats.spaceRepeat +
             subrepeatinfo->repeats.nextfreeRepeat; rep++)
  {
    /* motifpos1 is the first position after the left repeat */
    motifpos1 = rep->pos1 + rep->len;
    /* motifpos2 is two positions before the right repeat */
    motifpos2 = rep->pos1 + rep->offset - 2;

    for (back = 0; back < rep->len - subrepeatinfo->tsd_lmin + 1; back++)
    {
      for (forward = 0;
           forward < rep->len - subrepeatinfo->tsd_lmin + 1 - back;
           forward++)
      {
        tmp_mm.left = tmp_mm.right = 0;
        if (gt_encseq_get_encoded_char(/* Random access */ lo->encseq,
                                       motifpos1 - back, GT_READMODE_FORWARD)
            != lo->motif->firstleft)
        {
          tmp_mm.left++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ lo->encseq,
                                              motifpos1 + 1 - back,
                                              GT_READMODE_FORWARD)
            != lo->motif->secondleft)
        {
          tmp_mm.left++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ lo->encseq,
                                              motifpos2 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif->firstright)
        {
          tmp_mm.right++;
        }
        if (gt_encseq_get_encoded_char(/* Random access */ lo->encseq,
                                              motifpos2 + 1 + forward,
                                              GT_READMODE_FORWARD)
            != lo->motif->secondright)
        {
          tmp_mm.right++;
        }

        if (tmp_mm.left <= lo->motif->allowedmismatches
            && tmp_mm.right <= lo->motif->allowedmismatches)
        {
          const GtUword tsd_len = rep->len - back - forward;

          /* TSD length not too big */
          if (tsd_len <= subrepeatinfo->tsd_lmax)
          {
            if (!boundaries->motif_near_tsd)
            {
              GtUword max, min;

              /* save number of mismatches */
              *mismatches = tmp_mm;

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
            } else
            {
              GtUword max, min, difffromnewboundary1,
                   difffromnewboundary2;

              /* test if hit is nearer to old boundaries than previous hit */
              max = MAX(oldleftLTR_5, motifpos1 - back);
              min = MIN(oldleftLTR_5, motifpos1 - back);
              difffromnewboundary1 = max - min;
              max = MAX(oldrightLTR_3, motifpos2 + 1 + forward);
              min = MIN(oldrightLTR_3, motifpos2 + 1 + forward);
              difffromnewboundary2 = max - min;

              if (difffromnewboundary1 + difffromnewboundary2 <
                  difffromoldboundary1 + difffromoldboundary2)
              {
                /* save number of mismatches */
                *mismatches = tmp_mm;

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

/* The following function searches only for a specified palindromic motif
   at the 5'-border of left LTR and 3'-border of right LTR. */
static void searchformotifonlyborders(const GtLTRharvestStream *lo,
                                      LTRboundaries *boundaries,
                                      GtUword startleftLTR,
                                      GtUword endleftLTR,
                                      GtUword startrightLTR,
                                      GtUword endrightLTR,
                                      LTRMotifmismatches *motifmismatches)
{
  bool motif1 = false,
       motif2 = false;
  LTRMotifmismatches tmp_mm;
  unsigned int motifmismatches_frombestmatch = 0;
  GtUword idx,
         oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary = 0;

  /**** search for left motif around leftLTR_5 ****/
  for (idx = startleftLTR; idx < endleftLTR; idx++)
  {
    tmp_mm.left = 0;
    if (gt_encseq_get_encoded_char(lo->encseq, idx, GT_READMODE_FORWARD)
          != lo->motif->firstleft)
    {
      tmp_mm.left++;
    }
    if (gt_encseq_get_encoded_char(lo->encseq, idx + 1, GT_READMODE_FORWARD)
          != lo->motif->secondleft)
    {
      tmp_mm.left++;
    }
    if (tmp_mm.left + motifmismatches->left <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if (!motif1)
       {
         GtUword max, min;

         motifmismatches_frombestmatch = tmp_mm.left;
         boundaries->leftLTR_5 = idx;
         motif1 = true;
         max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
         min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
         difffromoldboundary = max - min;
       } else /* next hit */
       {
         GtUword maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_5, idx);
         minval = MIN(oldleftLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if (difffromnewboundary < difffromoldboundary)
         {
           motifmismatches_frombestmatch = tmp_mm.left;
           boundaries->leftLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  motifmismatches->left += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  for (idx = startrightLTR + 1; idx <= endrightLTR; idx++)
  {
    tmp_mm.right = 0;
    if (gt_encseq_get_encoded_char(lo->encseq, idx, GT_READMODE_FORWARD)
            != lo->motif->secondright)
    {
      tmp_mm.right++;
    }
    if (gt_encseq_get_encoded_char(lo->encseq, idx - 1, GT_READMODE_FORWARD)
            != lo->motif->firstright)
    {
      tmp_mm.right++;
    }
    if (tmp_mm.right + motifmismatches->right <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if (!motif2)
       {
         GtUword max, min;

         motifmismatches_frombestmatch = tmp_mm.right;
         boundaries->rightLTR_3 = idx;
         motif2 = true;
         max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
         min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
         difffromoldboundary = max - min;
       } else /* next hit */
       {
         GtUword maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_3, idx);
         minval = MIN(oldrightLTR_3, idx);
         difffromnewboundary = maxval - minval;
         if (difffromnewboundary < difffromoldboundary)
         {
           motifmismatches_frombestmatch = tmp_mm.right;
           boundaries->rightLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  motifmismatches->right += motifmismatches_frombestmatch;
  if (motif1 && motif2)
  {
    boundaries->motif_near_tsd = true;
  } else
  {
    boundaries->motif_near_tsd = false;
  }
}

/* The following function searches for a specified palindromic motif at the
   3' border of left LTR and the 5'border of right LTR. */
static void searchformotifonlyinside(const GtLTRharvestStream *lo,
                                     LTRboundaries *boundaries,
                                     LTRMotifmismatches *motifmismatches)
{
  bool motif1 = false,
       motif2 = false;
  GtUword startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         oldleftLTR_3  = boundaries->leftLTR_3,
         oldrightLTR_5 = boundaries->rightLTR_5,
         difffromoldboundary = 0,
         idx;
  LTRMotifmismatches tmp_mm;
  unsigned int motifmismatches_frombestmatch = 0;

  /* Search vicinity of 3' border of left LTR.
     Do not align over 5' border of left LTR,
     in case of need decrease alignment length */
  if ((boundaries->leftLTR_3 < lo->vicinityforcorrectboundaries)
      || (startleftLTR = boundaries->leftLTR_3 -
         lo->vicinityforcorrectboundaries + 1) <
      boundaries->leftLTR_5 + 2)
  {
    startleftLTR = boundaries->leftLTR_5 + 2;
  }
  /* do not align over 5' border of right LTR */
  if ((endleftLTR = boundaries->leftLTR_3 +
       lo->vicinityforcorrectboundaries - 1) >
      boundaries->rightLTR_5 - 1)
  {
    endleftLTR = boundaries->rightLTR_5 - 1;
  }
  /* Search vicinity of 5' border of right LTR.
     Do not align over 3' border of left LTR */
  if ((startrightLTR = boundaries->rightLTR_5 -
         lo->vicinityforcorrectboundaries + 1)
       < boundaries->leftLTR_3 + 1)
  {
    startrightLTR = boundaries->leftLTR_3 + 1;
  }
  /* do not align over 3' border of right LTR */
  if ((endrightLTR = boundaries->rightLTR_5 +
       lo->vicinityforcorrectboundaries - 1) >
      boundaries->rightLTR_3 - 2)
  {
    endrightLTR = boundaries->rightLTR_3 - 2;
  }

  /* Search for right motif around leftLTR_3 */
  for (idx = startleftLTR + 1; idx <= endleftLTR; idx++)
  {
    tmp_mm.left = 0;
    if (gt_encseq_get_encoded_char(lo->encseq, idx, GT_READMODE_FORWARD)
              != lo->motif->secondright)
    {
      tmp_mm.left++;
    }
    if (gt_encseq_get_encoded_char(lo->encseq, idx - 1, GT_READMODE_FORWARD)
              != lo->motif->firstright)
    {
      tmp_mm.left++;
    }
    if (tmp_mm.left + motifmismatches->left <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if (!motif1)
       {
         GtUword maxval, minval;

         motifmismatches_frombestmatch = tmp_mm.left;
         boundaries->leftLTR_3 = idx;
         motif1 = true;
         maxval = MAX(oldleftLTR_3, boundaries->leftLTR_3);
         minval = MIN(oldleftLTR_3, boundaries->leftLTR_3);
         difffromoldboundary = maxval - minval;
       } else /* next hit */
       {
         GtUword maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldleftLTR_3, idx);
         minval = MIN(oldleftLTR_3, idx);
         difffromnewboundary = maxval - minval;

         if (difffromnewboundary < difffromoldboundary)
         {
           motifmismatches_frombestmatch = tmp_mm.left;
           boundaries->leftLTR_3 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  motifmismatches->left += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  /* Search for left motif around rightLTR_5 */
  for (idx = startrightLTR ; idx < endrightLTR; idx++)
  {
    tmp_mm.right = 0;
    if (gt_encseq_get_encoded_char(lo->encseq, idx, GT_READMODE_FORWARD)
               != lo->motif->firstleft)
    {
      tmp_mm.right++;
    }
    if (gt_encseq_get_encoded_char(lo->encseq, idx + 1, GT_READMODE_FORWARD)
               != lo->motif->secondleft)
    {
      tmp_mm.right++;
    }
    if (tmp_mm.right + motifmismatches->right <= lo->motif->allowedmismatches)
    {
       /* first hit */
       if (!motif2)
       {
         GtUword maxval, minval;

         motifmismatches_frombestmatch = tmp_mm.right;
         boundaries->rightLTR_5 = idx;
         motif2 = true;
         maxval = MAX(oldrightLTR_5, boundaries->rightLTR_5);
         minval = MIN(oldrightLTR_5, boundaries->rightLTR_5);
         difffromoldboundary = maxval - minval;
       } else /* next hit */
       {
         GtUword maxval, minval, difffromnewboundary;

         /* test if hit is nearer to old boundaries than previous hit */
         maxval = MAX(oldrightLTR_5, idx);
         minval = MIN(oldrightLTR_5, idx);
         difffromnewboundary = maxval - minval;

         if (difffromnewboundary < difffromoldboundary)
         {
           motifmismatches_frombestmatch = tmp_mm.right;
           boundaries->rightLTR_5 = idx;
           difffromoldboundary = difffromnewboundary;
         }
       }
    }
  }
  motifmismatches->right += motifmismatches_frombestmatch;
  if (motif1 && motif2)
  {
    boundaries->motif_far_tsd = true;
  } else
  {
    boundaries->motif_far_tsd = false;
  }
}

/* The following function searches for TSDs and/or a specified palindromic motif
   at the 5' border of left LTR and 3' border of right LTR. */

static int searchforTSDandorMotifoutside(const GtLTRharvestStream *lo,
                                         LTRboundaries *boundaries,
                                         LTRMotifmismatches *motifmismatches,
                                         GtError *err)
{
  GtUword startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         leftlen,
         rightlen,
         sequenceendpos,
         seqstartpos,
         seqlength;
  SubRepeatInfo subrepeatinfo;
  bool haserr = false;

  gt_error_check(err);

  /* check border cases */
  /* vicinity of 5'-border of left LTR */
  seqstartpos = gt_encseq_seqstartpos(lo->encseq, boundaries->contignumber);
  seqlength = gt_encseq_seqlength(lo->encseq, boundaries->contignumber);
  if (boundaries->contignumber == 0)
  {
    /* do not align over left sequence boundary,
       in case of need decrease alignment length */
    if (boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries)
    {
      startleftLTR = seqstartpos;
    } else
    {
      startleftLTR = boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries;
    }
  } else
  {
    /* do not align over left separator symbol
       in case of need decrease alignment length */
    if (boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries)
    {
      startleftLTR = seqstartpos;
    } else
    {
      startleftLTR = boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries;
      if (startleftLTR < seqstartpos && boundaries->leftLTR_5 >= seqstartpos)
      {
        startleftLTR = seqstartpos;
      }
    }
  }
  /* do not align over 3'-border of left LTR */
  endleftLTR = boundaries->leftLTR_5 + lo->vicinityforcorrectboundaries;
  if (endleftLTR > boundaries->leftLTR_3 - 2) /* -2 because of possible motif */
  {
    endleftLTR = boundaries->leftLTR_3 - 2;
  }
  leftlen = endleftLTR - startleftLTR + 1;

  /* vicinity of 3'-border of right LTR
     do not align over 5'border of right LTR */
  startrightLTR = boundaries->rightLTR_3 - lo->vicinityforcorrectboundaries;
  if (startrightLTR < boundaries->rightLTR_5 + 2) /* +2 for possible motif */
  {
    startrightLTR = boundaries->rightLTR_5 + 2;
  }
  sequenceendpos = seqstartpos + seqlength - 1;
  /* do not align into next sequence in case of need decrease alignment
     length */
  endrightLTR = boundaries->rightLTR_3 + lo->vicinityforcorrectboundaries;
  if (endrightLTR > sequenceendpos && boundaries->rightLTR_3 <= sequenceendpos)
  {
    endrightLTR = sequenceendpos;
  }
  rightlen = endrightLTR - startrightLTR + 1;

  /* now, search for correct boundaries */

  /* search for TSDs and/or motif */
  if (lo->minlengthTSD > 1U)
  {
    GtUchar *dbseq = gt_malloc(sizeof (*dbseq) * leftlen),
            *query = gt_malloc(sizeof (*query) * rightlen);

    gt_encseq_extract_encoded(lo->encseq,dbseq,startleftLTR,endleftLTR);
    gt_encseq_extract_encoded(lo->encseq,query,startrightLTR,endrightLTR);
    GT_INITARRAY(&subrepeatinfo.repeats, Repeat);
    subrepeatinfo.tsd_lmin = (GtUword) lo->minlengthTSD;
    subrepeatinfo.tsd_lmax = (GtUword) lo->maxlengthTSD;
    gt_assert(startleftLTR < startrightLTR);
    subrepeatinfo.offset1 = startleftLTR;
    subrepeatinfo.offset2 = startrightLTR;

    if (gt_sarrquerysubstringmatch(dbseq,
                                   leftlen,
                                   query,
                                   rightlen,
                                   lo->minlengthTSD,
                                   gt_encseq_alphabet(lo->encseq),
                                   gt_subsimpleexactselfmatchstore,
                                   &subrepeatinfo,
                                   NULL,
                                   err) != 0)
    {
       haserr = true;
    }
    gt_free(dbseq);
    gt_free(query);

    if (!haserr)
    {
      searchforbestTSDandormotifatborders(&subrepeatinfo,
                                          lo,
                                          boundaries,
                                          motifmismatches);
    }
    GT_FREEARRAY (&subrepeatinfo.repeats, Repeat);
  } else /* no search for TSDs, search for motif only */
  {
    searchformotifonlyborders(lo,
                              boundaries,
                              startleftLTR,
                              endleftLTR,
                              startrightLTR,
                              endrightLTR,
                              motifmismatches);
  }
  return haserr ? -1 : 0;
}

/* The following function searches for TSD and/or a specified palindromic motif
   at the borders of left LTR and the right LTR, respectively. */
static int gt_findcorrectboundaries(const GtLTRharvestStream *lo,
                                    LTRboundaries *boundaries,
                                    GtError *err)
{
  LTRMotifmismatches motifmismatches;
  gt_error_check(err);

  motifmismatches.left = motifmismatches.right = 0;
  /* first: 5'-border of left LTR and 3'-border of right LTR */
  if (searchforTSDandorMotifoutside(lo, boundaries, &motifmismatches, err) != 0)
  {
    return -1;
  }
  /* second: 3'-border of left LTR and 5'-border of right LTR */
  if (lo->motif->allowedmismatches < 4U)
  {
    searchformotifonlyinside(lo,boundaries,&motifmismatches);
  }
  return 0;
}

static bool checklengthanddistanceconstraints(LTRboundaries *boundaries,
                                              const RepeatInfo *repeatinfo)
{
  GtUword ulen, vlen, dist_between_LTRs;

  ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1;
  vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
  dist_between_LTRs = boundaries->rightLTR_5 - boundaries->leftLTR_5;
  if (ulen > repeatinfo->lmax
        || vlen > repeatinfo->lmax
        || ulen < repeatinfo->lmin
        || vlen < repeatinfo->lmin
        || dist_between_LTRs > repeatinfo->dmax
        || dist_between_LTRs < repeatinfo->dmin
        || boundaries->leftLTR_3 >= boundaries->rightLTR_5)
  {
    boundaries->lengthdistconstraint = false;
    boundaries->similarity = 0.0;
    return false;
  } else
  {
    boundaries->lengthdistconstraint = true;
    return true;
  }
}

static void adjustboundariesfromXdropextension(GtXdropbest xdropbest_left,
                                               GtXdropbest xdropbest_right,
                                               GtUword seed1_startpos,
                                               GtUword seed2_startpos,
                                               GtUword seed1_endpos,
                                               GtUword seed2_endpos,
                                               LTRboundaries *boundaries)
{
  /* left alignment */
  boundaries->leftLTR_5  = seed1_startpos - xdropbest_left.ivalue;
  boundaries->rightLTR_5 = seed2_startpos - xdropbest_left.jvalue;

  /* right alignment */
  boundaries->leftLTR_3  = seed1_endpos + xdropbest_right.ivalue;
  boundaries->rightLTR_3 = seed2_endpos + xdropbest_right.jvalue;
}

/* The following function applies the filter algorithms one after another
   to all candidate pairs */
static int gt_searchforLTRs(GtLTRharvestStream *lo,
                            GtArrayLTRboundaries *arrayLTRboundaries,
                            GT_UNUSED GtMutex *rmutex,
                            GT_UNUSED GtMutex *wmutex,
                            GtUword *cur_seed,
                            GtError *err)
{
  GtUword my_seed;
  GtXdropresources *xdropresources;
  GtXdropbest xdropbest_left, xdropbest_right;
#undef GT_GREEDY_BUFFER
#ifdef GT_GREEDY_BUFFER
  GtUword maxulen = 0, maxvlen = 0;
  GtUchar *useq = NULL,
          *vseq = NULL;
#endif
  GtSeqabstract *sa_useq = gt_seqabstract_new_empty(),
                *sa_vseq = gt_seqabstract_new_empty();
  GtUword edist,
                alilen = 0;
  Repeat *repeatptr;
  LTRboundaries boundaries, *boundaries_ptr;
  GtFrontResource *frontresource = gt_frontresource_new(100UL);
  bool haserr = false;

  gt_error_check(err);
  xdropresources = gt_xdrop_resources_new(&lo->arbitscores);

  /* XXX: do thread-synchronized error checking */
  while (true) {
    GtUword ulen,
                  vlen,
                  seqend,
                  seqstart;
    gt_mutex_lock(rmutex);
    if (*cur_seed < lo->repeatinfo.repeats.nextfreeRepeat)
      my_seed = (*cur_seed)++;
    else {
      gt_mutex_unlock(rmutex);
      break;
    }
    gt_mutex_unlock(rmutex);

    repeatptr = &(lo->repeatinfo.repeats.spaceRepeat[my_seed]);
    seqstart = gt_encseq_seqstartpos(lo->encseq, repeatptr->contignumber);
    seqend = seqstart + gt_encseq_seqlength(lo->encseq,repeatptr->contignumber);
    gt_assert(lo->repeatinfo.lmax >= repeatptr->len);
    alilen = lo->repeatinfo.lmax - repeatptr->len;

    /**** left (reverse) xdrop alignment ****/
    if (seqstart < repeatptr->pos1)
    {
      /* excluding pos1 and pos1+offset from the alignment */
      gt_assert(seqstart <= repeatptr->pos1);
      if (alilen <= repeatptr->pos1 - seqstart)
      {
        gt_seqabstract_reinit_encseq(sa_useq,
                                     lo->encseq,
                                     alilen,
                                     repeatptr->pos1 - (alilen));
        gt_seqabstract_reinit_encseq(sa_vseq,
                                     lo->encseq,
                                     alilen,
                                     (repeatptr->pos1 +
                                      repeatptr->offset) -
                                       alilen);
      } else
      {
        gt_seqabstract_reinit_encseq(sa_useq,
                                     lo->encseq,
                                     repeatptr->pos1 - seqstart,
                                     seqstart);
        gt_seqabstract_reinit_encseq(sa_vseq,
                                     lo->encseq,
                                     (repeatptr->pos1 +
                                       repeatptr->offset) -
                                       seqstart,
                                     seqstart + repeatptr->offset);
      }
      gt_evalxdroparbitscoresextend(false,
                                    &xdropbest_left,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                    (GtXdropscore) lo->xdropbelowscore);
    } else
    {
      xdropbest_left.ivalue = 0;
      xdropbest_left.jvalue = 0;
      xdropbest_left.score = 0;
    }
    /**** right (forward) xdrop alignment ****/
    if (repeatptr->pos1 + repeatptr->offset + repeatptr->len < seqend)
    {
      if (alilen <= seqend - (repeatptr->pos1 +
                              repeatptr->offset +
                              repeatptr->len))
      {
        gt_seqabstract_reinit_encseq(sa_useq,
                                     lo->encseq,
                                     alilen,
                                     repeatptr->pos1 +
                                       repeatptr->len);
        gt_seqabstract_reinit_encseq(sa_vseq,
                                     lo->encseq,
                                     alilen,
                                     repeatptr->pos1 +
                                       repeatptr->offset +
                                       repeatptr->len);
      } else
      {
        gt_seqabstract_reinit_encseq(sa_useq,lo->encseq,
                                     seqend - (repeatptr->pos1 +
                                               repeatptr->offset +
                                               repeatptr->len),
                                     repeatptr->pos1 + repeatptr->len);
        gt_seqabstract_reinit_encseq(sa_vseq,lo->encseq,
                                     seqend - (repeatptr->pos1 +
                                               repeatptr->offset +
                                               repeatptr->len),
                                     repeatptr->pos1 +
                                       repeatptr->offset +
                                       repeatptr->len);
      }
      gt_evalxdroparbitscoresextend(true,
                                    &xdropbest_right,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                    (GtXdropscore) lo->xdropbelowscore);
    } else
    {
      xdropbest_right.ivalue = 0;
      xdropbest_right.jvalue = 0;
      xdropbest_right.score = 0;
    }

    boundaries.contignumber = repeatptr->contignumber;
    boundaries.leftLTR_5 = (GtUword) 0;
    boundaries.leftLTR_3 = (GtUword) 0;
    boundaries.rightLTR_5 = (GtUword) 0;
    boundaries.rightLTR_3 = (GtUword) 0;
    boundaries.lenleftTSD = (GtUword) 0;
    boundaries.lenrightTSD = (GtUword) 0;
    boundaries.tsd = false;
    boundaries.motif_near_tsd = false;
    boundaries.motif_far_tsd = false;
    boundaries.skipped = false;
    boundaries.similarity = 0.0;

    /* store new boundaries-positions in boundaries */
    adjustboundariesfromXdropextension(
                   xdropbest_left,
                   xdropbest_right,
                   repeatptr->pos1,  /*seed1 startpos*/
                   repeatptr->pos1 + repeatptr->offset, /*seed2 startpos*/
                   repeatptr->pos1 + repeatptr->len - 1, /*seed1 endpos*/
                   repeatptr->pos1 + repeatptr->offset + repeatptr->len - 1,
                                                           /*seed2 endpos*/
                   &boundaries);

    /* if search for motif and/or TSD */
    if (lo->motif->allowedmismatches < 4U || lo->minlengthTSD > 1U)
    {
      if (gt_findcorrectboundaries(lo, &boundaries, err) != 0)
      {
        haserr = true;
        break;
      }

      /* if search for TSDs and (not) motif */
      if (boundaries.tsd &&
          (lo->motif->allowedmismatches >= 4U ||
          (boundaries.motif_near_tsd && boundaries.motif_far_tsd)))
      {
        /* predicted as full LTR-pair, keep it */
      } else
      {
        /* if search for motif only (and not TSD) */
        if (lo->minlengthTSD <= 1U &&
            boundaries.motif_near_tsd &&
            boundaries.motif_far_tsd)
        {
          /* predicted as full LTR-pair, keep it */
        } else
        {
          continue;
        }
      }
    }

    /* check length and distance constraints again */
    if (!checklengthanddistanceconstraints(&boundaries, &lo->repeatinfo))
    {
      continue;
    }

    /* check similarity from candidate pair
       copy LTR sequences for greedyunitedist function */
    ulen = boundaries.leftLTR_3 - boundaries.leftLTR_5 + 1;
    vlen = boundaries.rightLTR_3 - boundaries.rightLTR_5 + 1;
#ifdef GT_GREEDY_BUFFER
    if (ulen > maxulen)
    {
      maxulen = ulen;
      useq = gt_realloc(useq, sizeof (*useq) * maxulen);
    }
    if (vlen > maxvlen)
    {
      maxvlen = vlen;
      vseq = gt_realloc(vseq, sizeof (*vseq) * maxvlen);
    }
    gt_encseq_extract_encoded(encseq, useq, boundaries.leftLTR_5,
                                            boundaries.leftLTR_3);
    gt_encseq_extract_encoded(encseq, vseq, boundaries.rightLTR_5,
                                            boundaries.rightLTR_3);
    gt_seqabstract_reinit_gtuchar(sa_useq, useq, ulen, 0);
    gt_seqabstract_reinit_gtuchar(sa_vseq, vseq, vlen, 0);
#else
    gt_seqabstract_reinit_encseq(sa_useq, lo->encseq, ulen,
                                 boundaries.leftLTR_5);
    gt_seqabstract_reinit_encseq(sa_vseq, lo->encseq, vlen,
                                 boundaries.rightLTR_5);
#endif
    edist = greedyunitedist(frontresource,sa_useq,sa_vseq);

    /* determine similarity */
    boundaries.similarity = 100.0 * (1.0 - (double) edist/MAX(ulen,vlen));

    if (!gt_double_smaller_double(boundaries.similarity,
                                  lo->similaritythreshold))
    {
      gt_mutex_lock(wmutex);
      GT_GETNEXTFREEINARRAY(boundaries_ptr,arrayLTRboundaries,LTRboundaries,5);
      *boundaries_ptr = boundaries;
      gt_mutex_unlock(wmutex);
    }
  }
#ifdef GT_GREEDY_BUFFER
  FREESPACE(useq);
  FREESPACE(vseq);
#endif
  gt_xdrop_resources_delete(xdropresources);
  gt_seqabstract_delete(sa_useq);
  gt_seqabstract_delete(sa_vseq);
  gt_frontresource_delete(frontresource);
  return haserr ? -1 : 0;
}

typedef struct {
  GtLTRharvestStream *lo;
  GtArrayLTRboundaries *arrayLTRboundaries;
  const GtEncseq *encseq;
  GtError *err;
  GtMutex *rmutex, *wmutex;
  GtUword cur_seed;
} GtLTRharvestThreadInfo;

static void* gt_searchforLTRs_threadfunc(void *data) {
  GtLTRharvestThreadInfo *info = (GtLTRharvestThreadInfo*) data;
  GT_UNUSED int rval;
  gt_assert(info);
  rval = gt_searchforLTRs(info->lo, info->arrayLTRboundaries, info->rmutex,
                          info->wmutex, &info->cur_seed, info->err);
  gt_assert(rval == 0);
  return NULL;
}

/* The following function removes exact duplicates from the (sorted!)
   array of predicted LTR elements. Exact duplicates occur when different seeds
   are extended to same boundary coordinates. */
static void gt_removeduplicates(GtArrayLTRboundaries *arrayLTRboundaries)
{
  GtUword i;
  LTRboundaries *oldboundaries,
                *boundaries;

  oldboundaries = arrayLTRboundaries->spaceLTRboundaries;
  for (i = 1UL; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
  {
    boundaries = arrayLTRboundaries->spaceLTRboundaries + i;
    gt_assert(!boundaries->skipped);
    if (oldboundaries->leftLTR_5 == boundaries->leftLTR_5
          && oldboundaries->rightLTR_3 == boundaries->rightLTR_3)
      boundaries->skipped = true;
    else
      oldboundaries = boundaries;
  }
}

static inline bool ltrboundaries_overlap(GtRange *a, LTRboundaries *b) {
  gt_assert(b);
  return (a->start <= b->rightLTR_3 && a->end >= b->leftLTR_5);
}

/* The following function removes overlaps and deletes the prediction with
   a lower similarity value. If "nooverlapallowed" is set, then all
   overlapping predictions are deleted completely. */
static void gt_removeoverlapswithlowersimilarity(GtArrayLTRboundaries
                                                            *arrayLTRboundaries,
                                                 bool nooverlapallowed)
{
  GtUword i;
  LTRboundaries *boundaries,
                *oldboundaries,
                *maxsimboundaries = NULL;
  GtRange refrng;
  gt_assert(arrayLTRboundaries != NULL);

  if (arrayLTRboundaries->spaceLTRboundaries == NULL)
    return;
  gt_assert(arrayLTRboundaries->spaceLTRboundaries != NULL);

  maxsimboundaries = oldboundaries = arrayLTRboundaries->spaceLTRboundaries;
  refrng.start = oldboundaries->leftLTR_5;
  refrng.end = oldboundaries->rightLTR_3;
  for (i = 1UL; i < arrayLTRboundaries->nextfreeLTRboundaries; i++) {
    boundaries = arrayLTRboundaries->spaceLTRboundaries + i;
    if (boundaries->skipped)
      continue;
    if (ltrboundaries_overlap(&refrng, boundaries)) {
      /* due to sortedness of the array, the left border shouldn't change... */
      gt_assert(refrng.start <= boundaries->leftLTR_5);
      /* ...so only update the right border */
      refrng.end = MAX(boundaries->rightLTR_3, refrng.end);
      if (nooverlapallowed) {
        oldboundaries->skipped = true;
        boundaries->skipped = true;
      } else {
        if (gt_double_smaller_double(maxsimboundaries->similarity,
                                     boundaries->similarity)) {
          maxsimboundaries->skipped = true;
          maxsimboundaries = boundaries;
        } else {
          boundaries->skipped = true;
        }
      }
    } else {
      /* no more overlapping with cluster, reset */
      oldboundaries = boundaries;
      refrng.start = boundaries->leftLTR_5;
      refrng.end = boundaries->rightLTR_3;
      maxsimboundaries = boundaries;
    }
  }
}

static int gt_ltrharvest_stream_next(GtNodeStream *ns,
                                     GtGenomeNode **gn,
                                     GtError *err)
{
  GtLTRharvestStream *ltrh_stream;
  GtLTRharvestThreadInfo threadinfo;
  int had_err = 0;
  gt_error_check(err);

  ltrh_stream = gt_ltrharvest_stream_cast(ns);
  if (ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_START) {
    GT_INITARRAY(&ltrh_stream->repeatinfo.repeats, Repeat);
    ltrh_stream->prevseqnum = GT_UNDEF_UWORD;
    if (!had_err && gt_enumeratemaxpairs(ltrh_stream->ssar,
                      (unsigned int) ltrh_stream->minseedlength,
                      gt_simpleexactselfmatchstore,
                      &ltrh_stream->repeatinfo,
                      err) != 0)
    {
      had_err = -1;
    }

    threadinfo.lo = ltrh_stream;
    threadinfo.encseq = ltrh_stream->encseq;
    threadinfo.arrayLTRboundaries = &ltrh_stream->arrayLTRboundaries;
    threadinfo.err = err;
    threadinfo.cur_seed = 0;
    threadinfo.rmutex = gt_mutex_new();
    threadinfo.wmutex = gt_mutex_new();
    /* apply the seed extension and filter algorithms */
    if (!had_err && gt_multithread(gt_searchforLTRs_threadfunc,
                                   &threadinfo, err) != 0)
    {
      had_err = -1;
    }
    gt_mutex_delete(threadinfo.rmutex);
    gt_mutex_delete(threadinfo.wmutex);

    /* not needed any longer */
    GT_FREEARRAY(&ltrh_stream->repeatinfo.repeats, Repeat);

    /* sort results after seed extension */
    if (!had_err && threadinfo.arrayLTRboundaries->spaceLTRboundaries) {
      qsort(threadinfo.arrayLTRboundaries->spaceLTRboundaries,
            (size_t) threadinfo.arrayLTRboundaries->nextfreeLTRboundaries,
             sizeof (LTRboundaries),  bdcompare);
    }

    /* remove exact duplicates */
    if (!had_err) {
      gt_removeduplicates(&ltrh_stream->arrayLTRboundaries);
    }

    /* remove overlapping predictions if desired */
    if (!had_err && (ltrh_stream->nooverlaps || ltrh_stream->bestoverlaps)) {
      gt_removeoverlapswithlowersimilarity(&ltrh_stream->arrayLTRboundaries,
                                            ltrh_stream->nooverlaps);
    }

    /* build array of non-skipped elements */
    if (!had_err) {
      ltrh_stream->bdptrtab = compactboundaries(&ltrh_stream->numofboundaries,
                                              &ltrh_stream->arrayLTRboundaries);
    }

    /* advance stream state to output regions */
    ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_REGIONS;
  }

  /* first stream out the region nodes */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_REGIONS) {
    bool skip = false;
    if (ltrh_stream->cur_elem_index < ltrh_stream->numofboundaries) {
      GtUword seqnum, seqlength;
      GtGenomeNode *rn;
      GtStr *seqid;
      seqnum = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
      if (ltrh_stream->prevseqnum == GT_UNDEF_UWORD)
      {
        ltrh_stream->prevseqnum = seqnum;
      } else
      {
        while (ltrh_stream->prevseqnum == seqnum)
        {
          ltrh_stream->cur_elem_index++;
          if (ltrh_stream->cur_elem_index >= ltrh_stream->numofboundaries)
          {
            skip = true;
            break;
          }
          seqnum
            = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
        }
      }
      if (!skip)
      {
        ltrh_stream->prevseqnum = seqnum;
        seqlength = gt_encseq_seqlength(ltrh_stream->encseq, seqnum);
        seqid = gt_str_new();
        if (gt_encseq_has_md5_support(ltrh_stream->encseq)
              && ltrh_stream->output_md5) {
          GtMD5Tab *md5_tab = NULL;
          md5_tab = gt_encseq_get_md5_tab(ltrh_stream->encseq, err);
          if (!md5_tab) {
            had_err = -1;
            gt_str_delete(seqid);
          }
          if (!had_err) {
            gt_str_append_cstr(seqid, GT_MD5_SEQID_PREFIX);
            gt_str_append_cstr(seqid, gt_md5_tab_get(md5_tab, seqnum));
            gt_str_append_char(seqid, GT_MD5_SEQID_SEPARATOR);
          }
          gt_md5_tab_delete(md5_tab);
        }

        if (!had_err) {
          if (gt_encseq_has_description_support(ltrh_stream->encseq)
                && ltrh_stream->output_seqids) {
            GtUword desclength = 0UL,
                          i = 0UL;
            const char *desc;
            desc = gt_encseq_description(ltrh_stream->encseq,
                                         &desclength, seqnum);
            while (*(desc+i) != ' ' && i != desclength)
              i++;
            gt_str_append_cstr_nt(seqid, desc, i);
          } else {
            gt_str_append_cstr(seqid, "seq");
            gt_str_append_ulong(seqid, seqnum);
          }
          rn = gt_region_node_new(seqid,
                                  1 + (GtUword) ltrh_stream->offset,
                                  seqlength
                                    + (GtUword) ltrh_stream->offset);
          gt_str_delete(seqid);
          *gn = rn;
          ltrh_stream->cur_elem_index++;
        }
      } else {
        ltrh_stream->cur_elem_index = 0;
        if (gt_encseq_has_description_support(ltrh_stream->encseq))
          ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_COMMENTS;
        else
          ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    } else
    {
      ltrh_stream->cur_elem_index = 0;
      if (gt_encseq_has_description_support(ltrh_stream->encseq))
        ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_COMMENTS;
      else
        ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* then stream out the comment nodes */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_COMMENTS)
  {
    bool skip = false;
    if (ltrh_stream->cur_elem_index < ltrh_stream->numofboundaries)
    {
      const char *desc;
      char buf[BUFSIZ];
      GtUword desclen, seqnum;
      GtGenomeNode *cn;
      seqnum = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
      if (ltrh_stream->prevseqnum == GT_UNDEF_UWORD)
      {
        ltrh_stream->prevseqnum = seqnum;
      } else
      {
        while (ltrh_stream->prevseqnum == seqnum)
        {
          ltrh_stream->cur_elem_index++;
          if (ltrh_stream->cur_elem_index >= ltrh_stream->numofboundaries)
          {
            skip = true;
            break;
          }
          seqnum
            = ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]->contignumber;
        }
      }
      if (!skip)
      {
        ltrh_stream->prevseqnum = seqnum;
        desc = gt_encseq_description(ltrh_stream->encseq, &desclen, seqnum);
        (void) strncpy(buf, desc, (size_t) (desclen * sizeof (char)));
        buf[desclen] = '\0';
        cn = gt_comment_node_new(buf);
        *gn = cn;
        ltrh_stream->cur_elem_index++;
      } else
      {
        ltrh_stream->cur_elem_index = 0;
        ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    } else
    {
      ltrh_stream->cur_elem_index = 0;
      ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* finally, stream out the features */
  if (!had_err && ltrh_stream->state == GT_LTRHARVEST_STREAM_STATE_FEATURES)
  {
    if (ltrh_stream->cur_elem_index <  ltrh_stream->numofboundaries)
    {
      GtStr *seqid, *source;
      char buf[BUFSIZ];
      GtGenomeNode *node, *parent;
      const LTRboundaries *elem =
        (ltrh_stream->bdptrtab[ltrh_stream->cur_elem_index]);
      GtUword seqstartpos;

      seqstartpos = gt_encseq_seqstartpos(ltrh_stream->encseq,
                                          elem->contignumber);
      source = gt_str_new_cstr(GT_LTRHARVEST_NAME);
      seqid = gt_str_new();
      if (gt_encseq_has_md5_support(ltrh_stream->encseq)
            && ltrh_stream->output_md5) {
        GtMD5Tab *md5_tab = NULL;
        md5_tab = gt_encseq_get_md5_tab(ltrh_stream->encseq, err);
        if (!md5_tab) {
          had_err = -1;
          gt_str_delete(seqid);
          gt_str_delete(source);
        }
        if (!had_err) {
          gt_str_append_cstr(seqid, GT_MD5_SEQID_PREFIX);
          gt_str_append_cstr(seqid,
                             gt_md5_tab_get(md5_tab, elem->contignumber));
          gt_str_append_char(seqid, GT_MD5_SEQID_SEPARATOR);
        }
        gt_md5_tab_delete(md5_tab);
      }

      if (!had_err && gt_encseq_has_description_support(ltrh_stream->encseq)
            && ltrh_stream->output_seqids) {
        GtUword desclength = 0UL,
                      i = 0UL;
        const char *desc;
        desc = gt_encseq_description(ltrh_stream->encseq,
                                     &desclength, elem->contignumber);
        while (*(desc+i) != ' ' && i != desclength)
          i++;
        gt_str_append_cstr_nt(seqid, desc, i);
      } else {
        gt_str_append_cstr(seqid, "seq");
        gt_str_append_ulong(seqid, elem->contignumber);
      }

      /* repeat_region */
      node = gt_feature_node_new(seqid,
                                 gt_ft_repeat_region,
                                 elem->leftLTR_5 - seqstartpos + 1
                                   - elem->lenleftTSD
                                   + ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + elem->lenrightTSD
                                   + ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;
      if (ltrh_stream->motif->allowedmismatches < 4U)
      {
        node = gt_feature_node_new(seqid,
                                   gt_ft_inverted_repeat,
                                   elem->leftLTR_5 - seqstartpos + 1
                                     + ltrh_stream->offset,
                                   elem->leftLTR_5 - seqstartpos + 2
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   gt_ft_inverted_repeat,
                                   elem->leftLTR_3 - seqstartpos
                                     + ltrh_stream->offset,
                                   elem->leftLTR_3 - seqstartpos + 1
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   gt_ft_inverted_repeat,
                                   elem->rightLTR_5 - seqstartpos + 1
                                     + ltrh_stream->offset,
                                   elem->rightLTR_5 - seqstartpos + 2
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   gt_ft_inverted_repeat,
                                   elem->rightLTR_3 - seqstartpos
                                     + ltrh_stream->offset,
                                   elem->rightLTR_3 - seqstartpos + 1
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
      }
      /* target_site_duplication */
      if (ltrh_stream->minlengthTSD > 1U)
      {
        node = gt_feature_node_new(seqid,
                                   gt_ft_target_site_duplication,
                                   elem->leftLTR_5 - seqstartpos + 1
                                     - elem->lenleftTSD
                                     + ltrh_stream->offset,
                                   elem->leftLTR_5 - seqstartpos
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
        node = gt_feature_node_new(seqid,
                                   gt_ft_target_site_duplication,
                                   elem->rightLTR_3 - seqstartpos + 2
                                     + ltrh_stream->offset,
                                   elem->rightLTR_3 - seqstartpos + 1
                                     + elem->lenrightTSD
                                     + ltrh_stream->offset,
                                   GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                                  (GtFeatureNode*) node);
      }
      /* LTR_retrotransposon */
      node = gt_feature_node_new(seqid,
                                 gt_ft_LTR_retrotransposon,
                                 elem->leftLTR_5 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      (void) snprintf(buf, BUFSIZ-1, "%.2f", elem->similarity);
      gt_feature_node_set_attribute((GtFeatureNode*) node,
                                    "ltr_similarity",
                                    buf);
      (void) snprintf(buf, BUFSIZ-1, ""GT_WU"", elem->contignumber);
      gt_feature_node_set_attribute((GtFeatureNode*) node,
                                    "seq_number",
                                    buf);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      parent = node;
      /* long_terminal_repeat */
      node = gt_feature_node_new(seqid,
                                 gt_ft_long_terminal_repeat,
                                 elem->leftLTR_5 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 elem->leftLTR_3 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      node = gt_feature_node_new(seqid,
                                 gt_ft_long_terminal_repeat,
                                 elem->rightLTR_5 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 elem->rightLTR_3 - seqstartpos + 1
                                   + ltrh_stream->offset,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      gt_feature_node_add_child((GtFeatureNode*) parent,
                                (GtFeatureNode*) node);
      gt_str_delete(seqid);
      gt_str_delete(source);
      ltrh_stream->cur_elem_index++;
    } else
    {
      *gn = NULL;
    }
  }
  return had_err;
}

static void gt_ltrharvest_stream_free(GtNodeStream *ns)
{
  GtLTRharvestStream *ltrh_stream = gt_ltrharvest_stream_cast(ns);
  GT_FREEARRAY(&ltrh_stream->arrayLTRboundaries, LTRboundaries);
  if (ltrh_stream->ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&ltrh_stream->ssar);
  if (ltrh_stream->bdptrtab != NULL)
    gt_free(ltrh_stream->bdptrtab);
}

const GtNodeStreamClass* gt_ltrharvest_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRharvestStream),
                                   gt_ltrharvest_stream_free,
                                   gt_ltrharvest_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

const GtEncseq* gt_ltrharvest_stream_get_encseq(GtNodeStream *ltrh_stream)
{
  GtLTRharvestStream *stream;
  gt_assert(ltrh_stream != NULL);
  stream = gt_ltrharvest_stream_cast(ltrh_stream);
  return stream->encseq;
}

void gt_ltrharvest_stream_disable_md5_seqids(GtLTRharvestStream *ltrh_stream)
{
  gt_assert(ltrh_stream != NULL);
  ltrh_stream->output_md5 = false;
}

void gt_ltrharvest_stream_enable_md5_seqids(GtLTRharvestStream *ltrh_stream)
{
  gt_assert(ltrh_stream != NULL);
  ltrh_stream->output_md5 = true;
}

void gt_ltrharvest_stream_disable_seqids(GtLTRharvestStream *ltrh_stream)
{
  gt_assert(ltrh_stream != NULL);
  ltrh_stream->output_seqids = false;
}

void gt_ltrharvest_stream_enable_seqids(GtLTRharvestStream *ltrh_stream)
{
  gt_assert(ltrh_stream != NULL);
  ltrh_stream->output_seqids = true;
}

GtNodeStream* gt_ltrharvest_stream_new(GtStr *str_indexname,
                                       GtRange searchrange,
                                       GtUword minseedlength,
                                       GtUword minltrlength,
                                       GtUword maxltrlength,
                                       GtUword mindistance,
                                       GtUword maxdistance,
                                       double similaritythreshold,
                                       int xdropbelowscore,
                                       GtXdropArbitraryscores arbitscores,
                                       GtLTRFourCharMotif *motif,
                                       bool verbosemode,
                                       bool nooverlaps,
                                       bool bestoverlaps,
                                       bool scanfile,
                                       GtUword offset,
                                       unsigned int minlengthTSD,
                                       unsigned int maxlengthTSD,
                                       GtUword vicinity,
                                       GtError *err)
{
  int had_err = 0;
  GtUword max_contiglength;
  GtNodeStream *ns = gt_node_stream_create(gt_ltrharvest_stream_class(), false);
  GtLTRharvestStream *ltrh_stream = gt_ltrharvest_stream_cast(ns);

  ltrh_stream->str_indexname = str_indexname;
  ltrh_stream->minseedlength = minseedlength;
  ltrh_stream->similaritythreshold = similaritythreshold;
  ltrh_stream->xdropbelowscore = xdropbelowscore;
  ltrh_stream->arbitscores = arbitscores;
  ltrh_stream->motif = motif;
  ltrh_stream->verbosemode = verbosemode;
  ltrh_stream->offset = offset;
  ltrh_stream->repeatinfo.ltrsearchseqrange = searchrange;
  ltrh_stream->minlengthTSD = minlengthTSD;
  ltrh_stream->maxlengthTSD = maxlengthTSD;
  ltrh_stream->nooverlaps = nooverlaps;
  ltrh_stream->bestoverlaps = bestoverlaps;
  ltrh_stream->vicinityforcorrectboundaries = vicinity;
  ltrh_stream->cur_elem_index = 0;
  ltrh_stream->output_md5 = false;
  ltrh_stream->output_seqids = false;
  ltrh_stream->state = GT_LTRHARVEST_STREAM_STATE_START;
  /* init array for maximal repeats */
  GT_INITARRAY(&ltrh_stream->arrayLTRboundaries, LTRboundaries);

  ltrh_stream->ssar =
    gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(str_indexname),
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                scanfile,
                                                NULL,
                                                err);
  if (ltrh_stream->ssar == NULL)
  {
    gt_node_stream_delete(ns);
    return NULL;
  }
  /* get encseq associated with suffix array */
  ltrh_stream->encseq = gt_encseqSequentialsuffixarrayreader(ltrh_stream->ssar);
  max_contiglength = gt_encseq_max_seq_length(ltrh_stream->encseq);
  ltrh_stream->repeatinfo.dmax = MIN(maxdistance,max_contiglength);
  ltrh_stream->repeatinfo.dmin = mindistance;
  ltrh_stream->repeatinfo.lmin = minltrlength;
  ltrh_stream->repeatinfo.lmax = maxltrlength;

  /* encode motif according to encseq alphabet */
  had_err = gt_ltr_four_char_motif_encode(motif, ltrh_stream->encseq, err);
  if (had_err)
  {
    gt_node_stream_delete(ns);
    return NULL;
  }

  return had_err ? NULL : ns;
}
