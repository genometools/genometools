#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/intbits.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "ft-front-prune.h"
#include "ft-trimstat.h"
#include "core/minmax.h"
#include "ft-polish.h"
#include "ft-front-generation.h"

#define UPDATE_MATCH_HISTORY(MC,MH)\
        if ((MH) & mask)\
        {\
          gt_assert((MC) > 0);\
          (MC)--;\
        }\
        (MH) <<= 1

typedef unsigned int Rowvaluetype;
typedef uint8_t Matchcounttype;
typedef uint8_t Backreferencetype;

typedef struct
{
  uint64_t matchhistory;
  Rowvaluetype row,
               localmatch_count;
  Matchcounttype matchhistory_count;
  Backreferencetype backreference;
} Frontvalue;

#ifndef OUTSIDE_OF_GT
typedef struct
{
  const GtTwobitencoding *twobitencoding;
  const GtEncseq *encseq;
  bool forward;
  GtEncseqReader *encseqreader;
  GtUchar *cache_ptr;
  GtAllocatedMemory *sequence_cache;
  GtUword substringlength,
          min_access_pos,
          cache_num_positions,
          cache_offset,
          startpos;
} Sequenceobject;

static void sequenceobject_init(Sequenceobject *seq,
                                GtExtendCharAccess extend_char_access_mode,
                                const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtUword startpos,
                                GtUword len,
                                GtEncseqReader *encseq_r,
                                GtAllocatedMemory *sequence_cache,
                                GtUword totallength
                                )
{
  gt_assert(seq != NULL);
  seq->encseq = NULL;
  seq->encseqreader = NULL;
  seq->twobitencoding = NULL;
  seq->cache_ptr = NULL;
  seq->sequence_cache = NULL;
  if (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY &&
      gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    seq->twobitencoding = gt_encseq_twobitencoding_export(encseq);
  }
  if (seq->twobitencoding == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER))
  {
    gt_encseq_reader_reinit_with_readmode(encseq_r, encseq, readmode, startpos);
    seq->encseqreader = encseq_r;
    gt_assert(seq->encseqreader != NULL);
    seq->sequence_cache = sequence_cache;
    gt_assert(sequence_cache != NULL);
    seq->cache_ptr = sequence_cache->space;
    seq->min_access_pos = GT_UWORD_MAX;
    seq->cache_num_positions = 0;
    seq->cache_offset = 0;
  }
  if (seq->twobitencoding == NULL && seq->encseqreader == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ))
  {
    seq->encseq = encseq;
  }
  seq->substringlength = len;
  if (readmode == GT_READMODE_FORWARD)
  {
    seq->startpos = startpos;
    seq->forward = true;
  } else
  {
    gt_assert(readmode == GT_READMODE_REVERSE);
    gt_assert(gt_encseq_total_length(encseq) == totallength);
    gt_assert(startpos + 1 <= totallength);
    seq->startpos = totallength - 1 - startpos;
    seq->forward = false;
  }
  gt_assert(seq->twobitencoding != NULL || seq->encseqreader != NULL ||
            seq->encseq != NULL);
}

static GtUchar gt_twobitencoding_char_at_pos(
                                      const GtTwobitencoding *twobitencoding,
                                      GtUword pos)
{
  return (twobitencoding[GT_DIVBYUNITSIN2BITENC(pos)] >>
          GT_MULT2(GT_UNITSIN2BITENC - 1 - GT_MODBYUNITSIN2BITENC(pos))) & 3;
}

static GtUchar sequenceobject_get_char(Sequenceobject *seq,GtUword pos)
{
  if (seq->twobitencoding != NULL)
  {
    return gt_twobitencoding_char_at_pos(seq->twobitencoding,
                                         seq->forward ? seq->startpos + pos
                                                      : seq->startpos - pos);
  }
  if (seq->encseqreader != NULL)
  {
    const GtUword addamount = 256UL;

    if (seq->min_access_pos != GT_UWORD_MAX &&
        seq->min_access_pos >= seq->cache_offset + addamount)
    {
      GtUword idx, end = MIN(seq->cache_num_positions,seq->substringlength);
      GtUchar *cs = ((GtUchar *) seq->sequence_cache->space)
                    - seq->min_access_pos;

      for (idx = seq->min_access_pos; idx < end; idx++)
      {
        cs[idx] = seq->cache_ptr[idx];
      }
      seq->cache_offset = seq->min_access_pos;
      seq->cache_ptr = ((GtUchar *) seq->sequence_cache->space)
                       - seq->cache_offset;
    }
    if (pos >= seq->cache_num_positions)
    {
      GtUword idx, tostore;

      tostore = MIN(seq->cache_num_positions + addamount,seq->substringlength);
      if (tostore > seq->cache_offset + seq->sequence_cache->allocated)
      {
        seq->sequence_cache->allocated += addamount;
        seq->sequence_cache->space
          = gt_realloc(seq->sequence_cache->space,
                       sizeof (GtUchar) * seq->sequence_cache->allocated);
        seq->cache_ptr = ((GtUchar *) seq->sequence_cache->space)
                         - seq->cache_offset;
      }
      gt_assert(pos >= seq->cache_offset);
      for (idx = seq->cache_num_positions; idx < tostore; idx++)
      {
        seq->cache_ptr[idx]
          = gt_encseq_reader_next_encoded_char(seq->encseqreader);
      }
      seq->cache_num_positions = tostore;
    }
    gt_assert(pos < seq->cache_offset + seq->sequence_cache->allocated);
    gt_assert(seq->cache_ptr != NULL);
    return seq->cache_ptr[pos];
  }
  gt_assert(seq->encseq != NULL);
  gt_assert(seq->forward || seq->startpos >= pos);
  return gt_encseq_get_encoded_char(seq->encseq,
                                    seq->forward ? seq->startpos + pos
                                                 : seq->startpos - pos,
                                    GT_READMODE_FORWARD);
}

#else
typedef struct
{
  const GtUchar *sequence_ptr;
  GtUword substringlength;
} Sequenceobject;

static void sequenceobject_init(Sequenceobject *seq,
                                const GtUchar *ptr,
                                GtUword startpos,
                                GtUword len)
{
  gt_assert(seq != NULL);
  seq->sequence_ptr = ptr + startpos;
  seq->substringlength = len;
}
#endif

#define FRONT_DIAGONAL(FRONTPTR) (GtWord) ((FRONTPTR) - midfront)

static bool sequenceobject_symbol_match(Sequenceobject *useq,
                                        GtUword upos,
                                        Sequenceobject *vseq,
                                        GtUword vpos)
{
#ifndef OUTSIDE_OF_GT
  GtUchar cu = sequenceobject_get_char(useq,upos);
  if (ISSPECIAL(cu))
  {
    return false;
  }
  return cu == sequenceobject_get_char(vseq,vpos) ? true : false;
#else
  GtUchar cu = useq->sequence_ptr[upos];
  return cu == vseq->sequence_ptr[vpos] ? true : false;
#endif
}

static void inline add_matches(Frontvalue *midfront,
                               Frontvalue *fv,
                               uint64_t mask,
                               Sequenceobject *useq,
                               Sequenceobject *vseq)
{
  GtUword upos, vpos;

  fv->localmatch_count = 0;
  for (upos = fv->row, vpos = fv->row + FRONT_DIAGONAL(fv);
       upos < useq->substringlength && vpos < vseq->substringlength &&
       sequenceobject_symbol_match(useq,upos,vseq,vpos);
       upos++, vpos++)
  {
    fv->localmatch_count++;
    if (!(fv->matchhistory & mask))
    {
      gt_assert(fv->matchhistory_count < INT8_MAX);
      fv->matchhistory_count++;
    }
    fv->matchhistory = (fv->matchhistory << 1) | (uint64_t) 1;
  }
  fv->row += fv->localmatch_count;
}

static GtUword front_next_inplace(Frontvalue *midfront,
                                  Frontvalue *lowfront,
                                  Frontvalue *highfront,
                                  GtUword history,
                                  Sequenceobject *useq,
                                  Sequenceobject *vseq)
{
  GtUword alignedlen, maxalignedlen;
  const uint64_t mask = ((uint64_t) 1) << (history-1);
  Frontvalue bestfront, insertion_value, replacement_value, *frontptr;

  insertion_value = *lowfront; /* from previous diag -(d-1) => -d => DELETION */
  bestfront = insertion_value;
  bestfront.row++;
  UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
  *lowfront = bestfront;
  lowfront->backreference = FT_EOP_DELETION;
  add_matches(midfront,lowfront,mask,useq,vseq);
  maxalignedlen = GT_MULT2(lowfront->row) + FRONT_DIAGONAL(lowfront);

  replacement_value = *(lowfront+1);
  if (bestfront.row < replacement_value.row + 1)
  {
    bestfront = replacement_value;
    bestfront.backreference = FT_EOP_DELETION;
    bestfront.row++;
    UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
  } else
  {
    bestfront.backreference = FT_EOP_REPLACEMENT;
    if (bestfront.row == replacement_value.row + 1)
    {
      bestfront.backreference |= FT_EOP_DELETION;
    }
  }
  *(lowfront+1) = bestfront;
  add_matches(midfront,lowfront + 1,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+1)->row) + FRONT_DIAGONAL(lowfront + 1);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }
  for (frontptr = lowfront+2; frontptr <= highfront; frontptr++)
  {
    bestfront = insertion_value;
    bestfront.backreference = FT_EOP_INSERTION;
    if (frontptr <= highfront - 1)
    {
      if (bestfront.row < replacement_value.row + 1)
      {
        bestfront = replacement_value;
        bestfront.backreference = FT_EOP_REPLACEMENT;
        bestfront.row++;
      } else
      {
        if (bestfront.row == replacement_value.row + 1)
        {
          bestfront.backreference |= FT_EOP_REPLACEMENT;
        }
      }
    }
    if (frontptr <= highfront - 2)
    {
      if (bestfront.row < frontptr->row + 1)
      {
        bestfront = *frontptr;
        bestfront.backreference = FT_EOP_DELETION;
        bestfront.row++;
      } else
      {
        if (bestfront.row == frontptr->row + 1)
        {
          bestfront.backreference |= FT_EOP_DELETION;
        }
      }
    }
    UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
    if (frontptr < highfront)
    {
      insertion_value = replacement_value;
      replacement_value = *frontptr;
    }
    *frontptr = bestfront;
    add_matches(midfront,frontptr,mask,useq,vseq);
    alignedlen = GT_MULT2(frontptr->row) + FRONT_DIAGONAL(frontptr);
    if (maxalignedlen < alignedlen)
    {
      maxalignedlen = alignedlen;
    }
  }
  return maxalignedlen;
}

static GtUword front_second_inplace(Frontvalue *midfront,
                                    Frontvalue *lowfront,
                                    GtUword history,
                                    Sequenceobject *useq,
                                    Sequenceobject *vseq)
{
  GtUword alignedlen, maxalignedlen;
  const uint64_t mask = ((uint64_t) 1) << (history-1);

  *(lowfront+1) = *(lowfront+2) = *lowfront;
  lowfront->row++;
  lowfront->backreference = FT_EOP_DELETION;
  UPDATE_MATCH_HISTORY(lowfront->matchhistory_count,lowfront->matchhistory);
  add_matches(midfront,lowfront,mask,useq,vseq);
  maxalignedlen = GT_MULT2(lowfront->row) + FRONT_DIAGONAL(lowfront);

  (lowfront+1)->row++;
  (lowfront+1)->backreference = FT_EOP_REPLACEMENT;
  UPDATE_MATCH_HISTORY((lowfront+1)->matchhistory_count,
                       (lowfront+1)->matchhistory);
  add_matches(midfront,lowfront + 1,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+1)->row) + FRONT_DIAGONAL(lowfront + 1);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }

  (lowfront+2)->backreference = FT_EOP_INSERTION;
  UPDATE_MATCH_HISTORY((lowfront+2)->matchhistory_count,
                       (lowfront+2)->matchhistory);
  add_matches(midfront,lowfront + 2,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+2)->row) + FRONT_DIAGONAL(lowfront + 2);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }
  return maxalignedlen;
}

#undef TRIM_INFO_OUT
#ifdef TRIM_INFO_OUT
static bool trimthisentry(GtUword distance,
                          Rowvaluetype row,
                          GtWord diagonal,
                          GtUword minlenforhistorycheck,
                          Matchcounttype matchhistory_count,
                          GtUword minmatchnum,
                          GtUword minlenfrommaxdiff)
{
  GtUword alignedlen = GT_MULT2(row) + diagonal;

  if (alignedlen >= minlenforhistorycheck && matchhistory_count < minmatchnum)
  {
    printf(GT_WD "&" GT_WU "&%u&1: matches=%d < " GT_WU "=minmatches\n",
              diagonal,distance,row,(int) matchhistory_count,minmatchnum);
    return true;
  }
  if (alignedlen < minlenfrommaxdiff)
  {
    printf(GT_WD "&" GT_WU "&%u&2: i'+j'=" GT_WU "<" GT_WU "=i+j-lag\n",
              diagonal,distance,row,alignedlen,minlenfrommaxdiff);
    return true;
  }
  printf(GT_WD "&" GT_WU "&%u\n", diagonal,distance,row);
  return false;
}
#else
static bool trimthisentry(Rowvaluetype row,
                          GtWord diagonal,
                          GtUword minlenforhistorycheck,
                          Matchcounttype matchhistory_count,
                          GtUword minmatchnum,
                          GtUword minlenfrommaxdiff)
{
  GtUword alignedlen = GT_MULT2(row) + diagonal;

  if (alignedlen >= minlenforhistorycheck && matchhistory_count < minmatchnum)
  {
    return true;
  }
  if (alignedlen < minlenfrommaxdiff)
  {
    return true;
  }
  return false;
}
#endif

static GtUword trim_front(bool upward,
#ifdef TRIM_INFO_OUT
                          GtUword distance,
#endif
                          GtUword ulen,
                          GtUword vlen,
                          GtUword minmatchnum,
                          GtUword minlenforhistorycheck,
                          GtUword minlenfrommaxdiff,
                          const Frontvalue *midfront,
                          const Frontvalue *from,
                          const Frontvalue *stop)
{
  const Frontvalue *frontptr;
  GtUword trim = 0;

  gt_assert ((upward && from < stop) || (!upward && stop < from));
  for (frontptr = from; frontptr != stop; frontptr = upward ? (frontptr + 1)
                                                            : (frontptr - 1))
  {
    if (frontptr->row > ulen ||
        frontptr->row + FRONT_DIAGONAL(frontptr) > vlen ||
        trimthisentry(
#ifdef TRIM_INFO_OUT
                        distance,
#endif
                        frontptr->row,
                        FRONT_DIAGONAL(frontptr),
                        minlenforhistorycheck,
                        frontptr->matchhistory_count,
                        minmatchnum,
                        minlenfrommaxdiff))
    {
      trim++;
    } else
    {
      break;
    }
  }
  return trim;
}

static void frontspace_check(GT_UNUSED const Frontvalue *from,
                             GT_UNUSED const Frontvalue *to,
                             GT_UNUSED const Frontvalue *ptr)
{
  gt_assert (ptr >= from && ptr <= to);
}

static Frontvalue *frontspace_allocate(GtUword minsizeforshift,
                                       GtUword trimleft,
                                       GtUword valid,
                                       GtAllocatedMemory *fs)
{
  if (trimleft - fs->offset + valid >= fs->allocated)
  {
    fs->allocated = 255UL + MAX(fs->allocated * 1.2,
                                           trimleft - fs->offset + valid);
    gt_assert(fs->allocated > trimleft - fs->offset + valid);
    fs->space = gt_realloc(fs->space,sizeof (Frontvalue) * fs->allocated);
    gt_assert(fs->space != NULL);
  }
  gt_assert(trimleft >= fs->offset);
  if (trimleft - fs->offset > MAX(valid,minsizeforshift))
  {
    memcpy(fs->space,((Frontvalue *) fs->space) + trimleft - fs->offset,
           sizeof (Frontvalue) * valid);
    fs->offset = trimleft;
  }
  return ((Frontvalue *) fs->space) - fs->offset;
}

static void update_trace_and_polished(Polished_point *best_polished_point,
#ifndef OUTSIDE_OF_GT
                                      GtUword *minrow,
                                      GtUword *mincol,
#endif
                                      Fronttrace *front_trace,
                                      const Polishing_info *pol_info,
                                      GtUword distance,
                                      GtUword trimleft,
                                      Frontvalue *midfront,
                                      Frontvalue *lowfront,
                                      Frontvalue *highfront)
{
  const Frontvalue *frontptr;
  uint64_t lsb;

#ifndef OUTSIDE_OF_GT
  *minrow = GT_UWORD_MAX;
  *mincol = GT_UWORD_MAX;
#endif
  for (frontptr = lowfront; frontptr <= highfront; frontptr++)
  {
    GtUword alignedlen = GT_MULT2(frontptr->row) + FRONT_DIAGONAL(frontptr);

#ifndef OUTSIDE_OF_GT
    GtUword currentcol;

    if (*minrow > frontptr->row)
    {
      *minrow = frontptr->row;
    }
    gt_assert(FRONT_DIAGONAL(frontptr) >= 0 ||
              frontptr->row >= -FRONT_DIAGONAL(frontptr));
    currentcol = frontptr->row + FRONT_DIAGONAL(frontptr);
    if (*mincol > currentcol)
    {
      *mincol = currentcol;
    }
#endif
    lsb = frontptr->matchhistory & pol_info->mask;
    if (HISTORY_IS_POLISHED(pol_info,frontptr->matchhistory,lsb) &&
        alignedlen > best_polished_point->alignedlen)
    {
      best_polished_point->alignedlen = alignedlen;
      best_polished_point->row = frontptr->row;
      best_polished_point->distance = distance;
      best_polished_point->trimleft = trimleft;
    }
    if (front_trace != NULL)
    {
      front_trace_add_trace(front_trace,frontptr->backreference,
                            frontptr->localmatch_count);
    }
  }
}

GtUword front_prune_edist_inplace(
#ifndef OUTSIDE_OF_GT
                         bool forward,
                         GtAllocatedMemory *frontspace,
#endif
                         Trimstat *trimstat,
                         Polished_point *best_polished_point,
                         Fronttrace *front_trace,
                         const Polishing_info *pol_info,
                         GtUword history,
                         GtUword minmatchnum,
                         GtUword maxalignedlendifference,
                         FTsequenceResources *ufsr,
                         GtUword ustart,
                         GtUword ulen,
                         FTsequenceResources *vfsr,
                         GtUword vstart,
                         GtUword vlen)
{
  const GtUword sumseqlength = ulen + vlen,
                minsizeforshift = sumseqlength/1000,
                minlenforhistorycheck = GT_MULT2(history);
  /* so the space for allocating the fronts is
     sizeof (Frontvalue) * ((m+n)/1000 + maxvalid), where maxvalid is a small
     constant. */
  GtUword distance, trimleft = 0, valid = 1UL, maxvalid = 0, sumvalid = 0;
  const uint64_t mask = ((uint64_t) 1) << (history-1);
  Frontvalue *validbasefront;
  bool diedout = false;
  Sequenceobject useq, vseq;

#ifdef OUTSIDE_OF_GT
  GtAllocatedMemory *frontspace = gt_malloc(sizeof *frontspace);
  frontspace->space = NULL;
  frontspace->allocated = 0;
  frontspace->offset = 0;
  sequenceobject_init(&useq,useqptr,ustart,ulen);
  sequenceobject_init(&vseq,vseqptr,vstart,vlen);
#else
  GtReadmode readmode = forward ? GT_READMODE_FORWARD : GT_READMODE_REVERSE;
  sequenceobject_init(&useq,ufsr->extend_char_access,ufsr->encseq,readmode,
                      ustart,ulen,ufsr->encseq_r,ufsr->sequence_cache,
                      ufsr->totallength);
  sequenceobject_init(&vseq,ufsr->extend_char_access,vfsr->encseq,readmode,
                      vstart,vlen,vfsr->encseq_r,vfsr->sequence_cache,
                      vfsr->totallength);
  frontspace->offset = 0;
#endif
#ifdef TRIM_INFO_OUT
  printf("regionalquality(minmatchnum)=" GT_WU "\n",minmatchnum);
#endif
  for (distance = 0, valid = 1UL; /* Nothing */; distance++, valid += 2)
  {
    GtUword trim, maxalignedlen, minlenfrommaxdiff;

#ifdef TRIM_INFO_OUT
    printf("distance=" GT_WU ",full=" GT_WU ",trimleft=" GT_WU
           ",valid=" GT_WU "\n",distance,
                  GT_MULT2(distance) + 1,
                  trimleft,valid);
#endif
    gt_assert(valid <= GT_MULT2(distance) + 1);
    sumvalid += valid;
    if (maxvalid < valid)
    {
      maxvalid = valid;
    }
    validbasefront = frontspace_allocate(minsizeforshift,trimleft,valid,
                                         frontspace);
    if (distance == 0)
    {
      validbasefront->row = 0;
      validbasefront->matchhistory = 0;
      validbasefront->matchhistory_count = 0;
      validbasefront->backreference = 0; /* No back reference */
      add_matches(validbasefront + distance,validbasefront,mask,&useq,&vseq);
      maxalignedlen = GT_MULT2(validbasefront->row);
    } else
    {
      gt_assert(valid >= 3UL);
      frontspace_check((const Frontvalue *) frontspace->space,
                       ((const Frontvalue *) frontspace->space)
                        + frontspace->allocated - 1,
                       validbasefront + trimleft);
      frontspace_check((const Frontvalue *) frontspace->space,
                       ((const Frontvalue *) frontspace->space)
                         + frontspace->allocated - 1,
                       validbasefront + trimleft + valid - 1);
      if (valid == 3UL)
      {
        maxalignedlen
          = front_second_inplace(validbasefront + distance,
                                 validbasefront + trimleft,
                                 history,
                                 &useq,
                                 &vseq);
      } else
      {
        maxalignedlen
          = front_next_inplace(validbasefront + distance,
                               validbasefront + trimleft,
                               validbasefront + trimleft + valid - 1,
                               history,
                               &useq,
                               &vseq);
      }
    }
    gt_assert(valid > 0);
    minlenfrommaxdiff = maxalignedlen >= maxalignedlendifference
                          ? maxalignedlen - maxalignedlendifference
                          : 0;
#ifdef TRIM_INFO_OUT
    printf("maxalignedlen=" GT_WU ",maxlenfrommaxdiff=" GT_WU "\n",
           maxalignedlen,minlenfrommaxdiff);
#endif
    trim = trim_front(true,
#ifdef TRIM_INFO_OUT
                      distance,
#endif
                      ulen,
                      vlen,
                      minmatchnum,
                      minlenforhistorycheck,
                      minlenfrommaxdiff,
                      validbasefront + distance,
                      validbasefront + trimleft,
                      validbasefront + trimleft + valid);
#ifdef TRIM_INFO_OUT
    printf("trim on left=" GT_WU "\n",trim);
#endif
    if (trim > 0)
    {
      trimleft += trim;
      gt_assert(valid >= trim);
      valid -= trim;
    }
    if (valid > 0)
    {
      trim = trim_front(false,
#ifdef TRIM_INFO_OUT
                        distance,
#endif
                        ulen,
                        vlen,
                        minmatchnum,
                        minlenforhistorycheck,
                        minlenfrommaxdiff,
                        validbasefront + distance,
                        validbasefront + trimleft + valid - 1,
                        validbasefront + trimleft - 1);
#ifdef TRIM_INFO_OUT
      printf("trim on right=" GT_WU "\n",trim);
#endif
      gt_assert(trim < valid);
      if (trim > 0)
      {
        gt_assert(valid >= trim);
        valid -= trim;
      }
    }
    if (valid == 0)
    {
      diedout = true;
      break;
    }
    if (front_trace != NULL)
    {
      front_trace_add_gen(front_trace,trimleft,valid);
    }
    update_trace_and_polished(best_polished_point,
#ifndef OUTSIDE_OF_GT
                              &useq.min_access_pos,
                              &vseq.min_access_pos,
#endif
                              front_trace,
                              pol_info,
                              distance,
                              trimleft,
                              validbasefront + distance,
                              validbasefront + trimleft,
                              validbasefront + trimleft + valid - 1);
    if ((vlen > ulen && vlen - ulen <= distance) ||
        (vlen <= ulen && ulen - vlen <= distance))
    {
      if (distance + vlen - ulen >= trimleft &&
          distance + vlen - ulen <= trimleft + valid - 1 &&
          validbasefront[distance + vlen - ulen].row == ulen)
      {
        break;
      }
    }
    if (distance >= sumseqlength)
    {
      break;
    }
  }
  trimstat_add(trimstat,diedout,sumvalid,maxvalid,distance,
               sizeof (Frontvalue) * frontspace->allocated,
#ifndef OUTSIDE_OF_GT
               useq.sequence_cache != NULL &&
               vseq.sequence_cache != NULL ? MAX(useq.sequence_cache->allocated,
                                                 vseq.sequence_cache->allocated)
                                           : 0
#else
               0
#endif
              );
  return diedout ? sumseqlength + 1 : distance;
}
