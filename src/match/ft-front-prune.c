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
#include "core/minmax.h"
#include "core/encseq.h"
#include "core/bitpackstringsimpleop.h"
#include "match/extend-offset.h"
#include "match/ft-front-prune.h"
#include "match/ft-trimstat.h"
#include "match/ft-polish.h"
#include "match/ft-front-generation.h"

#define GT_UPDATE_MATCH_HISTORY(FRONTVAL)\
        if ((FRONTVAL)->matchhistory_size < max_history)\
        {\
          gt_assert((FRONTVAL)->matchhistory_size < max_history);\
          (FRONTVAL)->matchhistory_size++;\
        }\
        (FRONTVAL)->matchhistory_bits <<= 1

typedef uint32_t GtFtRowvaluetype;
typedef uint8_t GtFtMatchcounttype;
typedef uint8_t GtFtBackreferencetype;

typedef struct
{
  uint64_t matchhistory_bits;
  GtFtRowvaluetype row,
                   localmatch_count;
  GtFtMatchcounttype matchhistory_size;
  GtFtBackreferencetype backreference;
  uint32_t max_mismatches; /* maximum number of mismatches in a path to this
                              Front-entry.*/
} GtFtFrontvalue;

typedef struct
{
  const GtTwobitencoding *twobitencoding;
  const GtEncseq *encseq;
  const GtUchar *bytesequenceptr;
  GtEncseqReader *encseqreader;
  GtUchar *cache_ptr;
  GtAllocatedMemory *sequence_cache;
  GtUword substringlength,
          totallength,
          cache_num_positions; /* number of positions in cache */
  GtUword offset,
          seqstartpos;
  bool read_seq_left2right,
       dir_is_complement;
} GtFtSequenceObject;

static void ft_sequenceobject_init(GtFtSequenceObject *seq,
                                   GtExtendCharAccess extend_char_access_mode,
                                   const GtEncseq *encseq,
                                   bool rightextension,
                                   GtReadmode readmode,
                                   GtUword seqstartpos,
                                   GtUword startpos,
                                   GtUword len,
                                   GtEncseqReader *encseq_r,
                                   GtAllocatedMemory *sequence_cache,
                                   const GtUchar *bytesequence,
                                   GtUword totallength
                                   )
{
  gt_assert(seq != NULL);
  seq->encseq = NULL;
  seq->encseqreader = NULL;
  seq->twobitencoding = NULL;
  seq->cache_ptr = NULL;
  seq->sequence_cache = NULL;
  seq->bytesequenceptr = NULL;
  seq->seqstartpos = seqstartpos;
  gt_assert(seqstartpos <= startpos);
  seq->offset = GT_EXTEND_OFFSET(rightextension,
                                 readmode,
                                 totallength,
                                 seqstartpos,
                                 startpos,
                                 len);
  seq->read_seq_left2right = GT_EXTEND_READ_SEQ_LEFT2RIGHT(rightextension,
                                                           readmode);
  if (encseq != NULL && extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY &&
      gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    seq->twobitencoding = gt_encseq_twobitencoding_export(encseq);
  }
  if (encseq != NULL && seq->twobitencoding == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER))
  {
    GtUword full_totallength = gt_encseq_total_length(encseq);
    gt_encseq_reader_reinit_with_readmode(encseq_r, encseq,
                                          seq->read_seq_left2right
                                            ? GT_READMODE_FORWARD
                                            : GT_READMODE_REVERSE,
                                          seq->read_seq_left2right
                                            ? seq->offset
                                            : GT_REVERSEPOS(full_totallength,
                                                            seq->offset));
    seq->encseqreader = encseq_r;
    gt_assert(seq->encseqreader != NULL && sequence_cache != NULL);
    seq->sequence_cache = sequence_cache;
    seq->cache_ptr = (GtUchar *) sequence_cache->space;
    seq->cache_num_positions = 0;
  }
  if (encseq != NULL && seq->twobitencoding == NULL &&
      seq->encseqreader == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ))
  {
    seq->encseq = encseq;
  }
  if (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_DIRECT)
  {
    gt_assert(seq->twobitencoding == NULL && seq->encseqreader == NULL &&
              seq->encseq == NULL);
    seq->bytesequenceptr = bytesequence;
  }
  seq->substringlength = len;
  seq->totallength = totallength;
  seq->dir_is_complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;
  gt_assert(seq->twobitencoding != NULL || seq->encseqreader != NULL ||
            seq->encseq != NULL || seq->bytesequenceptr != NULL);
}

static GtUchar gt_twobitencoding_char_at_pos(
                                      const GtTwobitencoding *twobitencoding,
                                      GtUword pos)
{
  return (twobitencoding[GT_DIVBYUNITSIN2BITENC(pos)] >>
          GT_MULT2(GT_UNITSIN2BITENC - 1 - GT_MODBYUNITSIN2BITENC(pos))) & 3;
}

static GtUchar ft_sequenceobject_get_char(GtFtSequenceObject *seq,GtUword idx)
{
  GtUchar cc;
  GtUword accesspos;

  if (seq->twobitencoding != NULL)
  {
    gt_assert (seq->read_seq_left2right || seq->offset >= idx);
    accesspos = seq->read_seq_left2right ? seq->offset + idx
                                         : seq->offset - idx;
    gt_assert(accesspos < seq->seqstartpos + seq->totallength);
    cc = gt_twobitencoding_char_at_pos(seq->twobitencoding, accesspos);
    return seq->dir_is_complement ? GT_COMPLEMENTBASE(cc) : cc;
  }
  if (seq->encseqreader != NULL)
  {
    gt_assert(idx < seq->substringlength);
    if (idx >= seq->cache_num_positions)
    {
      GtUword cidx, tostore;
      const GtUword addamount = 256UL;

      tostore = MIN(seq->cache_num_positions + addamount,
                    seq->substringlength);
      if (tostore > seq->sequence_cache->allocated)
      {
        seq->sequence_cache->allocated += addamount;
        seq->sequence_cache->space
          = gt_realloc(seq->sequence_cache->space,
                       sizeof (GtUchar) * seq->sequence_cache->allocated);
        seq->cache_ptr = (GtUchar *) seq->sequence_cache->space;
      }
      for (cidx = seq->cache_num_positions; cidx < tostore; cidx++)
      {
        seq->cache_ptr[cidx]
          = gt_encseq_reader_next_encoded_char(seq->encseqreader);
      }
      seq->cache_num_positions = tostore;
    }
    gt_assert(seq->cache_ptr != NULL && idx < seq->cache_num_positions);
    cc = seq->cache_ptr[idx];
  } else
  {
    gt_assert (seq->read_seq_left2right || seq->offset >= idx);
    accesspos = seq->read_seq_left2right ? seq->offset + idx
                                         : seq->offset - idx;
    if (seq->encseq != NULL)
    {
      cc = gt_encseq_get_encoded_char(seq->encseq,accesspos,
                                      GT_READMODE_FORWARD);
    } else
    {
      gt_assert(seq->bytesequenceptr != NULL);
      cc = seq->bytesequenceptr[accesspos];
    }
  }
  if (seq->dir_is_complement && !ISSPECIAL(cc))
  {
    return GT_COMPLEMENTBASE(cc);
  }
  return cc;
}

static bool ft_sequenceobject_symbol_match(GtFtSequenceObject *useq,
                                           GtUword upos,
                                           GtFtSequenceObject *vseq,
                                           GtUword vpos)
{
  GtUchar cu = ft_sequenceobject_get_char(useq,upos);
  return (!ISSPECIAL(cu) && cu == ft_sequenceobject_get_char(vseq,vpos))
           ? true
           : false;
}

#define GT_FRONT_DIAGONAL(FRONTPTR) (GtWord) ((FRONTPTR) - midfront)

static inline GtUword front_prune_longest_common(const GtFtFrontvalue *midfront,
                                                 const GtFtFrontvalue *fv,
                                                 GtFtSequenceObject *useq,
                                                 GtFtSequenceObject *vseq)
{
  GtUword upos, vpos;

  for (upos = fv->row, vpos = fv->row + GT_FRONT_DIAGONAL(fv);
       upos < useq->substringlength && vpos < vseq->substringlength &&
       ft_sequenceobject_symbol_match(useq,upos,vseq,vpos);
       upos++, vpos++)
    /* Nothing */ ;
  return upos - fv->row;
}

static void inline front_prune_add_matches(GtFtFrontvalue *midfront,
                                           GtFtFrontvalue *fv,
                                           uint64_t max_history_mask,
                                           GtUword max_history,
                                           GtFtSequenceObject *useq,
                                           GtFtSequenceObject *vseq,
                                           GtFtTrimstat *trimstat)
{
  uint64_t match_mask;

  fv->localmatch_count
    = (GtFtRowvaluetype) front_prune_longest_common(midfront,fv,useq,vseq);
  fv->row += fv->localmatch_count;
  match_mask = (fv->localmatch_count >= max_history)
                 ? max_history_mask
                 : ((((uint64_t) 1) << fv->localmatch_count) - 1);
  fv->matchhistory_bits = (fv->matchhistory_bits << fv->localmatch_count) |
                          match_mask;
  fv->matchhistory_size = MIN(fv->matchhistory_size + fv->localmatch_count,
                              max_history);
  if (trimstat != NULL)
  {
    gt_ft_trimstat_add_matchlength(trimstat,fv->localmatch_count);
  }
}

static GtUword front_next_inplace(GtFtFrontvalue *midfront,
                                  GtFtFrontvalue *lowfront,
                                  GtFtFrontvalue *highfront,
                                  GtUword max_history,
                                  uint64_t max_history_mask,
                                  GtFtSequenceObject *useq,
                                  GtFtSequenceObject *vseq,
                                  GtFtTrimstat *trimstat)
{
  GtUword alignedlen, maxalignedlen;
  GtFtFrontvalue bestfront, insertion_value, replacement_value, *frontptr;

  insertion_value = *lowfront; /* from previous diag -(d-1) => -d => DELETION */
  bestfront = insertion_value;
  bestfront.row++;
  GT_UPDATE_MATCH_HISTORY(&bestfront);
  *lowfront = bestfront;
  lowfront->backreference = FT_EOP_DELETION;
  front_prune_add_matches(midfront,lowfront,max_history_mask,
                          max_history,useq,vseq,
                          trimstat);
  maxalignedlen = GT_MULT2(lowfront->row) + GT_FRONT_DIAGONAL(lowfront);

  replacement_value = *(lowfront+1);
  if (bestfront.row < replacement_value.row + 1)
  {
    bestfront = replacement_value;
    bestfront.backreference = FT_EOP_DELETION;
    bestfront.row++;
    GT_UPDATE_MATCH_HISTORY(&bestfront);
  } else
  {
    bestfront.backreference = FT_EOP_MISMATCH;
    bestfront.max_mismatches++;
    if (bestfront.row == replacement_value.row + 1)
    {
      bestfront.backreference |= FT_EOP_DELETION;
      if (bestfront.max_mismatches < replacement_value.max_mismatches)
      {
        bestfront.max_mismatches = replacement_value.max_mismatches;
      }
    }
  }
  *(lowfront+1) = bestfront;
  front_prune_add_matches(midfront,lowfront + 1,max_history_mask,
                          max_history,useq,vseq,
                          trimstat);
  alignedlen = GT_MULT2((lowfront+1)->row) + GT_FRONT_DIAGONAL(lowfront + 1);
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
        bestfront.backreference = FT_EOP_MISMATCH;
        bestfront.max_mismatches++;
        bestfront.row++;
      } else
      {
        if (bestfront.row == replacement_value.row + 1)
        {
          bestfront.backreference |= FT_EOP_MISMATCH;
          if (bestfront.max_mismatches < replacement_value.max_mismatches + 1)
          {
            bestfront.max_mismatches = replacement_value.max_mismatches + 1;
          }
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
    GT_UPDATE_MATCH_HISTORY(&bestfront);
    if (frontptr < highfront)
    {
      insertion_value = replacement_value;
      replacement_value = *frontptr;
    }
    *frontptr = bestfront;
    front_prune_add_matches(midfront,frontptr,max_history_mask,
                            max_history,useq,vseq,
                            trimstat);
    alignedlen = GT_MULT2(frontptr->row) + GT_FRONT_DIAGONAL(frontptr);
    if (maxalignedlen < alignedlen)
    {
      maxalignedlen = alignedlen;
    }
  }
  return maxalignedlen;
}

static GtUword front_second_inplace(GtFtFrontvalue *midfront,
                                    GtFtFrontvalue *lowfront,
                                    GtUword max_history,
                                    uint64_t max_history_mask,
                                    GtFtSequenceObject *useq,
                                    GtFtSequenceObject *vseq,
                                    GtFtTrimstat *trimstat)
{
  GtUword alignedlen, maxalignedlen;

  *(lowfront+1) = *(lowfront+2) = *lowfront;
  lowfront->row++;
  lowfront->backreference = FT_EOP_DELETION;
  GT_UPDATE_MATCH_HISTORY(lowfront);
  front_prune_add_matches(midfront,lowfront,max_history_mask,
                          max_history,useq,vseq,
                          trimstat);
  maxalignedlen = GT_MULT2(lowfront->row) + GT_FRONT_DIAGONAL(lowfront);

  (lowfront+1)->row++;
  (lowfront+1)->backreference = FT_EOP_MISMATCH;
  (lowfront+1)->max_mismatches++;
  GT_UPDATE_MATCH_HISTORY(lowfront+1);
  front_prune_add_matches(midfront,lowfront + 1,max_history_mask,
                          max_history,useq,vseq,
                          trimstat);
  alignedlen = GT_MULT2((lowfront+1)->row) + GT_FRONT_DIAGONAL(lowfront + 1);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }

  (lowfront+2)->backreference = FT_EOP_INSERTION;
  GT_UPDATE_MATCH_HISTORY(lowfront+2);
  front_prune_add_matches(midfront,lowfront + 2,max_history_mask,
                          max_history,useq,vseq,
                          trimstat);
  alignedlen = GT_MULT2((lowfront+2)->row) + GT_FRONT_DIAGONAL(lowfront + 2);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }
  return maxalignedlen;
}

#if defined (__GNUC__) && defined (__POPCNT__)
static inline unsigned
bitCountUInt64(uint64_t v)
{
  return __builtin_popcountl(v);
}
#else
static inline unsigned
bitCountUInt64(uint64_t v)
{
  return bitCountUInt32((uint32_t) (v & UINT32_MAX)) +
         bitCountUInt32((uint32_t) (v >> 32));
}
#endif

static bool trimthisentry(GtFtRowvaluetype row,
                          GtWord diagonal,
                          const GtFtFrontvalue *fv,
                          GtUword minmatchpercentage128,
                          GtUword minlenfrommaxdiff,
                          uint64_t max_history_mask,
                          GT_UNUSED GtUword distance,
                          GT_UNUSED bool showfrontinfo)
{
  if (GT_MULT2(row) + diagonal < minlenfrommaxdiff)
  {
    return true;
  }
  if (bitCountUInt64(fv->matchhistory_bits & max_history_mask)
      < ((fv->matchhistory_size * minmatchpercentage128) >> 7))
  {
    return true;
  }
  return false;
}

static GtUword trim_front(bool upward,
                          GtUword distance,
                          GtUword ulen,
                          GtUword vlen,
                          GtUword minmatchpercentage128,
                          GtUword minlenfrommaxdiff,
                          GtTrimmingStrategy trimstrategy,
                          const GtFtPolished_point *best_polished_point,
                          const GtFtFrontvalue *midfront,
                          const GtFtFrontvalue *from,
                          const GtFtFrontvalue *stop,
                          uint64_t max_history_mask,
                          bool showfrontinfo)
{
  const GtFtFrontvalue *frontptr;
  const int step = upward ? +1 : -1;

  if (trimstrategy == GT_OUTSENSE_TRIM_NEVER ||
      (trimstrategy == GT_OUTSENSE_TRIM_ON_NEW_PP &&
       best_polished_point != NULL &&
       best_polished_point->distance + 1 < distance &&
       best_polished_point->distance + 30 >= distance))
  {
    return 0;
  }
  gt_assert ((upward && from < stop) || (!upward && stop < from));
  for (frontptr = from; frontptr != stop; frontptr += step)
  {
    if (frontptr->row <= ulen &&
        frontptr->row + GT_FRONT_DIAGONAL(frontptr) <= vlen &&
        !trimthisentry(frontptr->row,
                       GT_FRONT_DIAGONAL(frontptr),
                       frontptr,
                       minmatchpercentage128,
                       minlenfrommaxdiff,
                       max_history_mask,
                       distance,
                       showfrontinfo))
    {
      break;
    }
  }
  if (upward)
  {
    return (GtUword) (frontptr - from);
  }
  return (GtUword) (from - frontptr);
}

static void frontspace_check(GT_UNUSED const GtFtFrontvalue *from,
                             GT_UNUSED const GtFtFrontvalue *to,
                             GT_UNUSED const GtFtFrontvalue *ptr)
{
  gt_assert (ptr >= from && ptr <= to);
}

static GtFtFrontvalue *frontspace_allocate(GtUword minsizeforshift,
                                       GtUword trimleft,
                                       GtUword valid,
                                       GtAllocatedMemory *fs)
{
  if (trimleft - fs->offset + valid >= fs->allocated)
  {
    fs->allocated = 255UL + MAX(fs->allocated * 1.2,
                                           trimleft - fs->offset + valid);
    gt_assert(fs->allocated > trimleft - fs->offset + valid);
    fs->space = gt_realloc(fs->space,sizeof (GtFtFrontvalue) * fs->allocated);
    gt_assert(fs->space != NULL);
  }
  gt_assert(trimleft >= fs->offset);
  if (trimleft - fs->offset > MAX(valid,minsizeforshift))
  {
    memcpy(fs->space,((GtFtFrontvalue *) fs->space) + trimleft - fs->offset,
           sizeof (GtFtFrontvalue) * valid);
    fs->offset = trimleft;
  }
  return ((GtFtFrontvalue *) fs->space) - fs->offset;
}

static void ft_update_trace_and_polished(
                                      GtFtPolished_point *best_polished_point,
                                      GtFronttrace *front_trace,
                                      const GtFtPolishing_info *pol_info,
                                      GtUword distance,
                                      GtUword trimleft,
                                      GtFtFrontvalue *midfront,
                                      GtFtFrontvalue *lowfront,
                                      GtFtFrontvalue *highfront,
                                      GT_UNUSED bool showfrontinfo)
{
  const GtFtFrontvalue *frontptr;

  for (frontptr = lowfront; frontptr <= highfront; frontptr++)
  {
    GtUword alignedlen = GT_MULT2(frontptr->row) + GT_FRONT_DIAGONAL(frontptr);

    gt_assert(GT_FRONT_DIAGONAL(frontptr) >= 0 ||
              frontptr->row >= -GT_FRONT_DIAGONAL(frontptr));
    if (alignedlen > best_polished_point->alignedlen)
    {
      uint64_t filled_matchhistory_bits = frontptr->matchhistory_bits;

      if (frontptr->matchhistory_size < pol_info->pol_size)
      {
        const int shift = pol_info->pol_size - frontptr->matchhistory_size;
        const uint64_t fill_bits = ((uint64_t) 1 << shift) - 1;
        filled_matchhistory_bits |= (fill_bits << frontptr->matchhistory_size);
      }
      if  (GT_HISTORY_IS_POLISHED(pol_info,filled_matchhistory_bits))
      {
        best_polished_point->alignedlen = alignedlen;
        best_polished_point->row = frontptr->row;
        best_polished_point->distance = distance;
        best_polished_point->trimleft = trimleft;
        best_polished_point->max_mismatches
          = (GtUword) frontptr->max_mismatches;
      }
    }
    if (front_trace != NULL)
    {
      front_trace_add_trace(front_trace,frontptr->backreference,
                            frontptr->localmatch_count);
    }
  }
}

GtUword front_prune_edist_inplace(
                         bool rightextension,
                         GtAllocatedMemory *frontspace,
                         GtFtPolished_point *best_polished_point,
                         GtFronttrace *front_trace,
                         const GtFtPolishing_info *pol_info,
                         GtTrimmingStrategy trimstrategy,
                         GtUword max_history,
                         GtUword minmatchpercentage,
                         GtUword maxalignedlendifference,
                         bool showfrontinfo,
                         GtUword seedlength,
                         GtFTsequenceResources *ufsr,
                         GtUword ustart,
                         GtUword ulen,
                         GtUword vseqstartpos,
                         GtFTsequenceResources *vfsr,
                         GtUword vstart,
                         GtUword vlen,
                         GtFtTrimstat *trimstat)
{
  const GtUword sumseqlength = ulen + vlen,
                minsizeforshift = sumseqlength/1000;
  /* so the space for allocating the fronts is
     sizeof (GtFtFrontvalue) * ((m+n)/1000 + maxvalid), where maxvalid is a
     small constant. */
  GtUword distance, trimleft = 0, valid = 1UL, maxvalid = 0, sumvalid = 0;
  GtFtFrontvalue *validbasefront;
  bool diedout = false;
  GtFtSequenceObject useq, vseq;
  /* This transformation allows to compute the /100  by >> 7 */
  const GtUword minmatchpercentage128
    = (minmatchpercentage * 128)/100 +
      (((minmatchpercentage * 128) % 100 == 0) ? 0 : 1);
  const uint64_t max_history_mask
    = max_history == 64 ? (~((uint64_t) 0))
                        : ((((uint64_t) 1) << max_history) - 1);

  ft_sequenceobject_init(&useq,
                         ufsr->extend_char_access,
                         ufsr->encseq,
                         rightextension,
                         ufsr->readmode,
                         0,
                         ustart,
                         ulen,
                         ufsr->encseq_r,
                         ufsr->sequence_cache,
                         ufsr->bytesequence,
                         ufsr->totallength);
  ft_sequenceobject_init(&vseq,
                         vfsr->extend_char_access,
                         vfsr->encseq,
                         rightextension,
                         vfsr->readmode,
                         vseqstartpos,
                         vstart,
                         vlen,
                         vfsr->encseq_r,
                         vfsr->sequence_cache,
                         vfsr->bytesequence,
                         vfsr->totallength);
  frontspace->offset = 0;
  for (distance = 0, valid = 1UL; /* Nothing */; distance++, valid += 2)
  {
    GtUword trim, maxalignedlen, minlenfrommaxdiff;

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
      if (seedlength >= sizeof (validbasefront->matchhistory_bits) * CHAR_BIT)
      {
        validbasefront->matchhistory_bits = ~((uint64_t) 0);
      } else
      {
        validbasefront->matchhistory_bits = (((uint64_t) 1) << seedlength) - 1;
      }
      validbasefront->matchhistory_size = MIN(max_history,seedlength);
      validbasefront->backreference = 0; /* No back reference */
      validbasefront->max_mismatches = 0;
      front_prune_add_matches(validbasefront + distance,validbasefront,
                              max_history_mask,max_history,
                              &useq,&vseq,trimstat);
      maxalignedlen = GT_MULT2(validbasefront->row);
    } else
    {
      gt_assert(valid >= 3UL);
      frontspace_check((const GtFtFrontvalue *) frontspace->space,
                       ((const GtFtFrontvalue *) frontspace->space)
                        + frontspace->allocated - 1,
                       validbasefront + trimleft);
      frontspace_check((const GtFtFrontvalue *) frontspace->space,
                       ((const GtFtFrontvalue *) frontspace->space)
                         + frontspace->allocated - 1,
                       validbasefront + trimleft + valid - 1);
      if (valid == 3UL)
      {
        maxalignedlen = front_second_inplace(validbasefront + distance,
                                             validbasefront + trimleft,
                                             max_history,
                                             max_history_mask,
                                             &useq,
                                             &vseq,
                                             trimstat);
      } else
      {
        maxalignedlen
          = front_next_inplace(validbasefront + distance,
                               validbasefront + trimleft,
                               validbasefront + trimleft + valid - 1,
                               max_history,
                               max_history_mask,
                               &useq,
                               &vseq,
                               trimstat);
      }
    }
    gt_assert(valid > 0);
    minlenfrommaxdiff = maxalignedlen >= maxalignedlendifference
                          ? maxalignedlen - maxalignedlendifference
                          : 0;
    trim = trim_front(true,
                      distance,
                      ulen,
                      vlen,
                      minmatchpercentage128,
                      minlenfrommaxdiff,
                      trimstrategy,
                      best_polished_point,
                      validbasefront + distance,
                      validbasefront + trimleft,
                      validbasefront + trimleft + valid,
                      max_history_mask,
                      showfrontinfo);
    if (trim > 0)
    {
      trimleft += trim;
      gt_assert(valid >= trim);
      valid -= trim;
    }
    if (valid > 0)
    {
      trim = trim_front(false,
                        distance,
                        ulen,
                        vlen,
                        minmatchpercentage128,
                        minlenfrommaxdiff,
                        trimstrategy,
                        best_polished_point,
                        validbasefront + distance,
                        validbasefront + trimleft + valid - 1,
                        validbasefront + trimleft - 1,
                        max_history_mask,
                        showfrontinfo);
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
    ft_update_trace_and_polished(best_polished_point,
                                 front_trace,
                                 pol_info,
                                 distance,
                                 trimleft,
                                 validbasefront + distance,
                                 validbasefront + trimleft,
                                 validbasefront + trimleft + valid - 1,
                                 showfrontinfo);
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
  if (trimstat != NULL)
  {
    gt_ft_trimstat_add(trimstat,diedout,sumvalid,maxvalid,distance,
                       sizeof (GtFtFrontvalue) * frontspace->allocated,
                       (useq.sequence_cache != NULL &&
                        vseq.sequence_cache != NULL)
                         ? MAX(useq.sequence_cache->allocated,
                               vseq.sequence_cache->allocated)
                         : 0
                      );
  }
  return diedout ? sumseqlength + 1 : distance;
}
