#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "core/assert_api.h"
#include "core/arraydef.h"
#include "match/ft-polish.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/ma_api.h"
#include "ft-eoplist.h"

#define DELETION_CHAR    'D'
#define INSERTION_CHAR   'I'
#define MATCH_CHAR       '='
#define MISMATCH_CHAR    'X'
#define REPLACEMENT_CHAR 'M'

char gt_eoplist_pretty_print(GtEopType eoptype,bool distinguish_mismatch_match)
{
  switch (eoptype)
  {
    case GtDeletionOp:
      return DELETION_CHAR;
    case GtInsertionOp:
      return INSERTION_CHAR;
    case GtMismatchOp:
      return distinguish_mismatch_match ? MISMATCH_CHAR : REPLACEMENT_CHAR;
    case GtMatchOp:
      return distinguish_mismatch_match ? MATCH_CHAR : REPLACEMENT_CHAR;
    default:
      fprintf(stderr,"file %s, line %d: illegal eoptype = %d\n",
              __FILE__,__LINE__,(int) eoptype);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

struct GtEoplist
{
  GtUword nextfreeuint8_t, allocateduint8_t, countmismatches, countmatches,
                                             countdeletions, countinsertions;
  uint8_t *spaceuint8_t;
  const GtUchar *useq, *vseq;
  GtUword ulen, vlen;
  bool withpolcheck, seed_display;
  GtUword useedoffset, seedlen;
#ifndef OUTSIDE_OF_GT
  const Polishing_info *pol_info;
#endif
};

void gt_eoplist_reset(GtEoplist *eoplist)
{
  if (eoplist != NULL)
  {
    eoplist->nextfreeuint8_t = 0;
    eoplist->countmismatches = 0;
    eoplist->countmatches = 0;
    eoplist->countdeletions = 0;
    eoplist->countinsertions = 0;
    eoplist->useq = eoplist->vseq = NULL;
    eoplist->ulen = eoplist->vlen = 0;
  }
}

GtEoplist *gt_eoplist_new(void)
{
  GtEoplist *eoplist = gt_malloc(sizeof *eoplist);

  gt_assert(eoplist != NULL);
  eoplist->allocateduint8_t = 0;
  eoplist->spaceuint8_t = NULL;
  eoplist->seed_display = false;
  eoplist->withpolcheck = false;
  eoplist->useedoffset = eoplist->seedlen = 0;
#ifndef OUTSIDE_OF_GT
  eoplist->pol_info = NULL;
#endif
  gt_eoplist_reset(eoplist);
  return eoplist;
}

typedef struct
{
  char *space;
  size_t nextfree, allocated;
} GtStringBuffer;

static void stringbuffer_append_cigar(GtStringBuffer *sbuf,
                                      const GtCigarOp *co,
                                      bool distinguish_mismatch_match)
{
  const size_t gt_uword_maxwidth = sizeof ("18446744073709551615");

  if (sbuf->nextfree + gt_uword_maxwidth + 1 + 1 >= sbuf->allocated)
  {
    sbuf->allocated = sbuf->allocated * 1.2 + gt_uword_maxwidth + 1 + 1 + 1;
    sbuf->space = gt_realloc(sbuf->space,sizeof *sbuf->space * sbuf->allocated);
  }
  sbuf->nextfree +=
    sprintf(sbuf->space + sbuf->nextfree,"" GT_WU "%c",
            co->iteration,gt_eoplist_pretty_print(co->eoptype,
                                                  distinguish_mismatch_match));
}

char *gt_eoplist2cigar_string(const GtEoplist *eoplist,
                              bool distinguish_mismatch_match)
{
  GtStringBuffer sbuf = {NULL,0,0};
  GtCigarOp co;
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);

  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    stringbuffer_append_cigar(&sbuf,&co,distinguish_mismatch_match);
  }
  gt_eoplist_reader_delete(eoplist_reader);
  return sbuf.space;
}

GtEoplist *gt_eoplist_new_from_cigar(const char *cigarstring,GtUword length)
{
  const char *cptr;
  GtUword iteration = 0;
  GtEoplist *eoplist = gt_eoplist_new();

  for (cptr = cigarstring; cptr < cigarstring + length; cptr++)
  {
    if (isdigit(*cptr))
    {
      iteration = iteration * 10 + (GtUword) (*cptr - '0');
    } else
    {
      GtUword idx;

      switch (*cptr)
      {
        case DELETION_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_deletion_add(eoplist);
          }
          break;
        case INSERTION_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_insertion_add(eoplist);
          }
          break;
        case MATCH_CHAR:
        case REPLACEMENT_CHAR:
          gt_eoplist_match_add(eoplist,iteration);
          break;
        case MISMATCH_CHAR:
          for (idx = 0; idx < iteration; idx++)
          {
            gt_eoplist_mismatch_add(eoplist);
          }
          break;
        default:
          fprintf(stderr,"file %s, line %d: illegal symbol '%c' "
                         "in cigar string\n",__FILE__,__LINE__,*cptr);
          exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      iteration = 0;
    }
  }
  return eoplist;
}

void gt_eoplist_delete(GtEoplist *eoplist)
{
  if (eoplist != NULL)
  {
    gt_free(eoplist->spaceuint8_t);
    gt_free(eoplist);
  }
}

#define FT_EOPCODE_MAXMATCHES     253
#define FT_EOPCODE_MISMATCH       253
#define FT_EOPCODE_DELETION       254
#define FT_EOPCODE_INSERTION      255

#ifdef OUTSIDE_OF_GT
#define GT_CHECKARRAYSPACE(A,TYPE,L)\
        do {\
          if ((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
          {\
            (A)->allocated##TYPE += L;\
            (A)->space##TYPE = (TYPE *) gt_realloc((A)->space##TYPE,\
                                                   sizeof (TYPE) *\
                                                   (A)->allocated##TYPE);\
          }\
          gt_assert((A)->space##TYPE != NULL);\
        } while (false)

#define GT_STOREINARRAY(A,TYPE,L,VAL)\
        do {\
          GT_CHECKARRAYSPACE(A,TYPE,L);\
          (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL;\
        } while (false)
#endif

#define GT_EOPLIST_PUSH(EOPLIST,EOP)\
        gt_assert((EOPLIST) != NULL);\
        {\
          const GtUword addamount = (EOPLIST)->allocateduint8_t * 0.2 + 128;\
          GT_STOREINARRAY(EOPLIST,uint8_t,addamount,(uint8_t) (EOP));\
        }

void gt_eoplist_match_add(GtEoplist *eoplist,GtUword length)
{
  gt_assert(eoplist != NULL && length > 0);
  eoplist->countmatches += length;
  while (true)
  {
    if (length <= FT_EOPCODE_MAXMATCHES)
    {
      gt_assert(length > 0);
      GT_EOPLIST_PUSH(eoplist,(uint8_t) (length - 1)); /* R length */
      break;
    }
    GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_MAXMATCHES - 1); /* R max */
    length -= FT_EOPCODE_MAXMATCHES;
  }
}

void gt_eoplist_mismatch_add(GtEoplist *eoplist)
{
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_MISMATCH); /* R 1 */
  eoplist->countmismatches++;
}

void gt_eoplist_deletion_add(GtEoplist *eoplist)
{
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_DELETION);
  eoplist->countdeletions++;
}

void gt_eoplist_insertion_add(GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  GT_EOPLIST_PUSH(eoplist,FT_EOPCODE_INSERTION);
  eoplist->countinsertions++;
}

GtUword gt_eoplist_length(const GtEoplist *eoplist)
{
  if (eoplist == NULL)
  {
    return 0;
  }
  return eoplist->nextfreeuint8_t;
}

void gt_eoplist_reverse_end(GtEoplist *eoplist,GtUword firstindex)
{
  uint8_t *fwd, *bck;

  gt_assert(eoplist != NULL);
  if (firstindex + 1 >= eoplist->nextfreeuint8_t)
  {
    return;
  }
  for (fwd = eoplist->spaceuint8_t + firstindex,
       bck = eoplist->spaceuint8_t + eoplist->nextfreeuint8_t - 1; fwd < bck;
       fwd++, bck--)
  {
    uint8_t tmp = *fwd;
    *fwd = *bck;
    *bck = tmp;
  }
}

GtUword gt_eoplist_matches_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countmatches;
}

GtUword gt_eoplist_mismatches_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countmismatches;
}

GtUword gt_eoplist_deletions_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countdeletions;
}

GtUword gt_eoplist_insertions_count(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return eoplist->countinsertions;
}

GtUword gt_eoplist_unit_cost(const GtEoplist *eoplist)
{
  return gt_eoplist_insertions_count(eoplist) +
         gt_eoplist_deletions_count(eoplist) +
         gt_eoplist_mismatches_count(eoplist);
}

struct GtEoplistReader
{
  const uint8_t *endeoplist,
                *currenteop;
  bool distinguish_mismatch_match;
  GtUchar *outbuffer;
  unsigned int width;
  GtUword aligned_u, aligned_v, repcount;
};

void gt_eoplist_reader_reset(GtEoplistReader *eoplist_reader,
                             const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  gt_assert(eoplist_reader != NULL);
  if (eoplist->spaceuint8_t == NULL || eoplist->nextfreeuint8_t == 0)
  {
    eoplist_reader->currenteop = NULL;
    eoplist_reader->endeoplist = NULL;
  } else
  {
    eoplist_reader->currenteop = eoplist->spaceuint8_t;
    eoplist_reader->endeoplist = eoplist->spaceuint8_t +
                                 eoplist->nextfreeuint8_t;
  }
  eoplist_reader->aligned_u = eoplist_reader->aligned_v
                            = eoplist_reader->repcount = 0;
}

void gt_eoplist_reader_reset_width(GtEoplistReader *eoplist_reader,
                                   unsigned int width)
{
  if (eoplist_reader->width < width)
  {
    eoplist_reader->width = width;
    eoplist_reader->outbuffer = gt_realloc(eoplist_reader->outbuffer,
                                           sizeof *eoplist_reader->outbuffer *
                                           3 * (width+1));
    gt_assert(eoplist_reader->outbuffer != NULL);
  }
}

GtEoplistReader *gt_eoplist_reader_new(const GtEoplist *eoplist)
{
  GtEoplistReader *eoplist_reader;

  gt_assert(eoplist != NULL);
  eoplist_reader = gt_malloc(sizeof *eoplist_reader);
  gt_assert(eoplist_reader != NULL);
  eoplist_reader->width = 0;
  eoplist_reader->outbuffer = NULL;
  eoplist_reader->distinguish_mismatch_match = false;
  eoplist_reader->width = 70;
  eoplist_reader->outbuffer = gt_realloc(eoplist_reader->outbuffer,
                                         sizeof *eoplist_reader->outbuffer *
                                         3 * (eoplist_reader->width+1));
  gt_eoplist_reader_reset(eoplist_reader,eoplist);
  return eoplist_reader;
}

void gt_eoplist_reader_distinguish_mismatch_match(
      GtEoplistReader *eoplist_reader)
{
  eoplist_reader->distinguish_mismatch_match = true;
}

void gt_eoplist_reader_delete(GtEoplistReader *eoplist_reader)
{
  if (eoplist_reader != NULL)
  {
    if (eoplist_reader->outbuffer != NULL)
    {
      gt_free(eoplist_reader->outbuffer);
    }
    gt_free(eoplist_reader);
  }
}

bool gt_eoplist_reader_next_cigar(GtCigarOp *cigar_op,
                                  GtEoplistReader *eoplist_reader)
{
  if (eoplist_reader->currenteop == NULL ||
      eoplist_reader->currenteop >= eoplist_reader->endeoplist)
  {
    return false;
  }
  cigar_op->eoptype = GtUndefinedOp;
  cigar_op->iteration = 0;
  while (true)
  {
    if (cigar_op->iteration > 0)
    {
      switch (*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (cigar_op->eoptype == GtDeletionOp)
          {
            cigar_op->iteration++;
            eoplist_reader->currenteop++;
            break;
          }
          return true;
        case FT_EOPCODE_INSERTION:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (cigar_op->eoptype == GtInsertionOp)
          {
            cigar_op->iteration++;
            eoplist_reader->currenteop++;
            break;
          }
          return true;
       case FT_EOPCODE_MISMATCH:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (eoplist_reader->distinguish_mismatch_match)
          {
            if (cigar_op->eoptype == GtMismatchOp)
            {
              cigar_op->iteration++;
              eoplist_reader->currenteop++;
              break;
            } else
            {
              return true;
            }
          } else
          {
            if (cigar_op->eoptype == GtMatchOp)
            {
              cigar_op->iteration++;
              eoplist_reader->currenteop++;
              break;
            } else
            {
              return true;
            }
          }
          gt_assert(false);
       default:
          if (cigar_op->eoptype == GtMatchOp)
          {
            gt_assert(*eoplist_reader->currenteop < FT_EOPCODE_MAXMATCHES);
            cigar_op->iteration += (1UL + *eoplist_reader->currenteop);
            eoplist_reader->currenteop++;
          } else
          {
            return true;
          }
      }
    } else
    {
      switch (*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          cigar_op->eoptype = GtDeletionOp;
          cigar_op->iteration = 1UL;
          break;
        case FT_EOPCODE_INSERTION:
          cigar_op->eoptype = GtInsertionOp;
          cigar_op->iteration = 1UL;
          break;
        case FT_EOPCODE_MISMATCH:
          cigar_op->eoptype
            = eoplist_reader->distinguish_mismatch_match ? GtMismatchOp
                                                         : GtMatchOp;
          cigar_op->iteration = 1UL;
          break;
        default:
          cigar_op->eoptype = GtMatchOp;
          cigar_op->iteration = (1UL + *eoplist_reader->currenteop);
          break;
      }
      eoplist_reader->currenteop++;
    }
    if (eoplist_reader->currenteop >= eoplist_reader->endeoplist)
    {
      return true;
    }
  }
}

GtUword gt_eoplist_num_segments(const GtEoplist *eoplist,GtUword delta)
{
  GtUword length_u = gt_eoplist_matches_count(eoplist) +
                     gt_eoplist_mismatches_count(eoplist) +
                     gt_eoplist_deletions_count(eoplist);

  return length_u/delta + 1;
}

bool gt_eoplist_reader_next_segment(GtEoplistSegment *segment,
                                    GtEoplistReader *eoplist_reader,
                                    GtUword delta)
{
  while (true)
  {
    if (eoplist_reader->repcount > 0)
    {
      eoplist_reader->aligned_u++;
      eoplist_reader->aligned_v++;
      eoplist_reader->repcount--;
    } else
    {
      if (eoplist_reader->currenteop >= eoplist_reader->endeoplist)
      {
        break;
      }
      switch (*eoplist_reader->currenteop)
      {
        case FT_EOPCODE_DELETION:
          eoplist_reader->aligned_u++;
          break;
        case FT_EOPCODE_INSERTION:
          eoplist_reader->aligned_v++;
          break;
        case FT_EOPCODE_MISMATCH:
          eoplist_reader->aligned_u++;
          eoplist_reader->aligned_v++;
          break;
        default:
          eoplist_reader->aligned_u++;
          eoplist_reader->aligned_v++;
          eoplist_reader->repcount = (GtUword) *eoplist_reader->currenteop;
      }
      eoplist_reader->currenteop++;
    }
    if (eoplist_reader->aligned_u == delta)
    {
      segment->aligned_u = eoplist_reader->aligned_u;
      segment->aligned_v = eoplist_reader->aligned_v;
      eoplist_reader->aligned_u = eoplist_reader->aligned_v = 0;
      return true;
    }
  }
  if (eoplist_reader->aligned_v > 0 || eoplist_reader->aligned_u > 0)
  {
    gt_assert(eoplist_reader->repcount == 0);
    segment->aligned_u = eoplist_reader->aligned_u;
    segment->aligned_v = eoplist_reader->aligned_v;
    eoplist_reader->aligned_v = 0;
    eoplist_reader->aligned_u = 0;
    gt_assert(eoplist_reader->repcount == 0 &&
           eoplist_reader->currenteop >= eoplist_reader->endeoplist);
    return true;
  }
  return false;
}

double gt_eoplist_segments_entropy(const GtEoplist *eoplist,GtUword delta)
{
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);
  GtEoplistSegment segment;
  const GtUword max_value = 2 * delta + 1;
  GtUword segment_count = 0, idx,
          *segment_dist = gt_calloc(max_value + 1,sizeof *segment_dist);
  double entropy = 0.0;

  while (gt_eoplist_reader_next_segment(&segment,eoplist_reader,delta))
  {
    gt_assert(segment.aligned_v <= max_value);
    segment_dist[segment.aligned_v]++;
    segment_count++;
  }
  gt_assert(segment_count <= gt_eoplist_num_segments(eoplist,delta));
  gt_eoplist_reader_delete(eoplist_reader);
  for (idx = 0; idx <= max_value; idx++)
  {
    if (segment_dist[idx] > 0)
    {
      double prob = (double) segment_dist[idx]/segment_count;
      entropy += prob * log2(prob);
    }
  }
  gt_free(segment_dist);
  return entropy == 0.0 ? 0.0 : -entropy;
}

void gt_eoplist_show_plain(const GtEoplist *eoplist,FILE *fp)
{
  GtUword idx;

  fprintf(fp,"[");
  for (idx = 0; idx < eoplist->nextfreeuint8_t; idx++)
  {
    if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_DELETION)
    {
      fputc(DELETION_CHAR,fp);
    } else
    {
      if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_INSERTION)
      {
        fputc(INSERTION_CHAR,fp);
      } else
      {
        if (eoplist->spaceuint8_t[idx] == FT_EOPCODE_MISMATCH)
        {
          fputc(MISMATCH_CHAR,fp);
        } else
        {
          fprintf(fp,"%d",eoplist->spaceuint8_t[idx]);
        }
      }
    }
    fputc(idx + 1 < eoplist->nextfreeuint8_t ? ',' : ']',fp);
  }
  fputc('\n',fp);
}

static unsigned int gt_eoplist_show_advance(unsigned int pos,
                                            unsigned int width,
                                            const GtUchar *topbuf,
                                            FILE *fp)
{
  gt_assert(width > 0);
  if (pos < width - 1)
  {
    return pos + 1;
  }
  gt_assert(pos == width - 1);
  fwrite(topbuf,sizeof *topbuf,3 * (width+1),fp);
  return 0;
}

void gt_eoplist_set_sequences(GtEoplist *eoplist,
                              const GtUchar *useq,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vlen)
{
  gt_assert(eoplist != NULL);
  eoplist->useq = useq;
  eoplist->vseq = vseq;
  eoplist->ulen = ulen;
  eoplist->vlen = vlen;
}

#define EOPLIST_MATCHSYMBOL    '|'
#define EOPLIST_MISMATCHSYMBOL ' '
#define EOPLIST_GAPSYMBOL      '-'

#ifndef OUTSIDE_OF_GT

#define GT_UPDATE_POSITIVE_INFO(ISMATCH)\
        if (eoplist->pol_info != NULL)\
        {\
          if (prefix_positive < max_history && prefix_positive_sum >= 0)\
          {\
            if (ISMATCH)\
            {\
              prefix_positive_sum += eoplist->pol_info->match_score;\
            } else\
            {\
              prefix_positive_sum -= eoplist->pol_info->difference_score;\
            }\
            if (prefix_positive_sum >= 0)\
            {\
              prefix_positive++;\
            }\
          }\
          if (suffix_bits_used < max_history)\
          {\
            suffix_bits_used++;\
          }\
          suffix_bits >>= 1;\
          if (ISMATCH)\
          {\
            suffix_bits |= set_mask;\
          }\
        }
#else
#define GT_UPDATE_POSITIVE_INFO(ISMATCH) /* Nothing */
#define ISSPECIAL(CC) ((CC) >= 254)
#endif

void gt_eoplist_format_generic(FILE *fp,
                               const GtEoplist *eoplist,
                               GtEoplistReader *eoplist_reader,
                               bool distinguish_mismatch_match,
                               const GtUchar *characters,
                               GtUchar wildcardshow)
{
  GtCigarOp co;
  unsigned int pos = 0;
  GtUword idx_u = 0, idx_v = 0, alignmentlength = 0,
          firstseedcolumn = GT_UWORD_MAX;
  GtUchar *topbuf = eoplist_reader->outbuffer, *midbuf = NULL, *lowbuf = NULL;
#ifndef OUTSIDE_OF_GT
  uint64_t suffix_bits = 0, set_mask = 0;
  GtUword suffix_bits_used = 0, prefix_positive = 0, pol_size = 0,
          lastseedcolumn = GT_UWORD_MAX;
  const GtUword max_history = 64;
  GtWord prefix_positive_sum = 0;

  if (eoplist->pol_info != NULL)
  {
    pol_size = GT_MULT2(eoplist->pol_info->cut_depth);
    set_mask = ((uint64_t) 1) << (max_history - 1);
  }
#endif
  gt_assert(eoplist_reader != NULL);
  topbuf[eoplist_reader->width] = '\n';
  midbuf = topbuf + eoplist_reader->width + 1;
  midbuf[eoplist_reader->width] = '\n';
  lowbuf = midbuf + eoplist_reader->width + 1;
  lowbuf[eoplist_reader->width] = '\n';
  gt_eoplist_reader_reset(eoplist_reader,eoplist);
  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    switch (co.eoptype)
    {
      GtUword j;
      GtUchar cc_a, cc_b;

      case GtMatchOp:
      case GtMismatchOp:
        for (j = 0; j < co.iteration && idx_u < eoplist->ulen &&
                                        idx_v < eoplist->vlen; j++)
        {
          cc_a = eoplist->useq[idx_u];
          cc_b = eoplist->vseq[idx_v];
          bool is_match;

          if (characters != NULL)
          {
            topbuf[pos] = ISSPECIAL(cc_a) ? wildcardshow : characters[cc_a];
            lowbuf[pos] = ISSPECIAL(cc_b) ? wildcardshow : characters[cc_b];
            is_match = (cc_a == cc_b && !ISSPECIAL(cc_a)) ? true : false;
          } else
          {
            topbuf[pos] = cc_a;
            is_match = (cc_a == cc_b) ? true : false;
            lowbuf[pos] = cc_b;
          }
          if (is_match)
          {
            if (eoplist->useedoffset <= idx_u &&
                idx_u < eoplist->useedoffset + eoplist->seedlen)
            {
              if (eoplist->seed_display)
              {
                midbuf[pos] = (GtUchar) '+';
              } else
              {
                midbuf[pos] = (GtUchar) EOPLIST_MATCHSYMBOL;
              }
              if (firstseedcolumn == GT_UWORD_MAX)
              {
                firstseedcolumn = alignmentlength;
              }
#ifndef OUTSIDE_OF_GT
              lastseedcolumn = alignmentlength;
#endif
            } else
            {
              midbuf[pos] = (GtUchar) EOPLIST_MATCHSYMBOL;
            }
          } else
          {
            midbuf[pos] = (GtUchar) EOPLIST_MISMATCHSYMBOL;
          }
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(is_match);
          alignmentlength++;
          idx_u++;
          idx_v++;
        }
        break;
      case GtDeletionOp:
        for (j = 0; j < co.iteration && idx_u < eoplist->ulen; j++)
        {
          cc_a = eoplist->useq[idx_u++];
          if (characters != NULL)
          {
            topbuf[pos] = ISSPECIAL(cc_a) ? wildcardshow : characters[cc_a];
          } else
          {
            topbuf[pos] = cc_a;
          }
          midbuf[pos] = EOPLIST_MISMATCHSYMBOL;
          lowbuf[pos] = EOPLIST_GAPSYMBOL;
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
        }
        break;
      case GtInsertionOp:
        for (j = 0; j < co.iteration && idx_v < eoplist->vlen; j++)
        {
          cc_b = eoplist->vseq[idx_v++];

          topbuf[pos] = EOPLIST_GAPSYMBOL;
          midbuf[pos] = EOPLIST_MISMATCHSYMBOL;
          if (characters != NULL)
          {
            lowbuf[pos] = ISSPECIAL(cc_b) ? wildcardshow : characters[cc_b];
          } else
          {
            lowbuf[pos] = cc_b;
          }
          pos = gt_eoplist_show_advance(pos,eoplist_reader->width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
        }
        break;
      default:
        fprintf(stderr,"file %s, line %d: illegal eoptype %d\n",
                       __FILE__,__LINE__,co.eoptype);
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  if (pos > 0)
  {
    topbuf[pos] = '\n';
    fwrite(topbuf,sizeof *topbuf,pos+1,fp);
    midbuf[pos] = '\n';
    fwrite(midbuf,sizeof *midbuf,pos+1,fp);
    lowbuf[pos] = '\n';
    fwrite(lowbuf,sizeof *lowbuf,pos+1,fp);
  }
#ifndef OUTSIDE_OF_GT
  if (eoplist->pol_info != NULL)
  {
    GtUword suffix_positive;
    GtWord suffix_positive_sum = 0;
    bool startpolished = false, endpolished = false;

    for (suffix_positive = 0; suffix_positive < suffix_bits_used;
         suffix_positive++)
    {
      suffix_positive_sum += ((suffix_bits & set_mask)
                                ? eoplist->pol_info->match_score
                                : -eoplist->pol_info->difference_score);
      if (suffix_positive_sum < 0)
      {
        break;
      }
      set_mask >>= 1;
    }
    gt_assert(prefix_positive <= alignmentlength);
    if (prefix_positive >= pol_size || prefix_positive == alignmentlength ||
        firstseedcolumn < pol_size)
    {
      startpolished = true;
    }
    if (suffix_positive >= pol_size || suffix_positive == alignmentlength ||
        (lastseedcolumn != GT_UWORD_MAX &&
         lastseedcolumn + pol_size > alignmentlength))
    {
      endpolished = true;
    }
    fprintf(fp, "# polishing(m=" GT_WD ",d=" GT_WD ",p=" GT_WU
            "): " GT_WU "/" GT_WU,
            eoplist->pol_info->match_score,
            -eoplist->pol_info->difference_score,
            pol_size,
            prefix_positive,
            suffix_positive);
    if (firstseedcolumn < pol_size)
    {
      fprintf(fp, ", seed_on_start");
    }
    if (lastseedcolumn + pol_size > alignmentlength)
    {
      fprintf(fp, ", seed_on_end");
    }
    if (eoplist->withpolcheck)
    {
      fprintf(fp, "\n");
      gt_assert(startpolished);
      gt_assert(endpolished);
    } else
    {
      if (!startpolished)
      {
        fprintf(fp, ", start not polished");
      }
      if (!endpolished)
      {
        fprintf(fp, ", end not polished");
      }
      fprintf(fp, "\n");
    }
  }
#endif
}

void gt_eoplist_format_exact(FILE *fp,
                             const GtEoplist *eoplist,
                             GtEoplistReader *eoplist_reader,
                             const GtUchar *characters)
{
  GtUword idx;
  unsigned int pos = 0, width;
  GtUchar *topbuf = eoplist_reader->outbuffer, *midbuf = NULL, *lowbuf = NULL;

  width = MIN(eoplist->ulen, eoplist_reader->width);
  topbuf[width] = '\n';
  midbuf = topbuf + width + 1;
  for (idx = 0; idx < (GtUword) width; idx++)
  {
    midbuf[idx] = (GtUchar) EOPLIST_MATCHSYMBOL;
  }
  midbuf[width] = '\n';
  lowbuf = midbuf + width + 1;
  lowbuf[width] = '\n';
  for (idx = 0; idx < eoplist->ulen; idx++)
  {
    GtUchar cc_a = eoplist->useq[idx];
    if (characters != NULL)
    {
      cc_a = characters[cc_a];
    }
    lowbuf[pos] = topbuf[pos] = cc_a;
    pos = gt_eoplist_show_advance(pos,width,topbuf,fp);
  }
  if (pos > 0)
  {
    topbuf[pos] = '\n';
    fwrite(topbuf,sizeof *topbuf,pos+1,fp);
    midbuf[pos] = '\n';
    fwrite(midbuf,sizeof *midbuf,pos+1,fp);
    lowbuf[pos] = '\n';
    fwrite(lowbuf,sizeof *lowbuf,pos+1,fp);
  }
}

void gt_eoplist_verify(const GtEoplist *eoplist,
                       GtEoplistReader *eoplist_reader,
                       GtUword edist,
                       bool distinguish_mismatch_match)
{
  GtCigarOp co;
  GtUword sumulen = 0, sumvlen = 0, sumdist = 0;

  gt_assert(eoplist != NULL);
  gt_eoplist_reader_reset(eoplist_reader,eoplist);
  if (distinguish_mismatch_match)
  {
    gt_eoplist_reader_distinguish_mismatch_match(eoplist_reader);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    if (co.eoptype == GtDeletionOp)
    {
      sumulen += co.iteration;
      sumdist += co.iteration;
    } else
    {
      if (co.eoptype == GtInsertionOp)
      {
        sumvlen += co.iteration;
        sumdist += co.iteration;
      } else
      {
        if (co.eoptype == GtMismatchOp)
        {
          sumdist += co.iteration;
        }
        if (eoplist->useq == NULL && eoplist->vseq == NULL)
        {
          gt_assert(eoplist_reader->distinguish_mismatch_match);
        } else
        {
          GtUword idx;

          for (idx = 0; idx < co.iteration; idx++)
          {
            GtUchar a = eoplist->useq[sumulen+idx],
                    b = eoplist->vseq[sumvlen+idx];
            if (a == b && !ISSPECIAL(a))
            {
              gt_assert(co.eoptype == GtMatchOp);
            } else
            {
              gt_assert(!eoplist_reader->distinguish_mismatch_match ||
                        co.eoptype == GtMismatchOp);
            }
            if (!eoplist_reader->distinguish_mismatch_match &&
                (a != b || ISSPECIAL(a)))
            {
              sumdist++;
            }
          }
        }
        sumulen += co.iteration;
        sumvlen += co.iteration;
      }
    }
  }
  if (eoplist->ulen != sumulen)
  {
    fprintf(stderr,"ulen = " GT_WU " != " GT_WU " = sumulen\n",
            eoplist->ulen,sumulen);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (eoplist->vlen != sumvlen)
  {
    fprintf(stderr,"vlen = " GT_WU " != " GT_WU " = sumvlen\n",
            eoplist->vlen,sumvlen);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (edist != sumdist)
  {
    fprintf(stderr,"edist = " GT_WU " != " GT_WU " = sumdist\n",
            edist,sumdist);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(eoplist->ulen + eoplist->vlen
            == 2 * (gt_eoplist_matches_count(eoplist) + edist) -
               (gt_eoplist_deletions_count(eoplist) +
                gt_eoplist_insertions_count(eoplist)));
}

void gt_eoplist_seed_display_set(GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  eoplist->seed_display = true;
}

void gt_eoplist_set_seedoffset(GtEoplist *eoplist,
                               GtUword useedoffset,
                               GtUword seedlen)
{
  eoplist->useedoffset = useedoffset;
  eoplist->seedlen = seedlen;
}

#ifndef OUTSIDE_OF_GT
void gt_eoplist_polished_ends(GtEoplist *eoplist,
                              const Polishing_info *pol_info,
                              bool withpolcheck)
{
  gt_assert(eoplist != NULL);
  eoplist->pol_info = pol_info;
  eoplist->withpolcheck = withpolcheck;
}
#endif
