#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "core/assert_api.h"
#include "core/arraydef_api.h"
#include "core/minmax_api.h"
#include "core/chardef_api.h"
#include "core/divmodmul_api.h"
#include "core/readmode.h"
#include "match/ft-polish.h"
#include "match/ft-eoplist.h"
#include "match/ft-front-prune.h"

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
  GtUword nextfreeuint8_t, allocateduint8_t, countdeletions, countinsertions;
  uint8_t *spaceuint8_t;
  const GtUchar *useq, *vseq;
  GtUword ustart, ulen, vstart, vlen;
  bool withpolcheck, pol_info_out, display_seed_in_alignment;
  GtUword useedoffset, seedlen;
  GtArrayint trace;
  const GtFtPolishing_info *pol_info;
  GtFullFrontEdistTrace *fet_segment;
};

void gt_eoplist_reset(GtEoplist *eoplist)
{
  if (eoplist != NULL)
  {
    eoplist->nextfreeuint8_t = 0;
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
  eoplist->display_seed_in_alignment = false;
  eoplist->withpolcheck = false;
  eoplist->pol_info_out = false;
  eoplist->useedoffset = eoplist->seedlen = 0;
  eoplist->pol_info = NULL;
  GT_INITARRAY(&eoplist->trace,int);
  eoplist->fet_segment = gt_full_front_edist_trace_new();
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
    sprintf(sbuf->space + sbuf->nextfree,GT_WU "%c",
            co->iteration,gt_eoplist_pretty_print(co->eoptype,
                                                  distinguish_mismatch_match));
}

char *gt_eoplist2cigar_string(const GtEoplist *eoplist,
                              bool distinguish_mismatch_match)
{
  GtStringBuffer sbuf = {NULL,0,0};
  GtCigarOp co;
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new();

  gt_eoplist_reader_reset(eoplist_reader,eoplist,true);
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader,
                                      distinguish_mismatch_match))
  {
    stringbuffer_append_cigar(&sbuf,&co,distinguish_mismatch_match);
  }
  gt_eoplist_reader_delete(eoplist_reader);
  return sbuf.space;
}

void gt_eoplist_from_cigar(GtEoplist *eoplist,
                           const char *cigarstring,char sep)
{
  const char *cptr;
  GtUword iteration = 0;

  for (cptr = cigarstring; *cptr != '\0' && *cptr != sep && *cptr != '\n';
       cptr++)
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
          gt_eoplist_match_add(eoplist,iteration);
          break;
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
}

void gt_eoplist_delete(GtEoplist *eoplist)
{
  if (eoplist != NULL)
  {
    GT_FREEARRAY(&eoplist->trace,int);
    gt_full_front_edist_trace_delete(eoplist->fet_segment);
    gt_free(eoplist->spaceuint8_t);
    gt_free(eoplist);
  }
}

#define FT_EOPCODE_MAXMATCHES     253
#define FT_EOPCODE_MISMATCH       253
#define FT_EOPCODE_DELETION       254
#define FT_EOPCODE_INSERTION      255

#define GT_EOPLIST_PUSH(EOPLIST,EOP)\
        gt_assert((EOPLIST) != NULL);\
        {\
          const GtUword addamount = (EOPLIST)->allocateduint8_t * 0.2 + 128;\
          GT_STOREINARRAY(EOPLIST,uint8_t,addamount,(uint8_t) (EOP));\
        }

void gt_eoplist_match_add(GtEoplist *eoplist,GtUword length)
{
  gt_assert(eoplist != NULL && length > 0);
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

struct GtEoplistReader
{
  const uint8_t *currenteop,
                *endeoplist;
  int difference;
  GtUchar *outbuffer;
  unsigned int width;
  GtUword aligned_u, aligned_v, repcount;
};

void gt_eoplist_reader_reset(GtEoplistReader *eoplist_reader,
                             const GtEoplist *eoplist,bool forward)
{
  gt_assert(eoplist != NULL && eoplist_reader != NULL);
  if (eoplist->spaceuint8_t == NULL || eoplist->nextfreeuint8_t == 0)
  {
    eoplist_reader->currenteop = NULL;
    eoplist_reader->endeoplist = NULL;
  } else
  {
    if (forward)
    {
      eoplist_reader->currenteop = eoplist->spaceuint8_t;
      eoplist_reader->endeoplist = eoplist->spaceuint8_t +
                                   eoplist->nextfreeuint8_t;
      eoplist_reader->difference = 1;
    } else
    {
      eoplist_reader->currenteop = eoplist->spaceuint8_t +
                                   eoplist->nextfreeuint8_t - 1;
      eoplist_reader->endeoplist = eoplist->spaceuint8_t - 1;
      eoplist_reader->difference = -1;
    }
  }
  eoplist_reader->aligned_u = eoplist_reader->aligned_v
                            = eoplist_reader->repcount = 0;
}

static size_t gt_eoplist_outbuffer_size(unsigned int width)
{
  return 3 * width;
}

void gt_eoplist_reader_reset_width(GtEoplistReader *eoplist_reader,
                                   unsigned int width)
{
  if (eoplist_reader->width < width)
  {
    eoplist_reader->outbuffer = gt_realloc(eoplist_reader->outbuffer,
                                           sizeof *eoplist_reader->outbuffer *
                                           gt_eoplist_outbuffer_size(width));
    gt_assert(eoplist_reader->outbuffer != NULL);
  }
  eoplist_reader->width = width;
}

GtEoplistReader *gt_eoplist_reader_new(void)
{
  GtEoplistReader *eoplist_reader;

  eoplist_reader = gt_malloc(sizeof *eoplist_reader);
  gt_assert(eoplist_reader != NULL);
  eoplist_reader->width = 0;
  eoplist_reader->outbuffer = NULL;
  eoplist_reader->width = 70;
  eoplist_reader->outbuffer
    = gt_realloc(eoplist_reader->outbuffer,
                 sizeof *eoplist_reader->outbuffer *
                 gt_eoplist_outbuffer_size(eoplist_reader->width));
  eoplist_reader->difference = 0;
  eoplist_reader->currenteop = NULL;
  eoplist_reader->endeoplist = NULL;
  eoplist_reader->aligned_u = eoplist_reader->aligned_v
                            = eoplist_reader->repcount = 0;
  return eoplist_reader;
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
                                  GtEoplistReader *eoplist_reader,
                                  bool distinguish_mismatch_match)
{
  if (eoplist_reader->currenteop == NULL ||
      eoplist_reader->currenteop == eoplist_reader->endeoplist)
  {
    return false;
  }
  gt_assert(eoplist_reader->difference == 1 ||
            eoplist_reader->difference == -1);
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
            cigar_op->iteration++; /* Add another consecutive deletion */
            eoplist_reader->currenteop += eoplist_reader->difference;
            break;
          }
          return true;
        case FT_EOPCODE_INSERTION:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (cigar_op->eoptype == GtInsertionOp)
          {
            cigar_op->iteration++; /* Add another consecutive insertion */
            eoplist_reader->currenteop += eoplist_reader->difference;
            break;
          }
          return true;
       case FT_EOPCODE_MISMATCH:
          gt_assert(cigar_op->eoptype != GtUndefinedOp);
          if (distinguish_mismatch_match)
          {
            if (cigar_op->eoptype == GtMismatchOp)
            {
              cigar_op->iteration++; /* Add another consecutive mismatch */
              eoplist_reader->currenteop += eoplist_reader->difference;
              break;
            }
            return true;
          }
          if (cigar_op->eoptype == GtMatchOp)
          {
            cigar_op->iteration++; /* Add another consecutive match */
            eoplist_reader->currenteop += eoplist_reader->difference;
            break;
          }
          return true;
       default:
          if (cigar_op->eoptype == GtMatchOp)
          {
            gt_assert(*eoplist_reader->currenteop < FT_EOPCODE_MAXMATCHES);
            cigar_op->iteration += (1UL + *eoplist_reader->currenteop);
            eoplist_reader->currenteop += eoplist_reader->difference;
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
            = distinguish_mismatch_match ? GtMismatchOp : GtMatchOp;
          cigar_op->iteration = 1UL;
          break;
        default:
          cigar_op->eoptype = GtMatchOp;
          cigar_op->iteration = (1UL + *eoplist_reader->currenteop);
          break;
      }
      eoplist_reader->currenteop += eoplist_reader->difference;
    }
    if (eoplist_reader->currenteop == eoplist_reader->endeoplist)
    {
      return true;
    }
  }
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
      if (eoplist_reader->currenteop == eoplist_reader->endeoplist)
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
      eoplist_reader->currenteop += eoplist_reader->difference;
    }
    if (eoplist_reader->aligned_u == delta)
    {
      segment->aligned_u = delta;
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
    gt_assert(eoplist_reader->currenteop == eoplist_reader->endeoplist);
    return true;
  }
  return false;
}

double gt_eoplist_segments_entropy(const GtEoplist *eoplist,GtUword delta)
{
  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new();
  GtEoplistSegment segment;
  const GtUword max_value = 2 * delta + 1;
  GtUword segment_count = 0, idx,
          *segment_dist = gt_calloc(max_value + 1,sizeof *segment_dist);
  double entropy = 0.0;

  gt_eoplist_reader_reset(eoplist_reader,eoplist,true);
  while (gt_eoplist_reader_next_segment(&segment,eoplist_reader,delta))
  {
    gt_assert(segment.aligned_v <= max_value);
    segment_dist[segment.aligned_v]++;
    segment_count++;
  }
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

/*
  A trace is a sequence of
  differences with the delta value specifying the length of a substring in the
  query sequence aligned to a sequence of length delta in the reference
  sequence. Such a pair of substrings is a segment. From the segments
  one can reconstruct an alignment by computing alignments from the
  segments and concatenating them.  This is done by the following functions.
*/

void gt_eoplist_read_trace(GtEoplist *eoplist,
                           const char *trace,
                           char separator)
{
  if (eoplist != NULL)
  {
    eoplist->trace.nextfreeint = 0;
  }
  while (true)
  {
    int value;
    const char *ptr;

    if (sscanf(trace,"%d",&value) != 1)
    {
      fprintf(stderr,"cannot read number from trace %s\n",trace);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (eoplist != NULL)
    {
      GT_STOREINARRAY(&eoplist->trace,int,
                      256 + eoplist->trace.allocatedint * 0.2,value);
    }
    for (ptr = trace; *ptr != '\0' && *ptr != separator && *ptr != ','; ptr++)
      /* Nothing */;
    if (*ptr == '\0' || *ptr == separator)
    {
      break;
    }
    trace = ptr + 1;
  }
}

void gt_eoplist_trace2cigar(GtEoplist *eoplist,bool dtrace,GtUword trace_delta)
{
  GtUword idx, offset_u = 0, offset_v = 0;

  gt_assert(eoplist != NULL && eoplist->trace.nextfreeint > 0);
  for (idx = 0; idx < eoplist->trace.nextfreeint; idx++)
  {
    GtUword this_distance, aligned_u, aligned_v;

    if (dtrace)
    {
      const int value = -(eoplist->trace.spaceint[idx] - trace_delta);
      gt_assert(value >= 0);
      aligned_v = (GtUword) value;
    } else
    {
      aligned_v = (GtUword) eoplist->trace.spaceint[idx];
    }
    gt_assert(offset_u < eoplist->ulen);
    aligned_u = GT_MIN(trace_delta,eoplist->ulen - offset_u);
    this_distance = gt_full_front_edist_trace_distance(eoplist->fet_segment,
                                                       eoplist->useq + offset_u,
                                                       aligned_u,
                                                       eoplist->vseq + offset_v,
                                                       aligned_v);
    gt_front_trace2eoplist_full_front_directed(eoplist,
                                               gt_full_front_trace_get(
                                                   eoplist->fet_segment),
                                               this_distance,
                                               eoplist->useq + offset_u,
                                               aligned_u,
                                               eoplist->vseq + offset_v,
                                               aligned_v);
    offset_u += aligned_u;
    offset_v += aligned_v;
  }
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

void gt_eoplist_show_cigar(GtEoplistReader *eoplist_reader,
                           bool distinguish_mismatch_match,FILE *fp)
{
  GtCigarOp co;

  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader,
                                      distinguish_mismatch_match))
  {
    fprintf(fp,GT_WU "%c",co.iteration,
            gt_eoplist_pretty_print(co.eoptype,distinguish_mismatch_match));
  }
}

static void gt_eoplist_single_line(const char *tag,
                                   int numwidth,
                                   unsigned int width,
                                   const GtUchar *buf,
                                   GtUword start_pos,
                                   GtUword end_pos,
                                   FILE *fp)
{
  fprintf(fp,"%s  %-*" GT_WUS "  ",tag,numwidth,start_pos);
  fwrite(buf,sizeof *buf,width,fp);
  fprintf(fp,"  " GT_WU "\n",end_pos);
}

static void gt_eoplist_middle_line(int numwidth,unsigned int width,
                                   const GtUchar *midbuf,
                                   FILE *fp)
{
  /* 5 is the length of the strings Sbjct and Query */
  fprintf(fp,"%*s",(int) (numwidth + 5 + 4),"");
  fwrite(midbuf,sizeof *midbuf,width,fp);
  fputc('\n',fp);
}

static void gt_eoplist_write_lines(GtUword one_off,
                                   bool subject_first,
                                   int numwidth,
                                   unsigned int width,
                                   const GtUchar *subject_buf,
                                   GtUword subject_seqlength,
                                   GtUword subject_start_pos,
                                   GtUword subject_end_pos,
                                   const GtUchar *midbuf,
                                   const GtUchar *query_buf,
                                   GtUword query_start_pos,
                                   GtUword query_end_pos,
                                   FILE *fp)
{
  gt_assert(numwidth > 0);
  if (subject_first)
  {
    gt_eoplist_single_line("Sbjct",numwidth,width,subject_buf,
                           subject_start_pos + one_off,
                           subject_end_pos + one_off,fp);
    gt_eoplist_middle_line(numwidth,width,midbuf,fp);
    gt_eoplist_single_line("Query",numwidth,width,query_buf,
                           query_start_pos + one_off,
                           query_end_pos + one_off,fp);
  } else
  {
    gt_eoplist_single_line("Query",numwidth,width,query_buf,
                           query_start_pos + one_off,
                           query_end_pos + one_off,fp);
    gt_eoplist_middle_line(numwidth,width,midbuf,fp);
    if (subject_seqlength == 0)
    {
      gt_eoplist_single_line("Sbjct",numwidth,width,subject_buf,
                             subject_start_pos + one_off,
                             subject_end_pos + one_off,fp);
    } else
    {
      gt_assert(subject_seqlength > subject_start_pos &&
                subject_seqlength >= subject_end_pos);
      gt_eoplist_single_line("Sbjct",numwidth,width,subject_buf,
                             subject_seqlength - 1 - subject_start_pos
                               + one_off,
                             one_off +
                             (subject_seqlength > subject_end_pos
                               ? subject_seqlength - 1 - subject_end_pos
                               : 0),
                             fp);
    }
  }
  fputc('\n',fp);
}

static unsigned int gt_eoplist_show_advance(GtUword one_off,
                                            bool subject_first,
                                            int numwidth,
                                            unsigned int pos,
                                            unsigned int width,
                                            const GtUchar *topbuf,
                                            GtUword top_seqlength,
                                            GtUword top_start_pos,
                                            GtUword top_end_pos,
                                            const GtUchar *midbuf,
                                            const GtUchar *lowbuf,
                                            GtUword low_start_pos,
                                            GtUword low_end_pos,
                                            FILE *fp)
{
  gt_assert(width > 0);
  if (pos + 1 < width)
  {
    return pos + 1;
  }
  gt_assert(pos == width - 1);
  gt_eoplist_write_lines(one_off,subject_first,
                         numwidth, width, topbuf,
                         top_seqlength, top_start_pos, top_end_pos,
                         midbuf, lowbuf, low_start_pos, low_end_pos, fp);
  return 0;
}

void gt_eoplist_set_sequences(GtEoplist *eoplist,
                              const GtUchar *useq,
                              GtUword ustart,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vstart,
                              GtUword vlen)
{
  gt_assert(eoplist != NULL);
  eoplist->useq = useq;
  eoplist->ustart = ustart;
  eoplist->ulen = ulen;
  eoplist->vseq = vseq;
  eoplist->vstart = vstart;
  eoplist->vlen = vlen;
}

#define EOPLIST_MATCHSYMBOL    '|'
#define EOPLIST_MISMATCHSYMBOL ' '
#define EOPLIST_GAPSYMBOL      '-'

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

static int gt_eoplist_numwidth(const GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  return 1 + log10((double) GT_MAX(eoplist->ustart + eoplist->ulen - 1,
                                eoplist->vstart + eoplist->vlen - 1));
}

void gt_eoplist_format_generic(FILE *fp,
                               const GtEoplist *eoplist,
                               GtEoplistReader *eoplist_reader,
                               const GtUchar *characters,
                               GtUword top_seqlength,
                               GtUword low_reference,
                               GtUword one_off,
                               bool distinguish_mismatch_match,
                               bool subject_first,
                               bool alignment_show_forward,
                               bool show_complement_characters,
                               GtUchar wildcardshow)
{
  GtCigarOp co;
  unsigned int pos = 0;
  GtUword idx_u = 0, idx_v = 0, alignmentlength = 0,
          firstseedcolumn = GT_UWORD_MAX;
  GtUchar *topbuf = eoplist_reader->outbuffer, *midbuf = NULL, *lowbuf = NULL;
  const int numwidth = gt_eoplist_numwidth(eoplist);
  const GtUword low_start_base
                  = low_reference == 0 ? eoplist->vstart
                                       : low_reference - eoplist->vstart;
  GtUword top_start_pos = eoplist->ustart,
          low_start_pos = low_start_base;

  gt_assert(alignment_show_forward || top_seqlength > 0);
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
  gt_assert(eoplist_reader != NULL);
  midbuf = topbuf + eoplist_reader->width;
  lowbuf = midbuf + eoplist_reader->width;
  gt_eoplist_reader_reset(eoplist_reader,eoplist,alignment_show_forward);
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader,
                                      distinguish_mismatch_match))
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
          bool is_match = true;

          if (alignment_show_forward)
          {
            cc_a = eoplist->useq[idx_u];
            cc_b = eoplist->vseq[idx_v];
          } else
          {
            cc_a = eoplist->useq[eoplist->ulen - 1 - idx_u];
            cc_b = eoplist->vseq[eoplist->vlen - 1 - idx_v];
          }
          if (characters != NULL)
          {
            if (GT_ISSPECIAL(cc_a))
            {
              cc_a = wildcardshow;
              is_match = false;
            } else
            {
              if (show_complement_characters)
              {
                gt_assert(cc_a < 4);
                cc_a = GT_COMPLEMENTBASE(cc_a);
              }
              cc_a = characters[cc_a];
            }
            if (GT_ISSPECIAL(cc_b))
            {
              cc_b = wildcardshow;
              is_match = false;
            } else
            {
              if (show_complement_characters)
              {
                gt_assert(cc_b < 4);
                cc_b = GT_COMPLEMENTBASE(cc_b);
              }
              cc_b = characters[cc_b];
            }
          }
          topbuf[pos] = cc_a;
          if (is_match)
          {
            is_match = cc_a == cc_b ? true : false;
          }
          lowbuf[pos] = cc_b;
          if (is_match)
          {
            if (eoplist->useedoffset <= idx_u &&
                idx_u < eoplist->useedoffset + eoplist->seedlen)
            {
              if (eoplist->display_seed_in_alignment)
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
              lastseedcolumn = alignmentlength;
            } else
            {
              midbuf[pos] = (GtUchar) EOPLIST_MATCHSYMBOL;
            }
          } else
          {
            midbuf[pos] = (GtUchar) EOPLIST_MISMATCHSYMBOL;
          }
          pos = gt_eoplist_show_advance(one_off,
                                        subject_first,
                                        numwidth,
                                        pos,
                                        eoplist_reader->width,
                                        topbuf,
                                        top_seqlength,
                                        top_start_pos,
                                        eoplist->ustart + idx_u,
                                        midbuf,
                                        lowbuf,
                                        low_start_pos,
                                        low_start_base + idx_v,
                                        fp);
          if (pos == 0)
          {
            top_start_pos = eoplist->ustart + idx_u + 1;
            low_start_pos = low_start_base + idx_v + 1;
          }
          GT_UPDATE_POSITIVE_INFO(is_match);
          alignmentlength++;
          idx_u++;
          idx_v++;
        }
        break;
      case GtDeletionOp:
        for (j = 0; j < co.iteration && idx_u < eoplist->ulen; j++)
        {
          cc_a = eoplist->useq[alignment_show_forward ? idx_u
                                                      : eoplist->ulen-1-idx_u];
          if (characters != NULL)
          {
            if (GT_ISSPECIAL(cc_a))
            {
              topbuf[pos] = wildcardshow;
            } else
            {
              if (show_complement_characters)
              {
                gt_assert(cc_a < 4);
                cc_a = GT_COMPLEMENTBASE(cc_a);
              }
              topbuf[pos] = characters[cc_a];
            }
          } else
          {
            topbuf[pos] = cc_a;
          }
          midbuf[pos] = EOPLIST_MISMATCHSYMBOL;
          lowbuf[pos] = EOPLIST_GAPSYMBOL;
          pos = gt_eoplist_show_advance(one_off,
                                        subject_first,
                                        numwidth,
                                        pos,
                                        eoplist_reader->width,
                                        topbuf,
                                        top_seqlength,
                                        top_start_pos,
                                        eoplist->ustart + idx_u,
                                        midbuf,
                                        lowbuf,
                                        low_start_pos,
                                        low_start_base + idx_v,
                                        fp);
          if (pos == 0)
          {
            top_start_pos = eoplist->ustart + idx_u + 1;
            low_start_pos = low_start_base + idx_v + 1;
          }
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
          idx_u++;
        }
        break;
      case GtInsertionOp:
        for (j = 0; j < co.iteration && idx_v < eoplist->vlen; j++)
        {
          cc_b = eoplist->vseq[alignment_show_forward ? idx_v
                                                      : eoplist->vlen-1-idx_v];
          topbuf[pos] = EOPLIST_GAPSYMBOL;
          midbuf[pos] = EOPLIST_MISMATCHSYMBOL;
          if (characters != NULL)
          {
            if (GT_ISSPECIAL(cc_b))
            {
              lowbuf[pos] = wildcardshow;
            } else
            {
              if (show_complement_characters)
              {
                gt_assert(cc_b < 4);
                cc_b = GT_COMPLEMENTBASE(cc_b);
              }
              lowbuf[pos] = characters[cc_b];
            }
          } else
          {
            lowbuf[pos] = cc_b;
          }
          pos = gt_eoplist_show_advance(one_off,
                                        subject_first,
                                        numwidth,
                                        pos,
                                        eoplist_reader->width,
                                        topbuf,
                                        top_seqlength,
                                        top_start_pos,
                                        eoplist->ustart + idx_u,
                                        midbuf,
                                        lowbuf,
                                        low_start_pos,
                                        low_start_base + idx_v,
                                        fp);
          if (pos == 0)
          {
            top_start_pos = eoplist->ustart + idx_u + 1;
            low_start_pos = low_start_base + idx_v + 1;
          }
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
          idx_v++;
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
    gt_eoplist_write_lines(one_off,
                           subject_first,
                           numwidth,
                           pos,
                           topbuf,
                           top_seqlength,
                           top_start_pos,
                           eoplist->ustart + GT_MIN(idx_u,eoplist->ulen - 1),
                           midbuf,
                           lowbuf,
                           low_start_pos,
                           low_start_base + GT_MIN(idx_v,eoplist->vlen - 1),
                           fp);
  }
  if (eoplist->pol_info != NULL && eoplist->pol_info_out)
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
}

void gt_eoplist_format_exact(FILE *fp,
                             const GtEoplist *eoplist,
                             GtEoplistReader *eoplist_reader,
                             GtUword top_seqlength,
                             GtUword low_reference,
                             GtUword one_off,
                             bool subject_first,
                             bool alignment_show_forward,
                             bool show_complement_characters,
                             const GtUchar *characters)
{
  GtUword idx;
  unsigned int pos = 0, width;
  GtUchar *topbuf = eoplist_reader->outbuffer, *midbuf = NULL, *lowbuf = NULL;
  const int numwidth = gt_eoplist_numwidth(eoplist);
  GtUword top_start_pos = eoplist->ustart,
          low_start_pos = low_reference == 0 ? eoplist->vstart
                                             : low_reference - eoplist->vstart;

  gt_assert(alignment_show_forward || top_seqlength > 0);
  width = GT_MIN(eoplist->ulen, eoplist_reader->width);
  midbuf = topbuf + width;
  for (idx = 0; idx < (GtUword) width; idx++)
  {
    midbuf[idx] = (GtUchar) EOPLIST_MATCHSYMBOL;
  }
  lowbuf = midbuf + width;
  for (idx = 0; idx < eoplist->ulen; idx++)
  {
    GtUchar cc_a = eoplist->useq[alignment_show_forward ? idx
                                                        : eoplist->ulen-1-idx];

    if (characters != NULL)
    {
      if (show_complement_characters)
      {
        cc_a = GT_COMPLEMENTBASE(cc_a);
      }
      cc_a = characters[cc_a];
    }
    lowbuf[pos] = topbuf[pos] = cc_a;
    pos = gt_eoplist_show_advance(one_off,
                                  subject_first,
                                  numwidth,
                                  pos,
                                  width,
                                  topbuf,
                                  top_seqlength,
                                  top_start_pos,
                                  eoplist->ustart + idx,
                                  midbuf,
                                  lowbuf,
                                  low_start_pos,
                                  eoplist->vstart + idx,fp);
    if (pos == 0)
    {
      top_start_pos = eoplist->ustart + idx + 1;
      low_start_pos = eoplist->vstart + idx + 1;
    }
  }
  if (pos > 0)
  {
    gt_eoplist_write_lines(one_off,
                           subject_first,
                           numwidth,
                           pos,
                           topbuf,
                           top_seqlength,
                           top_start_pos,
                           eoplist->ustart + GT_MIN(idx,eoplist->ulen - 1),
                           midbuf,
                           lowbuf,
                           low_start_pos,
                           eoplist->vstart + GT_MIN(idx,eoplist->vlen - 1),
                           fp);
  }
}

void gt_eoplist_verify(const GtEoplist *eoplist,
                       GtEoplistReader *eoplist_reader,
                       GtUword edist)
{
  GtCigarOp co;
  GtUword sumulen = 0, sumvlen = 0, sumdist = 0;
  const bool distinguish_mismatch_match = true;

  gt_assert(eoplist != NULL);
  gt_eoplist_reader_reset(eoplist_reader,eoplist,true);
  if (eoplist->useq == NULL)
  {
    gt_assert(eoplist->vseq == NULL && distinguish_mismatch_match);
  } else
  {
    gt_assert(eoplist->vseq != NULL);
  }
  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader,
                                      distinguish_mismatch_match))
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
        if (eoplist->useq != NULL)
        {
          GtUword idx;

          for (idx = 0; idx < co.iteration; idx++)
          {
            const GtUchar a = eoplist->useq[sumulen+idx],
                          b = eoplist->vseq[sumvlen+idx];
            if (a == b && !GT_ISSPECIAL(a))
            {
              gt_assert(co.eoptype == GtMatchOp);
            } else
            {
              gt_assert(!distinguish_mismatch_match ||
                        co.eoptype == GtMismatchOp);
            }
            if (!distinguish_mismatch_match && (a != b || GT_ISSPECIAL(a)))
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
}

void gt_eoplist_display_seed_in_alignment_set(GtEoplist *eoplist)
{
  gt_assert(eoplist != NULL);
  eoplist->display_seed_in_alignment = true;
}

void gt_eoplist_set_seedoffset(GtEoplist *eoplist,
                               GtUword useedoffset,
                               GtUword seedlen)
{
  eoplist->useedoffset = useedoffset;
  eoplist->seedlen = seedlen;
}

void gt_eoplist_polished_ends(GtEoplist *eoplist,
                              const GtFtPolishing_info *pol_info,
                              bool withpolcheck,
                              bool pol_info_out)
{
  gt_assert(eoplist != NULL);
  eoplist->pol_info = pol_info;
  eoplist->withpolcheck = withpolcheck;
  eoplist->pol_info_out = pol_info_out;
}
