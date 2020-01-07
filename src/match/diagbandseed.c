/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2015-2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015-2017 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <float.h>
#include <math.h>
#include "core/arraydef_api.h"
#include "core/codetype.h"
#include "core/complement.h"
#include "core/cstr_api.h"
#include "core/encseq.h"
#include "core/fa_api.h"
#include "core/fileutils_api.h"
#include "core/ma_api.h"
#include "core/minmax_api.h"
#include "core/radix_sort.h"
#include "core/timer_api.h"
#include "core/spacecalc.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "core/intbits.h"
#include "core/qsort-ulong.h"
#include "core/log_api.h"
#include "core/bittab_api.h"
#include "match/chain2dim.h"
#include "match/declare-readfunc.h"
#include "match/kmercodes.h"
#include "match/querymatch.h"
#include "match/querymatch-align.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"
#include "match/rectangle-store.h"
#include "match/diagband-struct.h"
#include "match/dbs_spaced_seeds.h"
#include "match/diagbandseed.h"

#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

/* We need to use 6 digits for the micro seconds */
#define GT_DIAGBANDSEED_FMT          "in " GT_WD ".%06ld seconds.\n"

typedef uint32_t GtDiagbandseedSeqnum;

typedef struct { /* 8 + 4 + 4 bytes */
  GtCodetype code;              /* only sort criterion */
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
} GtDiagbandseedKmerPos;

typedef struct
{
  GtBitcount_type bits_kmerpos; /* sum of bits_code, bits_seqnum, bits_endpos */
  GtUword first_seqnum;
  GtCodetype mask_seqnum, mask_endpos;
  int shift_code, shift_seqnum, shift_endpos;
} GtKmerPosListEncodeInfo;

typedef struct
{
  GtDiagbandseedKmerPos *spaceGtDiagbandseedKmerPos;
  GtUword *spaceGtUword;
  GtUword allocated, nextfree, longest_code_run;
  const GtKmerPosListEncodeInfo *encode_info;
} GtKmerPosList;

GT_DECLAREBufferedfiletype(GtDiagbandseedKmerPos);
GT_DECLAREREADFUNCTION(GtDiagbandseedKmerPos);

typedef struct
{ /* 4 + 4 + 4 + 4 bytes */
  GtDiagbandseedSeqnum bseqnum, /*  2nd important sort criterion */
                       aseqnum; /* most important sort criterion */
  GtDiagbandseedPosition apos,
                         bpos; /* 3rd important sort criterion */
} GtDiagbandseedSeedPair;

#define GT_DIAGBANDSEED_SETPOS(SP,APOS,BPOS)\
        (SP)->bpos = BPOS;\
        (SP)->apos = APOS

#define GT_DIAGBANDSEED_GETPOS_A(SP) (SP)->apos
#define GT_DIAGBANDSEED_GETPOS_B(SP) (SP)->bpos

#define GT_DIAGBANDSEED_DIAGONAL2BPOS(AMAXLEN,APOS,DIAG)\
        ((GtUword) (DIAG) + (GtUword) (APOS) - (AMAXLEN))

#define GT_DIAGBANDSEED_CONV_B(APOS,BPOS)\
        (seedpairlist->maxmat_compute \
          ? GT_DIAGBANDSEED_DIAGONAL2BPOS(seedpairlist->amaxlen,APOS,BPOS) \
          : (BPOS))

struct GtDiagbandseedInfo
{
  const GtEncseq *aencseq,
                 *bencseq;
  const GtDiagbandseedExtendParams *extp;
  const GtRange *seedpairdistance;
  const GtStr *chainarguments,
              *diagband_statistics_arg;
  GtUword maxfreq,
          memlimit,
          maxmat;
  unsigned int spacedseedweight,
               seedlength;
  GtSpacedSeedSpec *spaced_seed_spec;
  GtDiagbandseedBaseListType splt,
                             kmplt;
  bool norev,
       nofwd,
       verify,
       verbose,
       debug_kmer,
       debug_seedpair,
       use_kmerfile,
       trimstat_on;
};

struct GtDiagbandseedExtendParams
{
  GtUword logdiagbandwidth,
          mincoverage,
          maxalignedlendifference,
          history_size,
          perc_mat_history,
          sensitivity,
          userdefinedleastlength,
          errorpercentage;
  double evalue_threshold;
  GtXdropscore xdropbelowscore;
  GtExtendCharAccess a_extend_char_access,
                     b_extend_char_access;
  const GtSeedExtendDisplayFlag *out_display_flag;
  double matchscore_bias;
  GtUword use_apos;
  GtAniAccumulate *ani_accumulate;
  bool extendgreedy,
       extendxdrop,
       weakends,
       benchmark,
       always_polished_ends,
       verify_alignment,
       only_selected_seqpairs,
       cam_generic;
};

typedef struct
{
  GtKmerPosList *kmerpos_list_ref;
  GtDiagbandseedSeqnum current_seqnum;
  GtDiagbandseedPosition current_endpos;
  const GtEncseq *encseq;
  GtSpecialrangeiterator *sri;
  GtRange *specialrange;
  GtUword last_specialpos,
          prev_separator,
          next_separator;
  unsigned int spacedseedweight,
               seedlength;
  const GtSpacedSeedSpec *spaced_seed_spec;
  GtReadmode readmode;
} GtDiagbandseedProcKmerInfo;

typedef struct
{
  GtUword b_firstseq, b_numsequences;
  GtBitsequence *b_bitsequence;
} GtSegmentRejectInfo;

typedef bool (*GtSegmentRejectFunc)(const GtSegmentRejectInfo *,GtUword);

static GtSegmentRejectInfo *gt_segment_reject_info_new(GtUword b_firstseq,
                                                       GtUword b_numsequences)
{
  GtSegmentRejectInfo *segment_reject_info
    = gt_malloc(sizeof *segment_reject_info);

  gt_assert(b_numsequences > 0);
  GT_INITBITTAB(segment_reject_info->b_bitsequence,b_numsequences);
  segment_reject_info->b_firstseq = b_firstseq;
  segment_reject_info->b_numsequences = b_numsequences;
  return segment_reject_info;
}

static void gt_segment_reject_info_delete(GtSegmentRejectInfo
                                            *segment_reject_info)
{
  if (segment_reject_info != NULL)
  {
    gt_free(segment_reject_info->b_bitsequence);
    gt_free(segment_reject_info);
  }
}

static void gt_segment_reject_register_match(GtSegmentRejectInfo
                                                *segment_reject_info,
                                             GtUword bseqnum)
{
  GtUword idx;

  gt_assert(segment_reject_info !=NULL);
  idx = bseqnum - segment_reject_info->b_firstseq;
  gt_assert(bseqnum >= segment_reject_info->b_firstseq &&
            idx < segment_reject_info->b_numsequences &&
            !GT_ISIBITSET(segment_reject_info->b_bitsequence,idx));
  GT_SETIBIT(segment_reject_info->b_bitsequence,idx);
}

static bool gt_segment_reject_check(const GtSegmentRejectInfo
                                      *segment_reject_info,
                                       GtUword bseqnum)
{
  const GtUword idx = bseqnum - segment_reject_info->b_firstseq;

  gt_assert(bseqnum >= segment_reject_info->b_firstseq &&
            idx < segment_reject_info->b_numsequences);
  return GT_ISIBITSET(segment_reject_info->b_bitsequence,idx) ? true : false;
}

/* * * * * CONSTRUCTORS AND DESTRUCTORS * * * * */

GtDiagbandseedInfo *gt_diagbandseed_info_new(const GtEncseq *aencseq,
                                             const GtEncseq *bencseq,
                                             GtUword maxfreq,
                                             GtUword memlimit,
                                             unsigned int spacedseedweight,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             const GtRange *seedpairdistance,
                                             GtDiagbandseedBaseListType splt,
                                             GtDiagbandseedBaseListType kmplt,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             bool trimstat_on,
                                             GtUword maxmat,
                                             const GtStr *chainarguments,
                                             const GtStr
                                               *diagband_statistics_arg,
                                             const GtDiagbandseedExtendParams
                                               *extp)
{
  GtDiagbandseedInfo *info = gt_malloc(sizeof *info);

  info->aencseq = aencseq;
  info->bencseq = bencseq;
  info->maxfreq = maxfreq;
  info->memlimit = memlimit;
  if (spacedseedweight > 0)
  {
    info->spacedseedweight = spacedseedweight;
    gt_assert(spacedseedweight < seedlength);
    info->spaced_seed_spec
      = gt_spaced_seed_spec_new_from_ws((int) info->spacedseedweight,
                                        (int) seedlength);
  } else
  {
    info->spacedseedweight = seedlength;
    info->spaced_seed_spec = NULL;
  }
  info->seedlength = seedlength;
  info->norev = norev;
  info->nofwd = nofwd;
  info->seedpairdistance = seedpairdistance;
  info->splt = splt;
  info->kmplt = kmplt;
  info->verify = verify;
  info->verbose = verbose;
  info->debug_kmer = debug_kmer;
  info->debug_seedpair = debug_seedpair;
  info->use_kmerfile = use_kmerfile;
  info->trimstat_on = trimstat_on;
  info->maxmat = maxmat;
  info->chainarguments = chainarguments;
  info->diagband_statistics_arg = diagband_statistics_arg;
  info->extp = extp;
  return info;
}

void gt_diagbandseed_info_delete(GtDiagbandseedInfo *info)
{
  if (info != NULL) {
    gt_spaced_seed_spec_delete(info->spaced_seed_spec);
    gt_free(info);
  }
}

GtDiagbandseedExtendParams *gt_diagbandseed_extend_params_new(
                                GtUword userdefinedleastlength,
                                GtUword errorpercentage,
                                double evalue_threshold,
                                GtUword logdiagbandwidth,
                                GtUword mincoverage,
                                const GtSeedExtendDisplayFlag *out_display_flag,
                                GtUword use_apos,
                                GtXdropscore xdropbelowscore,
                                bool extendgreedy,
                                bool extendxdrop,
                                GtUword maxalignedlendifference,
                                GtUword history_size,
                                GtUword perc_mat_history,
                                GtExtendCharAccess a_extend_char_access,
                                GtExtendCharAccess b_extend_char_access,
                                bool cam_generic,
                                GtUword sensitivity,
                                double matchscore_bias,
                                bool weakends,
                                bool benchmark,
                                bool always_polished_ends,
                                bool verify_alignment,
                                bool only_selected_seqpairs,
                                GtAniAccumulate *ani_accumulate)
{
  GtDiagbandseedExtendParams *extp = gt_malloc(sizeof *extp);
  extp->userdefinedleastlength = userdefinedleastlength;
  extp->errorpercentage = errorpercentage;
  extp->evalue_threshold = evalue_threshold;
  extp->logdiagbandwidth = logdiagbandwidth;
  extp->mincoverage = mincoverage;
  extp->out_display_flag = out_display_flag;
  extp->use_apos = use_apos;
  extp->xdropbelowscore = xdropbelowscore;
  extp->extendgreedy = extendgreedy;
  extp->extendxdrop = extendxdrop;
  extp->maxalignedlendifference = maxalignedlendifference;
  extp->history_size = history_size;
  extp->perc_mat_history = perc_mat_history;
  extp->a_extend_char_access = a_extend_char_access;
  extp->b_extend_char_access = b_extend_char_access;
  extp->cam_generic = cam_generic;
  extp->sensitivity = sensitivity;
  extp->matchscore_bias = matchscore_bias;
  extp->weakends = weakends;
  extp->benchmark = benchmark;
  extp->always_polished_ends = always_polished_ends;
  extp->verify_alignment = verify_alignment;
  extp->only_selected_seqpairs = only_selected_seqpairs;
  extp->ani_accumulate = ani_accumulate;
  return extp;
}

void gt_diagbandseed_extend_params_delete(GtDiagbandseedExtendParams *extp)
{
  if (extp != NULL) {
    gt_free(extp);
  }
}

/* * * * * K-MER LIST CREATION * * * * */

/* Estimate the number of k-mers in the given encseq. */
static GtUword gt_seed_extend_numofkmers(const GtEncseq *encseq,
                                         unsigned int seedlength,
                                         GtUword seqrange_start,
                                         GtUword seqrange_end)
{
  GtUword lastpos, numofpos, subtract, ratioofspecial;
  const GtUword totalnumofspecial = gt_encseq_specialcharacters(encseq),
                totalnumofpos = gt_encseq_total_length(encseq),
                firstpos = gt_encseq_seqstartpos(encseq, seqrange_start),
                numofseq = seqrange_end - seqrange_start + 1;

  lastpos = (seqrange_end + 1 == gt_encseq_num_of_sequences(encseq)
             ? totalnumofpos
             : gt_encseq_seqstartpos(encseq, seqrange_end + 1) - 1);
  gt_assert(lastpos >= firstpos);
  numofpos = lastpos - firstpos;

  subtract = GT_MIN(seedlength - 1, gt_encseq_min_seq_length(encseq)) + 1;
  gt_assert(numofpos + 1 >= numofseq * subtract);
  ratioofspecial = GT_MIN(totalnumofspecial * numofpos / totalnumofpos,
                          numofpos);
  gt_assert(numofpos >= ratioofspecial);
  return numofpos - GT_MAX(numofseq * subtract - 1, ratioofspecial);
}

static GtKmerPosList *gt_kmerpos_list_new(GtUword init_size,
                                          const GtKmerPosListEncodeInfo
                                            *encode_info)
{
  GtKmerPosList *kmerpos_list = gt_malloc(sizeof *kmerpos_list);

  gt_assert(init_size > 0);
  kmerpos_list->nextfree = 0;
  kmerpos_list->allocated = init_size;
  if (encode_info != NULL)
  {
    kmerpos_list->spaceGtDiagbandseedKmerPos = NULL;
    kmerpos_list->spaceGtUword
      = gt_malloc(sizeof *kmerpos_list->spaceGtUword * init_size);
  } else
  {
    kmerpos_list->spaceGtUword = NULL;
    kmerpos_list->spaceGtDiagbandseedKmerPos
      = gt_malloc(sizeof *kmerpos_list->spaceGtDiagbandseedKmerPos *
                  init_size);
  }
  kmerpos_list->encode_info = encode_info;
  return kmerpos_list;
}

static void gt_kmerpos_list_delete(GtKmerPosList *kmerpos_list)
{
  if (kmerpos_list != NULL)
  {
    gt_free(kmerpos_list->spaceGtDiagbandseedKmerPos);
    gt_free(kmerpos_list->spaceGtUword);
    gt_free(kmerpos_list);
  }
}

static GtCodetype gt_kmerpos_entry_code(
                      const GtKmerPosListEncodeInfo *encode_info,
                      GtUword value)
{
  return value >> encode_info->shift_code;
}

static void gt_kmerpos_entry_decode(GtDiagbandseedKmerPos *dec,
                      const GtKmerPosListEncodeInfo *encode_info,
                      GtUword value)
{
  gt_assert(dec != NULL && encode_info != NULL);
  dec->code = gt_kmerpos_entry_code(encode_info,value);
  dec->seqnum = encode_info->first_seqnum +
                ((value >> encode_info->shift_seqnum) &
                encode_info->mask_seqnum);
  dec->endpos = (value >> encode_info->shift_endpos) & encode_info->mask_endpos;
}

static void gt_kmerpos_list_add(GtKmerPosList *kmerpos_list,
                                const GtDiagbandseedKmerPos *kmerpos_entry)
{
  gt_assert(kmerpos_list != NULL);
  if (kmerpos_list->nextfree >= kmerpos_list->allocated)
  {
    kmerpos_list->allocated = kmerpos_list->allocated * 1.2 + 256;
    if (kmerpos_list->encode_info != NULL)
    {
      kmerpos_list->spaceGtUword
        = gt_realloc(kmerpos_list->spaceGtUword,
                     sizeof *kmerpos_list->spaceGtUword *
                     kmerpos_list->allocated);
    } else
    {
      kmerpos_list->spaceGtDiagbandseedKmerPos
        = gt_realloc(kmerpos_list->spaceGtDiagbandseedKmerPos,
                     sizeof *kmerpos_list->spaceGtDiagbandseedKmerPos *
                     kmerpos_list->allocated);
    }
  }
  if (kmerpos_list->encode_info != NULL)
  {
    GtUword seqnum;
    const GtKmerPosListEncodeInfo *encode_info = kmerpos_list->encode_info;

    gt_assert(encode_info->first_seqnum <= kmerpos_entry->seqnum);
    seqnum = kmerpos_entry->seqnum - encode_info->first_seqnum;
    kmerpos_list->spaceGtUword[kmerpos_list->nextfree++]
      = (kmerpos_entry->code << encode_info->shift_code) |
        (seqnum << encode_info->shift_seqnum) |
        ((GtUword) kmerpos_entry->endpos << encode_info->shift_endpos);
  } else
  {
    kmerpos_list->spaceGtDiagbandseedKmerPos[kmerpos_list->nextfree++]
      = *kmerpos_entry;
  }
}

static size_t gt_kmerpos_list_elem_size(const GtKmerPosList *kmerpos_list)
{
  if (kmerpos_list->encode_info != NULL)
  {
    return sizeof *kmerpos_list->spaceGtUword;
  }
  return sizeof *kmerpos_list->spaceGtDiagbandseedKmerPos;
}

static void gt_kmerpos_list_reduce_size(GtKmerPosList *kmerpos_list)
{
  gt_assert(kmerpos_list != NULL &&
            kmerpos_list->nextfree <= kmerpos_list->allocated);
  if (kmerpos_list->nextfree < kmerpos_list->allocated)
  {
    if (kmerpos_list->encode_info != NULL)
    {
      kmerpos_list->spaceGtUword
        = gt_realloc(kmerpos_list->spaceGtUword,
                     kmerpos_list->nextfree *
                     sizeof *kmerpos_list->spaceGtUword);
    } else
    {
      kmerpos_list->spaceGtDiagbandseedKmerPos
        = gt_realloc(kmerpos_list->spaceGtDiagbandseedKmerPos,
                     kmerpos_list->nextfree *
                     sizeof *kmerpos_list->spaceGtDiagbandseedKmerPos);
    }
    kmerpos_list->allocated = kmerpos_list->nextfree;
  }
}

static void gt_kmerpos_list_sort(GtKmerPosList *kmerpos_list)
{

  if (kmerpos_list->encode_info != NULL)
  {
    gt_radixsort_inplace_ulong(kmerpos_list->spaceGtUword,
                               kmerpos_list->nextfree);
  } else
  {
    gt_radixsort_inplace_GtUwordPair((GtUwordPair *)
                                     kmerpos_list->spaceGtDiagbandseedKmerPos,
                                     kmerpos_list->nextfree);
  }
}

static void gt_kmerpos_list_show(FILE *stream,const GtKmerPosList *kmerpos_list)
{
  GtUword idx;

  for (idx = 0; idx < kmerpos_list->nextfree; idx++)
  {
    GtDiagbandseedKmerPos dec;

    if (kmerpos_list->encode_info != NULL)
    {
      gt_kmerpos_entry_decode(&dec,kmerpos_list->encode_info,
                              kmerpos_list->spaceGtUword[idx]);
    } else
    {
      dec = kmerpos_list->spaceGtDiagbandseedKmerPos[idx];
    }
    fprintf(stream, "# Kmer (" GT_LX ",%"PRIu32",%"PRIu32")\n",dec.code,
            dec.endpos, dec.seqnum);
  }
}

static GtUword gt_kmerpos_list_num_entries(const GtKmerPosList *kmerpos_list)
{
  gt_assert(kmerpos_list != NULL);
  return kmerpos_list->nextfree;
}

/* Returns the position of the next separator following specialrange.start.
   If the end of the encseq is reached, the position behind is returned. */
static GtUword gt_diagbandseed_update_separatorpos(GtRange *specialrange,
                                                   GtSpecialrangeiterator *sri,
                                                   const GtEncseq *encseq,
                                                   GtUword totallength,
                                                   GtReadmode readmode)
{
  gt_assert(sri != NULL && specialrange != NULL && encseq != NULL);
  do {
    GtUword idx;
    for (idx = specialrange->start; idx < specialrange->end; idx++) {
      if (gt_encseq_position_is_separator(encseq, idx, readmode)) {
        specialrange->start = idx + 1;
        return idx;
      }
    }
  } while (gt_specialrangeiterator_next(sri, specialrange));
  return totallength;
}

/* Add given code and its seqnum and position to a kmer list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword startpos,
                                            GtCodetype code)
{
  GtDiagbandseedProcKmerInfo *pkinfo;
  GtDiagbandseedKmerPos kmerpos_entry;

  gt_assert(prockmerinfo != NULL);
  pkinfo = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  /* check separator positions and determine next seqnum and endpos */
  if (firstinrange)
  {
    const GtUword current_endpos = startpos + pkinfo->seedlength - 1;

    while (current_endpos >= pkinfo->next_separator)
    {
      pkinfo->current_seqnum++;
      pkinfo->prev_separator = pkinfo->next_separator + 1;
      pkinfo->next_separator
        = gt_diagbandseed_update_separatorpos(pkinfo->specialrange,
                                              pkinfo->sri,
                                              pkinfo->encseq,
                                              pkinfo->last_specialpos,
                                              pkinfo->readmode);
      gt_assert(pkinfo->next_separator >= pkinfo->prev_separator);
    }
    gt_assert(current_endpos >= pkinfo->prev_separator &&
              startpos < pkinfo->next_separator);
    if (pkinfo->readmode == GT_READMODE_FORWARD)
    {
      pkinfo->current_endpos = (GtDiagbandseedPosition)
                               (current_endpos - pkinfo->prev_separator);
    } else
    {
      pkinfo->current_endpos = (GtDiagbandseedPosition)
                               (pkinfo->next_separator - 1 - startpos);
    }
  }

  /* save k-mer code */
  kmerpos_entry.code = (pkinfo->readmode == GT_READMODE_FORWARD)
                         ? code
                         : gt_kmercode_reverse(code, pkinfo->seedlength);
  if (pkinfo->spaced_seed_spec != NULL)
  {
    kmerpos_entry.code
      = gt_spaced_seed_extract_generic(pkinfo->spaced_seed_spec,
                                       kmerpos_entry.code);
  }
  /* save endpos and seqnum */
  gt_assert(pkinfo->current_endpos != UINT32_MAX);
  kmerpos_entry.endpos = pkinfo->current_endpos;
  pkinfo->current_endpos = (pkinfo->readmode == GT_READMODE_FORWARD)
                             ? pkinfo->current_endpos + 1
                             : pkinfo->current_endpos - 1;
  kmerpos_entry.seqnum = pkinfo->current_seqnum;
  gt_kmerpos_list_add(pkinfo->kmerpos_list_ref,&kmerpos_entry);
}

/* Uses GtKmercodeiterator for fetching the kmers. */
static void gt_diagbandseed_get_kmers_kciter(GtDiagbandseedProcKmerInfo *pkinfo)
{
  gt_assert(pkinfo != NULL);
  if (pkinfo->seedlength <= pkinfo->last_specialpos)
  {
    const GtKmercode *kmercode = NULL;
    bool firstinrange = true;
    GtKmercodeiterator *kc_iter = NULL;
    const GtUword maxpos = pkinfo->last_specialpos + 1 - pkinfo->seedlength;
    GtUword position = gt_encseq_seqstartpos(pkinfo->encseq,
                                             pkinfo->current_seqnum);

    kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq,
                                             pkinfo->readmode,
                                             pkinfo->seedlength,
                                             position);
    while (position < maxpos)
    {
      kmercode = gt_kmercodeiterator_encseq_next(kc_iter);
      if (!kmercode->definedspecialposition)
      {
        gt_diagbandseed_processkmercode((void *) pkinfo,
                                        firstinrange,
                                        position,
                                        kmercode->code);
        firstinrange = false;
      } else
      {
        firstinrange = true;
      }
      position++;
    }
    gt_kmercodeiterator_delete(kc_iter);
  }
}

static GtKmerPosListEncodeInfo *gt_kmerpos_encode_info_new(
                                       GtDiagbandseedBaseListType kmplt,
                                       const GtEncseq *encseq,
                                       GtUword spacedseedweight,
                                       const GtSequencePartsInfo *seqranges,
                                       GtUword idx)
{
  if (kmplt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    return NULL;
  } else
  {
    if (kmplt == GT_DIAGBANDSEED_BASE_LIST_ULONG ||
        kmplt == GT_DIAGBANDSEED_BASE_LIST_UNDEFINED)
    {
      GtBitcount_type bits_code, bits_seqnum, bits_endpos;
      const int allbits = sizeof (GtUword) * CHAR_BIT;
      GtUword seqrange_start = gt_sequence_parts_info_start_get(seqranges,idx),
              seqrange_end = gt_sequence_parts_info_end_get(seqranges,idx),
              max_endpos = gt_sequence_parts_info_max_length_get(seqranges,idx),
              numofsequences = seqrange_end - seqrange_start + 1;

      if (spacedseedweight >= 32)
      {
        bits_code = 64;
      } else
      {
        GtUword num_different_kmers
          = ceil(pow((double) gt_encseq_alphabetnumofchars(encseq),
                     (double) spacedseedweight));
        bits_code = (GtBitcount_type) gt_radixsort_bits(num_different_kmers);
      }
      bits_seqnum = (GtBitcount_type) gt_radixsort_bits(numofsequences);
      bits_endpos = (GtBitcount_type) gt_radixsort_bits(max_endpos);
      if (bits_code + bits_seqnum + bits_endpos <= 64)
      {
        GtKmerPosListEncodeInfo *encode_info = gt_malloc(sizeof *encode_info);

        encode_info->bits_kmerpos = bits_code + bits_seqnum + bits_endpos;
        encode_info->shift_code = allbits - bits_code;
        encode_info->shift_seqnum = encode_info->shift_code - bits_seqnum;
        encode_info->first_seqnum = seqrange_start;
        encode_info->shift_endpos = encode_info->shift_seqnum - bits_endpos;
        encode_info->mask_seqnum = (((GtCodetype) 1) << bits_seqnum) - 1;
        encode_info->mask_endpos = (((GtCodetype) 1) << bits_endpos) - 1;
        return encode_info;
      }
      return NULL;
    } else
    {
      gt_assert(kmplt == GT_DIAGBANDSEED_BASE_LIST_BYTESTRING);
      return NULL;
    }
  }
}

static void gt_kmerpos_encode_info_delete(GtKmerPosListEncodeInfo *encode_info)
{
  if (encode_info != NULL)
  {
    gt_free(encode_info);
  }
}

typedef GtUword GtLongestCodeRunType;

static GtUword gt_diagbandseed_longest_code_run(const GtKmerPosList
                                                   *kmerpos_list)
{
  GtUword idx;
  GtLongestCodeRunType longest_code_run = 1, current_code_run = 1;
  GtCodetype previouscode;

  gt_assert(kmerpos_list != NULL);
  if (kmerpos_list->encode_info != NULL)
  {
    gt_assert(kmerpos_list->spaceGtUword != NULL);
    previouscode = gt_kmerpos_entry_code(kmerpos_list->encode_info,
                                         kmerpos_list->spaceGtUword[0]);
    for (idx = 1; idx < kmerpos_list->nextfree; idx++)
    {
      const GtCodetype currentcode
        = gt_kmerpos_entry_code(kmerpos_list->encode_info,
                                kmerpos_list->spaceGtUword[idx]);

      if (previouscode == currentcode)
      {
        current_code_run++;
      } else
      {
        if (current_code_run > longest_code_run)
        {
          longest_code_run = current_code_run;
        }
        current_code_run = 1;
        previouscode = currentcode;
      }
    }
  } else
  {
    gt_assert(kmerpos_list->spaceGtDiagbandseedKmerPos != NULL);
    previouscode = kmerpos_list->spaceGtDiagbandseedKmerPos[0].code;
    for (idx = 1; idx < kmerpos_list->nextfree; idx++)
    {
      const GtCodetype currentcode
        = kmerpos_list->spaceGtDiagbandseedKmerPos[idx].code;

      if (previouscode == currentcode)
      {
        current_code_run++;
      } else
      {
        if (current_code_run > longest_code_run)
        {
          longest_code_run = current_code_run;
        }
        current_code_run = 1;
        previouscode = currentcode;
      }
    }
  }
  if (current_code_run > longest_code_run)
  {
    longest_code_run = current_code_run;
  }
  return longest_code_run;
}

/* Return a sorted list of k-mers of given seedlength from specified encseq.
 * Only sequences in seqrange will be taken into account.
 * The caller is responsible for freeing the result. */
static GtKmerPosList *gt_diagbandseed_get_kmers(
                                   const GtEncseq *encseq,
                                   unsigned int spacedseedweight,
                                   unsigned int seedlength,
                                   const GtSpacedSeedSpec *spaced_seed_spec,
                                   GtReadmode readmode,
                                   GtUword seqrange_start,
                                   GtUword seqrange_end,
                                   const GtKmerPosListEncodeInfo *encode_info,
                                   bool debug_kmer,
                                   bool verbose,
                                   GtUword known_size,
                                   FILE *stream)
{
  GtKmerPosList *kmerpos_list;
  GtDiagbandseedProcKmerInfo pkinfo;
  GtRange specialrange;
  GtTimer *timer = NULL;
  GtUword kmerpos_list_len, totallength;

  gt_assert(encseq != NULL);
  totallength = gt_encseq_total_length(encseq);
  if (known_size > 0)
  {
    kmerpos_list_len = known_size;
  } else
  {
    kmerpos_list_len = gt_seed_extend_numofkmers(encseq, seedlength,
                                                 seqrange_start, seqrange_end);
    gt_assert(kmerpos_list_len > 0);
  }
  kmerpos_list = gt_kmerpos_list_new(kmerpos_list_len,encode_info);
  if (verbose) {
    GtBitcount_type bits_kmerpos;

    if (encode_info != NULL)
    {
      bits_kmerpos = encode_info->bits_kmerpos;
    } else
    {
      bits_kmerpos = sizeof (GtDiagbandseedKmerPos) * CHAR_BIT;
    }
    fprintf(stream, "# start fetching %u-mers (%hu bits/" GT_WU " bytes each, "
                    "expect " GT_WU ", allocate %.0f MB) ...\n",
            seedlength,
            bits_kmerpos,
            (GtUword) gt_kmerpos_list_elem_size(kmerpos_list),
            kmerpos_list_len,
            GT_MEGABYTES(kmerpos_list_len *
                         gt_kmerpos_list_elem_size(kmerpos_list)));
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  pkinfo.kmerpos_list_ref = kmerpos_list;
  pkinfo.current_seqnum = seqrange_start;
  pkinfo.current_endpos = 0;
  pkinfo.encseq = encseq;
  pkinfo.spacedseedweight = spacedseedweight;
  pkinfo.seedlength = seedlength;
  pkinfo.spaced_seed_spec = spaced_seed_spec;
  pkinfo.readmode = readmode;
  if (seqrange_end + 1 == gt_encseq_num_of_sequences(encseq))
  {
    pkinfo.last_specialpos = totallength;
  } else
  {
    /* start position of following sequence, minus separator position */
    pkinfo.last_specialpos
      = gt_encseq_seqstartpos(encseq, seqrange_end + 1) - 1;
  }
  pkinfo.prev_separator = gt_encseq_seqstartpos(encseq, seqrange_start);
  if (gt_encseq_has_specialranges(encseq)) {
    pkinfo.sri = gt_specialrangeiterator_new(encseq, true);
    while (gt_specialrangeiterator_next(pkinfo.sri, &specialrange) &&
           specialrange.end < pkinfo.prev_separator)
       /* Nothing */;
    specialrange.start = pkinfo.prev_separator;
    pkinfo.specialrange = &specialrange;
    pkinfo.next_separator
      = gt_diagbandseed_update_separatorpos(pkinfo.specialrange,
                                            pkinfo.sri,
                                            pkinfo.encseq,
                                            totallength,
                                            pkinfo.readmode);
  } else {
    pkinfo.sri = NULL;
    pkinfo.specialrange = NULL;
    pkinfo.next_separator = pkinfo.last_specialpos;
  }

  if (gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    /* Use fast access to encseq, requires 2bit-enc and absence of wildcards. */
    gt_getencseqkmers_twobitencoding_slice(encseq,
                                           readmode,
                                           seedlength,
                                           seedlength,
                                           false,
                                           gt_diagbandseed_processkmercode,
                                           (void *) &pkinfo,
                                           NULL,
                                           NULL,
                                           pkinfo.prev_separator,
                                           pkinfo.last_specialpos);
  } else {
    /* Use GtKmercodeiterator for encseq access */
    gt_diagbandseed_get_kmers_kciter(&pkinfo);
  }
  if (gt_encseq_has_specialranges(encseq)) {
    gt_specialrangeiterator_delete(pkinfo.sri);
  }
  /* reduce size of array to number of entries */
  gt_kmerpos_list_reduce_size(kmerpos_list);
  if (debug_kmer)
  {
    gt_kmerpos_list_show(stream,kmerpos_list);
  }
  if (verbose)
  {
    fprintf(stream, "# ... collected " GT_WU " %u-mers ",
            gt_kmerpos_list_num_entries(kmerpos_list),seedlength);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_start(timer);
  }
  gt_kmerpos_list_sort(kmerpos_list);
  kmerpos_list->longest_code_run
    = gt_diagbandseed_longest_code_run(kmerpos_list);
  if (verbose) {
    fprintf(stream, "# ... sorted " GT_WU " %u-mers ",
            gt_kmerpos_list_num_entries(kmerpos_list),seedlength);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }

  return kmerpos_list;
}

/* * * * * SEEDPAIR LIST CREATION * * * * */

typedef struct {
  /* common for list/file based iterator */
  GtKmerPosList section;
  bool at_list_end;
  /* for list based iterator */
  const GtKmerPosList *original;
  const GtDiagbandseedKmerPos *listend_struct;
  GtDiagbandseedKmerPos *listptr_struct;
  const GtUword *listend_uword;
  GtUword *listptr_uword;
  /* for file based iterator */
  GtBufferedfile_GtDiagbandseedKmerPos kmerstream_struct;
  GtBufferedfile_GtUword kmerstream_uword;
  GtUword buffer_uword;
} GtDiagbandseedKmerIterator;

static void gt_diagbandseed_kmer_iter_reset(GtDiagbandseedKmerIterator *ki)
{
  gt_assert(ki != NULL);
  ki->at_list_end = false;
  if (ki->original != NULL) /* list based */
  {
    if (ki->section.encode_info != NULL)
    {
      ki->listptr_uword = ki->original->spaceGtUword;
      ki->section.spaceGtUword = ki->original->spaceGtUword;
      ki->listptr_struct = ki->section.spaceGtDiagbandseedKmerPos;
    } else
    {
      ki->listptr_struct = ki->original->spaceGtDiagbandseedKmerPos;
      ki->section.spaceGtDiagbandseedKmerPos
        = ki->original->spaceGtDiagbandseedKmerPos;
    }
    if (gt_kmerpos_list_num_entries(ki->original) == 0)
    {
      ki->at_list_end = true;
    }
  } else /* file based */
  {
    if (ki->section.encode_info != NULL)
    {
      ki->kmerstream_uword.nextread = ki->kmerstream_uword.nextfree = 0;
      fseek(ki->kmerstream_uword.fp,sizeof (GtLongestCodeRunType),SEEK_SET);
      if (gt_readnextfromstream_GtUword(&ki->buffer_uword,
                                        &ki->kmerstream_uword) != 1)
      {
        ki->at_list_end = true;
      }
    } else
    {
      ki->kmerstream_struct.nextread = ki->kmerstream_struct.nextfree = 0;
      fseek(ki->kmerstream_struct.fp,sizeof (GtLongestCodeRunType),SEEK_SET);
      if (gt_readnextfromstream_GtDiagbandseedKmerPos(ki->listptr_struct,
                                                      &ki->kmerstream_struct)
                                                      != 1)
      {
        ki->at_list_end = true;
      }
    }
  }
}

static GtDiagbandseedKmerIterator *gt_diagbandseed_kmer_iter_new_list(
                               const GtKmerPosList *original)
{
  GtDiagbandseedKmerIterator *ki = gt_malloc(sizeof *ki);

  gt_assert(original != NULL);
  ki->original = original;
  if (original->encode_info != NULL)
  {
    gt_assert(original->spaceGtUword != NULL);
    ki->listend_uword = original->spaceGtUword + original->nextfree;
    ki->section.spaceGtDiagbandseedKmerPos
      = gt_malloc(sizeof *ki->section.spaceGtDiagbandseedKmerPos
                  * original->longest_code_run);
    ki->listend_struct = ki->section.spaceGtDiagbandseedKmerPos +
                         original->longest_code_run;
  } else
  {
    gt_assert(original->spaceGtDiagbandseedKmerPos != NULL);
    ki->listend_uword = NULL;
    ki->listend_struct = original->spaceGtDiagbandseedKmerPos +
                         original->nextfree;
  }
  ki->section.encode_info = original->encode_info;
  gt_diagbandseed_kmer_iter_reset(ki);
  return ki;
}

static GtDiagbandseedKmerIterator *gt_diagbandseed_kmer_iter_new_file(FILE *fp,
                                   const GtKmerPosListEncodeInfo *encode_info)
{
  GtDiagbandseedKmerIterator *ki = gt_malloc(sizeof *ki);
  GtLongestCodeRunType longest_code_run;

  ki->original = NULL;
  ki->listend_uword = ki->listptr_uword = NULL;
  gt_assert(fp != NULL);
  gt_xfread(&longest_code_run,sizeof longest_code_run,1,fp);
  if (encode_info != NULL)
  {
    ki->kmerstream_uword.fp = fp;
    ki->kmerstream_uword.bufferedfilespace
      = gt_malloc(GT_FILEBUFFERSIZE *
                  sizeof *ki->kmerstream_uword.bufferedfilespace);
    ki->section.spaceGtUword = NULL;
    ki->section.allocated = longest_code_run;
  } else
  {
    ki->kmerstream_struct.fp = fp;
    ki->kmerstream_struct.bufferedfilespace
      = gt_malloc(GT_FILEBUFFERSIZE *
                  sizeof *ki->kmerstream_struct.bufferedfilespace);
    ki->section.allocated = longest_code_run + 1; /* one larger for reading
                                                     first of next run */
  }
  ki->section.spaceGtDiagbandseedKmerPos
    = gt_malloc(sizeof *ki->section.spaceGtDiagbandseedKmerPos *
                ki->section.allocated);
  ki->listptr_struct = ki->section.spaceGtDiagbandseedKmerPos;
  ki->listend_struct = ki->section.spaceGtDiagbandseedKmerPos +
                       ki->section.allocated;
  ki->section.encode_info = encode_info;
  gt_diagbandseed_kmer_iter_reset(ki);
  return ki;
}

static void gt_diagbandseed_kmer_iter_delete(GtDiagbandseedKmerIterator *ki)
{
  if (ki != NULL) {
    if (ki->original == NULL)
    { /* file based */
      if (ki->section.encode_info != NULL)
      {
        gt_free(ki->kmerstream_uword.bufferedfilespace);
        gt_fa_fclose(ki->kmerstream_uword.fp);
        gt_free(ki->section.spaceGtUword);
      } else
      {
        gt_free(ki->kmerstream_struct.bufferedfilespace);
        gt_fa_fclose(ki->kmerstream_struct.fp);
      }
      gt_free(ki->section.spaceGtDiagbandseedKmerPos);
      ki->section.allocated = 0;
    } else
    {
      if (ki->section.encode_info != NULL)
      {
        gt_free(ki->section.spaceGtDiagbandseedKmerPos);
      }
    }
    gt_free(ki);
  }
}

static const GtKmerPosList *gt_diagbandseed_kmer_iter_next(
                                              GtDiagbandseedKmerIterator *ki)
{
  GtCodetype code;

  if (ki->at_list_end)
  {
    return NULL;
  }
  /* reset section list */
  if (ki->original != NULL)
  { /* list based */
    if (ki->section.encode_info != NULL)
    {
      code = gt_kmerpos_entry_code(ki->section.encode_info,*ki->listptr_uword);
      ki->listptr_struct = ki->section.spaceGtDiagbandseedKmerPos;
      do
      {
        gt_assert(ki->listptr_struct < ki->listend_struct);
        gt_kmerpos_entry_decode(ki->listptr_struct++,
                                ki->section.encode_info,*ki->listptr_uword);
        ki->listptr_uword++;
      } while (ki->listptr_uword < ki->listend_uword &&
               code == gt_kmerpos_entry_code(ki->section.encode_info,
                                             *ki->listptr_uword));
      if (ki->listptr_uword >= ki->listend_uword)
      {
        ki->at_list_end = true;
      }
    } else
    {
      code = ki->listptr_struct->code;
      ki->section.spaceGtDiagbandseedKmerPos = ki->listptr_struct;
      /* add element to section until code differs */
      do
      {
        ki->listptr_struct++;
      } while (ki->listptr_struct < ki->listend_struct &&
               code == ki->listptr_struct->code);
      if (ki->listptr_struct >= ki->listend_struct)
      {
        ki->at_list_end = true;
      }
    }
  } else
  { /* file based */
    int rval;

    if (ki->section.encode_info != NULL)
    {
      ki->listptr_struct = ki->section.spaceGtDiagbandseedKmerPos;
      code = gt_kmerpos_entry_code(ki->section.encode_info,ki->buffer_uword);
      do
      {
        gt_assert(ki->listptr_struct < ki->listend_struct);
        gt_kmerpos_entry_decode(ki->listptr_struct++,
                                ki->section.encode_info,ki->buffer_uword);
        rval = gt_readnextfromstream_GtUword(&ki->buffer_uword,
                                             &ki->kmerstream_uword);
      } while (rval == 1 &&
               code == gt_kmerpos_entry_code(ki->section.encode_info,
                                             ki->buffer_uword));
    } else
    {
      if (ki->section.spaceGtDiagbandseedKmerPos < ki->listptr_struct)
      {
        ki->section.spaceGtDiagbandseedKmerPos[0] = *ki->listptr_struct;
        ki->listptr_struct = ki->section.spaceGtDiagbandseedKmerPos;
      }
      code = ki->listptr_struct->code;
      do
      {
        ki->listptr_struct++;
        gt_assert(ki->listptr_struct < ki->listend_struct);
        rval = gt_readnextfromstream_GtDiagbandseedKmerPos(ki->listptr_struct,
                                                           &ki->
                                                             kmerstream_struct);
      } while (rval == 1 && code == ki->listptr_struct->code);
    }
    if (rval != 1)
    {
      ki->at_list_end = true;
    }
  }
  ki->section.nextfree = (GtUword) (ki->listptr_struct -
                                    ki->section.spaceGtDiagbandseedKmerPos);
  return &ki->section;
}

/* Evaluate the results of the seed count histogram */
static GtUword gt_diagbandseed_processhistogram(GtUword *histogram,
                                                GtUword maxfreq,
                                                GtUword maxgram,
                                                GtUword memlimit,
                                                GtUword mem_used,
                                                bool alist_blist_id,
                                                size_t sizeofunit)
{
  /* calculate available memory, take 98% of memlimit */
  GtUword count = 0, frequency = 0, mem_avail = 0.98 * memlimit;

  if (mem_avail > mem_used) {
    mem_avail = (mem_avail - mem_used) / sizeofunit;
  } else {
    mem_avail = 0;
    maxfreq = 0;
  }

  /* there is enough free memory */
  if (mem_avail > 0) {
    /* count seeds until available memory reached */
    for (frequency = 1; frequency <= maxgram && count < mem_avail;
         frequency++) {
      count += histogram[frequency - 1];
    }
    if (count > mem_avail) {
      gt_assert(frequency >= 2);
      frequency -= 2;
      gt_assert(count >= histogram[frequency]);
      count -= histogram[frequency];
    } else if (frequency == maxgram + 1) {
      frequency = GT_UWORD_MAX;
    }
    maxfreq = GT_MIN(maxfreq, frequency);
  }

  /* determine minimum required memory for error message */
  if (maxfreq <= 1 && alist_blist_id) {
    count = (histogram[0] + histogram[1]) * sizeofunit;
    count = (count + mem_used) / 0.98;
  } else if (maxfreq == 0) {
    count = histogram[0] * sizeofunit;
    count = (count + mem_used) / 0.98;
  }
  histogram[maxgram] = count;
  return maxfreq;
}

const char *gt_diagbandseed_splt_comment(void)
{
  return "specify type of pairlist, possible values are struct, bytestring, "
         "and ulong";
}

const char *gt_diagbandseed_kmplt_comment(void)
{
  return "specify type of kmerpos, possible values are struct and ulong";
}

static const char *gt_base_list_arguments[]
  = {"struct","ulong","bytestring",""};

GtDiagbandseedBaseListType gt_diagbandseed_base_list_get(
                                   bool with_splt,
                                   const char *base_list_string,GtError *err)
{
  size_t idx;
  for (idx = 0; idx < sizeof gt_base_list_arguments/
                      sizeof gt_base_list_arguments[0]; idx++)
  {
    if (strcmp(base_list_string,gt_base_list_arguments[idx]) == 0)
    {
      return (GtDiagbandseedBaseListType) idx;
    }
  }
  gt_error_set(err,"illegal parameter for option -%s: %s",
                    with_splt ? "splt" : "kmplt",
                    with_splt ? gt_diagbandseed_splt_comment()
                              : gt_diagbandseed_kmplt_comment());
  return -1;
}

const int idx_aseqnum = 0, idx_bseqnum = 1, idx_bpos = 2, idx_apos = 3;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

typedef struct
{
  GtArrayGtDiagbandseedSeedPair *mlist_struct;
  GtArrayGtUword *mlist_ulong;
  GtArrayuint8_t *mlist_bytestring;
  GtDiagbandseedBaseListType splt;
  GtUword mask_tab[4], transfer_mask,
          aseqrange_start, bseqrange_start,
          aseqrange_end, bseqrange_end,
          aseqrange_max_length, bseqrange_max_length;
  GtBitcount_type bits_seedpair,
                  bits_values[4],
                  bits_units[2],
                  bits_left_adjust[4];
  size_t bytes_seedpair;
  int shift_tab[4],
      transfer_shift,
      bits_unused_in2GtUwords;
  bool maxmat_compute, maxmat_show;
  GtUword amaxlen;
} GtSeedpairlist;

#define GT_DIAGBANDSEED_ENCODE_SEQNUMS(ASEQNUM,BSEQNUM)\
             (((GtUword) (ASEQNUM) << seedpairlist->shift_tab[idx_aseqnum]) | \
              ((GtUword) (BSEQNUM) << seedpairlist->shift_tab[idx_bseqnum]))

#define GT_DIAGBANDSEED_ENCODE_POSITIONS(BPOS,APOS)\
             (((GtUword) (BPOS) << seedpairlist->shift_tab[idx_bpos]) | \
              ((GtUword) (APOS) << seedpairlist->shift_tab[idx_apos]))

#define GT_DIAGBANDSEED_ENCODE_SEEDPAIR(ASEQNUM,BSEQNUM,BPOS,APOS)\
        (GT_DIAGBANDSEED_ENCODE_SEQNUMS(ASEQNUM,BSEQNUM) | \
         GT_DIAGBANDSEED_ENCODE_POSITIONS(BPOS,APOS))

static bool gt_diagbandseed_derive_maxmat_show(GtUword maxmat)
{
  return maxmat == 1 ? true : false;
}

static GtSeedpairlist *gt_seedpairlist_new(GtDiagbandseedBaseListType splt,
                        const GtSequencePartsInfo *aseqranges,
                        GtUword aidx,
                        const GtSequencePartsInfo *bseqranges,
                        GtUword bidx,
                        GtUword maxmat,
                        GtUword amaxlen)
{
  GtSeedpairlist *seedpairlist = gt_malloc(sizeof *seedpairlist);
  int idx;
  const int allbits = sizeof (GtWord) * CHAR_BIT;
  const GtUword
    anumofseq = gt_sequence_parts_info_numofsequences_get(aseqranges,aidx),
    bnumofseq = gt_sequence_parts_info_numofsequences_get(bseqranges,bidx);

  gt_assert(maxmat <= 2);
  seedpairlist->maxmat_show = gt_diagbandseed_derive_maxmat_show(maxmat);
  seedpairlist->maxmat_compute = maxmat > 0 ? true : false;
  seedpairlist->amaxlen = amaxlen;
  seedpairlist->aseqrange_max_length
    = gt_sequence_parts_info_max_length_get(aseqranges,aidx);
  seedpairlist->bseqrange_max_length
    = gt_sequence_parts_info_max_length_get(bseqranges,bidx);
  seedpairlist->bits_values[idx_aseqnum]
    = (GtBitcount_type) gt_radixsort_bits(anumofseq);
  seedpairlist->bits_values[idx_bseqnum]
    = (GtBitcount_type) gt_radixsort_bits(bnumofseq);
  if (seedpairlist->maxmat_compute)
  {
    seedpairlist->bits_values[idx_bpos]
      = (GtBitcount_type) gt_radixsort_bits(
                             seedpairlist->aseqrange_max_length +
                             seedpairlist->bseqrange_max_length + 1);
    /* in this case we store diagonal values in SeedPair and we have
        make sure that they fit into the space for a position */
    gt_assert(seedpairlist->bits_values[idx_bpos] <=
              sizeof (GtDiagbandseedPosition) * CHAR_BIT);
  } else
  {
    seedpairlist->bits_values[idx_bpos]
      = (GtBitcount_type) gt_radixsort_bits(
                                         seedpairlist->bseqrange_max_length);
  }
  seedpairlist->bits_values[idx_apos]
    = (GtBitcount_type) gt_radixsort_bits(seedpairlist->aseqrange_max_length);
  seedpairlist->bits_units[0] = seedpairlist->bits_values[idx_aseqnum] +
                                seedpairlist->bits_values[idx_bseqnum];
  seedpairlist->bits_units[1] = seedpairlist->bits_values[idx_apos] +
                                seedpairlist->bits_values[idx_bpos];
  seedpairlist->bits_left_adjust[idx_aseqnum]
    = allbits - seedpairlist->bits_values[idx_aseqnum];
  seedpairlist->bits_left_adjust[idx_bseqnum]
    = allbits - seedpairlist->bits_units[0];
  seedpairlist->bits_left_adjust[idx_bpos]
    = allbits - seedpairlist->bits_values[idx_bpos];
  seedpairlist->bits_left_adjust[idx_apos]
    = allbits - seedpairlist->bits_units[1];
  seedpairlist->transfer_mask
    = (((GtUword) 1) << seedpairlist->bits_left_adjust[idx_bseqnum]) - 1;
  gt_assert(seedpairlist->bits_values[idx_apos] > 0);
  gt_assert(seedpairlist->bits_values[idx_bpos] > 0);
  for (idx = 0, seedpairlist->bits_seedpair = 0; idx < 4; idx++)
  {
    seedpairlist->bits_seedpair += seedpairlist->bits_values[idx];
    seedpairlist->mask_tab[idx]
      = (((GtUword) 1) << seedpairlist->bits_values[idx]) - 1;
  }
  /* the following is only used for bits_seedpair > allbits */
  seedpairlist->transfer_shift = seedpairlist->bits_seedpair - allbits;
  seedpairlist->bits_unused_in2GtUwords = 2 * allbits -
                                          seedpairlist->bits_seedpair;
  seedpairlist->bytes_seedpair
    = gt_radixsort_bits2bytes(seedpairlist->bits_seedpair);
  if (seedpairlist->bytes_seedpair <= sizeof (GtUword))
  {
    int shift = seedpairlist->bits_seedpair;
    for (idx = 0; idx < 4; idx++)
    {
      shift -= seedpairlist->bits_values[idx];
      seedpairlist->shift_tab[idx] = shift;
    }
  }
  seedpairlist->aseqrange_start
    = gt_sequence_parts_info_start_get(aseqranges,aidx);
  seedpairlist->bseqrange_start
    = gt_sequence_parts_info_start_get(bseqranges,bidx);
  seedpairlist->aseqrange_end
    = gt_sequence_parts_info_end_get(aseqranges,aidx);
  seedpairlist->bseqrange_end
    = gt_sequence_parts_info_end_get(bseqranges,bidx);
  seedpairlist->mlist_struct = NULL;
  seedpairlist->mlist_ulong = NULL;
  seedpairlist->mlist_bytestring = NULL;
  if (splt == GT_DIAGBANDSEED_BASE_LIST_UNDEFINED)
  {
    if (seedpairlist->bytes_seedpair <= sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_BASE_LIST_ULONG;
    } else
    {
      splt = GT_DIAGBANDSEED_BASE_LIST_BYTESTRING;
    }
  }
  if (splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
  {
    if (seedpairlist->bytes_seedpair > sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_BASE_LIST_BYTESTRING;
    }
  } else
  {
    if (splt == GT_DIAGBANDSEED_BASE_LIST_BYTESTRING &&
        seedpairlist->bytes_seedpair <= sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_BASE_LIST_ULONG;
    }
  }
  if (splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
  {
    gt_assert(seedpairlist->bytes_seedpair <= sizeof (GtUword));
    seedpairlist->mlist_ulong = gt_malloc(sizeof *seedpairlist->mlist_ulong);
    GT_INITARRAY(seedpairlist->mlist_ulong, GtUword);
  } else
  {
    if (splt == GT_DIAGBANDSEED_BASE_LIST_BYTESTRING)
    {
      gt_assert(seedpairlist->bytes_seedpair > sizeof (GtUword));
      seedpairlist->mlist_bytestring
        = gt_malloc(sizeof *seedpairlist->mlist_bytestring);
      GT_INITARRAY(seedpairlist->mlist_bytestring, uint8_t);
    } else
    {
      seedpairlist->mlist_struct
        = gt_malloc(sizeof *seedpairlist->mlist_struct);
      GT_INITARRAY(seedpairlist->mlist_struct, GtDiagbandseedSeedPair);
    }
  }
  seedpairlist->splt = splt;
  return seedpairlist;
}

static void gt_seedpairlist_reset(GtSeedpairlist *seedpairlist)
{
  if (seedpairlist->mlist_ulong != NULL)
  {
    seedpairlist->mlist_ulong->nextfreeGtUword = 0;
  }
  if (seedpairlist->mlist_bytestring != NULL)
  {
    seedpairlist->mlist_bytestring->nextfreeuint8_t = 0;
  }
  if (seedpairlist->mlist_struct != NULL)
  {
    seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair = 0;
  }
}

static void gt_seedpairlist_show_bits(FILE *stream,
                                      const GtSeedpairlist *seedpairlist)
{
  fprintf(stream,"# splt=%s, bits_seedpair=%d, bytes_seedpair=%u with ",
                    gt_base_list_arguments[seedpairlist->splt],
                    (int) seedpairlist->bits_seedpair,
                    (int) seedpairlist->bytes_seedpair);
  fprintf(stream,"aseqnum=%hu bits, ",seedpairlist->bits_values[idx_aseqnum]);
  fprintf(stream,"bseqnum=%hu bits, ",seedpairlist->bits_values[idx_bseqnum]);
  fprintf(stream,"bpos=%hu bits, ",seedpairlist->bits_values[idx_bpos]);
  fprintf(stream,"apos=%hu bits\n",seedpairlist->bits_values[idx_apos]);
}

static size_t gt_seedpairlist_sizeofunit(const GtSeedpairlist *seedpairlist)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    return sizeof (GtDiagbandseedSeedPair);
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
    {
      return sizeof (GtUword);
    } else
    {
      return seedpairlist->bytes_seedpair;
    }
  }
}

static void gt_seedpairlist_init(GtSeedpairlist *seedpairlist,
                                 GtUword known_size)
{
  if (known_size > 0) {
    gt_assert(seedpairlist != NULL);
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
    {
      gt_assert(seedpairlist->mlist_struct != NULL);
      GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_struct,
                              GtDiagbandseedSeedPair, known_size);
    } else
    {
      if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
      {
        gt_assert(seedpairlist->mlist_ulong != NULL);
        GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_ulong,GtUword, known_size);
      } else
      {
        gt_assert(seedpairlist->mlist_bytestring != NULL);
        GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_bytestring,uint8_t,
                                seedpairlist->bytes_seedpair * known_size);
      }
    }
  }
}

static void gt_seedpairlist_delete(GtSeedpairlist *seedpairlist)
{
  if (seedpairlist != NULL)
  {
    if (seedpairlist->mlist_struct != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_struct, GtDiagbandseedSeedPair);
      gt_free(seedpairlist->mlist_struct);
    }
    if (seedpairlist->mlist_ulong != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_ulong, GtUword);
      gt_free(seedpairlist->mlist_ulong);
    }
    if (seedpairlist->mlist_bytestring != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_bytestring, uint8_t);
      gt_free(seedpairlist->mlist_bytestring);
    }
    gt_free(seedpairlist);
  }
}

static GtUword gt_seedpairlist_length(const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL);
  if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    gt_assert(seedpairlist->mlist_struct != NULL);
    return seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair;
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
    {
      gt_assert(seedpairlist->mlist_ulong != NULL);
      return seedpairlist->mlist_ulong->nextfreeGtUword;
    } else
    {
      gt_assert(seedpairlist->mlist_bytestring != NULL &&
                (seedpairlist->mlist_bytestring->nextfreeuint8_t %
                 seedpairlist->bytes_seedpair == 0));
      return seedpairlist->mlist_bytestring->nextfreeuint8_t/
             seedpairlist->bytes_seedpair;
    }
  }
  return 0;
}

static const GtDiagbandseedSeedPair *gt_seedpairlist_mlist_struct(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL &&
            seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT &&
            seedpairlist->mlist_struct != NULL);
  return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair;
}

static const GtUword *gt_seedpairlist_mlist_ulong(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL && seedpairlist->mlist_ulong != NULL);
  return seedpairlist->mlist_ulong->spaceGtUword;
}

const uint8_t *gt_seedpairlist_mlist_bytestring(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL && seedpairlist->mlist_bytestring != NULL);
  return seedpairlist->mlist_bytestring->spaceuint8_t;
}

static void gt_diagbandseed_one_GtUword2bytestring(uint8_t *bytestring,
                                                   GtUword bytestring_length,
                                                   GtUword value)
{
  int idx, shift;

  for (idx = 0, shift = (sizeof (GtUword) - 1) * CHAR_BIT;
       idx < bytestring_length; idx++, shift -= CHAR_BIT)
  {
    bytestring[idx] = (uint8_t) (value >> shift);
  }
}

static void gt_diagbandseed_encode_seedpair(uint8_t *bytestring,
                                            const GtSeedpairlist *seedpairlist,
                                            GtDiagbandseedSeqnum aseqnum,
                                            GtDiagbandseedSeqnum bseqnum,
                                            GtDiagbandseedPosition bpos,
                                            GtDiagbandseedPosition apos)
{
  GtUword value_seqnums, value_positions;

  value_seqnums
    = (((GtUword) aseqnum) << seedpairlist->bits_left_adjust[idx_aseqnum]) |
      (((GtUword) bseqnum) << seedpairlist->bits_left_adjust[idx_bseqnum]);
  value_positions
    = (((GtUword) bpos) << seedpairlist->bits_left_adjust[idx_bpos]) |
      (((GtUword) apos) << seedpairlist->bits_left_adjust[idx_apos]);
  value_seqnums |= (value_positions >> seedpairlist->bits_units[0]);
  gt_diagbandseed_one_GtUword2bytestring(bytestring,sizeof (GtUword),
                                         value_seqnums);
  value_positions <<= seedpairlist->bits_left_adjust[idx_bseqnum];
  gt_diagbandseed_one_GtUword2bytestring(bytestring + sizeof (GtUword),
                                         seedpairlist->bytes_seedpair
                                         - sizeof (GtUword),
                                         value_positions);
}

static void gt_seedpairlist_add(GtSeedpairlist *seedpairlist,
                                bool knownsize,
                                GtDiagbandseedSeqnum aseqnum,
                                GtDiagbandseedSeqnum bseqnum,
                                GtDiagbandseedPosition bpos,
                                GtDiagbandseedPosition apos)
{
  gt_assert(seedpairlist != NULL);
  gt_assert(aseqnum >= seedpairlist->aseqrange_start &&
            aseqnum <= seedpairlist->aseqrange_end &&
            bseqnum >= seedpairlist->bseqrange_start &&
            bseqnum <= seedpairlist->bseqrange_end &&
            apos < seedpairlist->aseqrange_max_length &&
            bpos < seedpairlist->bseqrange_max_length);
  if (seedpairlist->maxmat_compute)
  {
    bpos = GT_DIAGBANDSEED_DIAGONAL(seedpairlist->amaxlen,apos,bpos);
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    GtDiagbandseedSeedPair *seedpair = NULL;
    if (knownsize)
    {
      gt_assert(seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair <
                seedpairlist->mlist_struct->allocatedGtDiagbandseedSeedPair);
      seedpair = seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair +
                 seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair++;
    } else
    {
      GT_GETNEXTFREEINARRAY(seedpair,
                            seedpairlist->mlist_struct,
                            GtDiagbandseedSeedPair,
                            256 + 0.2 *
                            seedpairlist->
                               mlist_struct->allocatedGtDiagbandseedSeedPair);
    }
    seedpair->aseqnum = aseqnum;
    seedpair->bseqnum = bseqnum;
    GT_DIAGBANDSEED_SETPOS(seedpair,apos,bpos); /* set splt struct */
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
    {
      const GtUword encoding
        = GT_DIAGBANDSEED_ENCODE_SEEDPAIR(aseqnum -
                                          seedpairlist->aseqrange_start,
                                          bseqnum -
                                          seedpairlist->bseqrange_start,
                                          bpos,apos); /*set splt bytestring */
      if (knownsize)
      {
        gt_assert(seedpairlist->mlist_ulong->nextfreeGtUword <
                  seedpairlist->mlist_ulong->allocatedGtUword);
        seedpairlist->mlist_ulong->spaceGtUword[
                   seedpairlist->mlist_ulong->nextfreeGtUword++] = encoding;
      } else
      {
        GT_STOREINARRAY(seedpairlist->mlist_ulong,
                        GtUword,
                        256 + 0.2 * seedpairlist->mlist_ulong->allocatedGtUword,
                        encoding);
      }
    } else
    {
       uint8_t *bytestring;

       if (knownsize)
       {
         gt_assert(seedpairlist->mlist_bytestring->nextfreeuint8_t +
                   seedpairlist->bytes_seedpair <=
                   seedpairlist->mlist_bytestring->allocateduint8_t);
       } else
       {
         GT_CHECKARRAYSPACE_GENERIC(seedpairlist->mlist_bytestring,
                                    uint8_t,
                                    seedpairlist->bytes_seedpair,
                                    256 + 0.2 *
                                    seedpairlist->
                                    mlist_bytestring->allocateduint8_t);
       }
       bytestring = seedpairlist->mlist_bytestring->spaceuint8_t +
                    seedpairlist->mlist_bytestring->nextfreeuint8_t;
       seedpairlist->mlist_bytestring->nextfreeuint8_t
         += seedpairlist->bytes_seedpair;
       gt_assert(aseqnum >= seedpairlist->aseqrange_start);
       gt_assert(bseqnum >= seedpairlist->bseqrange_start);
       gt_diagbandseed_encode_seedpair(bytestring,
                                       seedpairlist,
                                       aseqnum - seedpairlist->aseqrange_start,
                                       bseqnum - seedpairlist->bseqrange_start,
                                       bpos,
                                       apos);
    }
  }
}

static void gt_diagbandseed_seedpairlist_sort(GtSeedpairlist *seedpairlist)
{
  GtUword mlistlen = gt_seedpairlist_length(seedpairlist);

  if (mlistlen > 0)
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
    {
      gt_assert(seedpairlist->mlist_struct != NULL);
      gt_radixsort_inplace_Gtuint64keyPair(
            (Gtuint64keyPair *) seedpairlist->mlist_struct
                                            ->spaceGtDiagbandseedSeedPair,
            mlistlen);
    } else
    {
      if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
      {
        gt_assert(seedpairlist->mlist_ulong != NULL);
        gt_radixsort_inplace_ulong(seedpairlist->mlist_ulong->spaceGtUword,
                                   mlistlen);
      } else
      {
        gt_assert(seedpairlist->mlist_bytestring != NULL);
        gt_radixsort_inplace_flba(seedpairlist->mlist_bytestring->spaceuint8_t,
                                  gt_seedpairlist_length(seedpairlist),
                                  seedpairlist->bytes_seedpair);
      }
    }
  }
}

static GtUword gt_seedpairlist_a_bseqnum_ulong(
                                          const GtSeedpairlist *seedpairlist,
                                          GtUword encoding)
{
  return encoding >> seedpairlist->shift_tab[idx_bseqnum];
}

static GtUword gt_seedpairlist_extract_ulong(const GtSeedpairlist *seedpairlist,
                                             GtUword encoding,
                                             int compidx)
{
  return (encoding >> seedpairlist->shift_tab[compidx]) &
         seedpairlist->mask_tab[compidx];
}

static GtUword gt_seedpairlist_extract_ulong_at(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx,
                                   int compidx)
{
  return gt_seedpairlist_extract_ulong(seedpairlist,
                                       seedpairlist->mlist_ulong->
                                                     spaceGtUword[spidx],
                                       compidx);
}

static GtUword gt_diagbandseed_bytestring2GtUword(const uint8_t *bytestring,
                                                  GtUword bytestring_length)
{
  int idx, shift;
  GtUword value = 0;

  for (idx = 0, shift = (sizeof (GtUword) - 1) * CHAR_BIT;
       idx < bytestring_length; idx++, shift -= CHAR_BIT)
  {
    value |= (((GtUword) bytestring[idx]) << shift);
  }
  return value;
}

static void gt_diagbandseed_decode_seedpair(GtDiagbandseedSeedPair *seedpair,
                                            const GtSeedpairlist *seedpairlist,
                                            GtUword offset)
{
  GtUword transfer, value_seqnums, value_positions;
  GtDiagbandseedPosition apos, bpos;
  const uint8_t *bytestring = seedpairlist->mlist_bytestring->spaceuint8_t +
                              offset;

  value_seqnums
    = gt_diagbandseed_bytestring2GtUword(bytestring,sizeof (GtUword));
  seedpair->aseqnum = (GtDiagbandseedSeqnum)
                      ((value_seqnums >>
                        seedpairlist->bits_left_adjust[idx_aseqnum]) &
                       seedpairlist->mask_tab[idx_aseqnum]);
  seedpair->bseqnum = (GtDiagbandseedSeqnum)
                      ((value_seqnums >>
                        seedpairlist->bits_left_adjust[idx_bseqnum]) &
                       seedpairlist->mask_tab[idx_bseqnum]);
  transfer = value_seqnums & seedpairlist->transfer_mask;
  value_positions
    = gt_diagbandseed_bytestring2GtUword(bytestring + sizeof (GtUword),
                                         seedpairlist->bytes_seedpair -
                                         sizeof (GtUword));
  value_positions >>= seedpairlist->bits_unused_in2GtUwords;
  value_positions |= (transfer << seedpairlist->transfer_shift);
  apos = (GtDiagbandseedPosition) (value_positions &
                                   seedpairlist->mask_tab[idx_apos]);
  bpos = (GtDiagbandseedPosition) (value_positions >>
                                   seedpairlist->bits_values[idx_apos]);
  GT_DIAGBANDSEED_SETPOS(seedpair,apos,bpos); /* for splt bytestring */
}

static void gt_seedpairlist_at(GtDiagbandseedSeedPair *seedpair,
                               const GtSeedpairlist *seedpairlist,
                               GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    *seedpair = seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx];
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
    {
      GtDiagbandseedPosition apos, bpos;

      seedpair->aseqnum
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_aseqnum);
      seedpair->bseqnum
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_bseqnum);
      apos = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_apos),
      bpos = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_bpos);
      GT_DIAGBANDSEED_SETPOS(seedpair,apos,bpos); /* set for ulong */
    } else
    {
      gt_diagbandseed_decode_seedpair(seedpair,seedpairlist,
                                      spidx * seedpairlist->bytes_seedpair);
    }
  }
}

static void gt_diagbandseed_show_seed(FILE *stream,
                                      const GtSeedpairlist *seedpairlist,
                                      GtUword spidx)
{
  GtDiagbandseedSeedPair seedpair;
  GtDiagbandseedPosition apos;

  gt_seedpairlist_at(&seedpair,seedpairlist,spidx);
  apos = GT_DIAGBANDSEED_GETPOS_A(&seedpair);
  fprintf(stream, "(%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" GT_WUS ")",
                  seedpair.aseqnum,
                  seedpair.bseqnum,
                  apos,
                  GT_DIAGBANDSEED_CONV_B(apos,
                                         GT_DIAGBANDSEED_GETPOS_B(&seedpair)));
}

#ifndef NDEBUG
static int gt_diagbandseed_seeds_compare(const GtSeedpairlist *seedpairlist,
                                         const GtUword current)
{
  GtDiagbandseedSeedPair p_seedpair, c_seedpair;
  GtDiagbandseedPosition p_apos, c_apos, p_bpos, c_bpos;

  gt_assert(current > 0);
  gt_seedpairlist_at(&p_seedpair,seedpairlist,current - 1);
  gt_seedpairlist_at(&c_seedpair,seedpairlist,current);
  if (p_seedpair.aseqnum < c_seedpair.aseqnum)
  {
    return -1;
  }
  if (p_seedpair.aseqnum > c_seedpair.aseqnum)
  {
    return 1;
  }
  if (p_seedpair.bseqnum < c_seedpair.bseqnum)
  {
    return -1;
  }
  if (p_seedpair.bseqnum > c_seedpair.bseqnum)
  {
    return 1;
  }
  p_apos = GT_DIAGBANDSEED_GETPOS_A(&p_seedpair);
  c_apos = GT_DIAGBANDSEED_GETPOS_A(&c_seedpair);
  p_bpos = GT_DIAGBANDSEED_CONV_B(p_apos,GT_DIAGBANDSEED_GETPOS_B(&p_seedpair));
  c_bpos = GT_DIAGBANDSEED_CONV_B(c_apos,GT_DIAGBANDSEED_GETPOS_B(&c_seedpair));
  if (p_bpos < c_bpos)
  {
    return -1;
  }
  if (p_bpos > c_bpos)
  {
    return 1;
  }
  if (p_apos < c_apos)
  {
    return -1;
  }
  if (p_apos > c_apos)
  {
    return 1;
  }
  return 0;
}
#endif

static void gt_diagbandseed_seedpairlist_out(FILE *stream,
                                             const GtSeedpairlist *seedpairlist)
{
  GtUword spidx, mlistlen = gt_seedpairlist_length(seedpairlist);

  for (spidx = 0; spidx < mlistlen; spidx++)
  {
    gt_assert(spidx == 0 ||
              gt_diagbandseed_seeds_compare(seedpairlist,spidx) < 0);
    fprintf(stream, "# SeedPair ");
    gt_diagbandseed_show_seed(stream,seedpairlist,spidx);
    fprintf(stream, "\n");
  }
}

/* Fill a GtDiagbandseedSeedPair list of equal kmers from the iterators. */
static void gt_diagbandseed_merge(GtSeedpairlist *seedpairlist,
                                  GtUword *histogram,
                                  bool knowthesize,
                                  GtDiagbandseedKmerIterator *aiter,
                                  GtDiagbandseedKmerIterator *biter,
                                  GtUword maxfreq,
                                  GtUword maxgram,
                                  const GtRange *seedpairdistance,
                                  bool selfcomp)
{
  const GtKmerPosList *alist, *blist;
  const bool count_cartesian = (histogram != NULL && !selfcomp) ? true : false;

  gt_assert(aiter != NULL && biter != NULL &&
            ((histogram == NULL && seedpairlist != NULL) ||
            (histogram != NULL && seedpairlist == NULL)));
  alist = gt_diagbandseed_kmer_iter_next(aiter);
  blist = gt_diagbandseed_kmer_iter_next(biter);
  while (alist != NULL && blist != NULL)
  {
    const GtDiagbandseedKmerPos *asegment = alist->spaceGtDiagbandseedKmerPos,
                                *bsegment = blist->spaceGtDiagbandseedKmerPos;
    GtUword alen = alist->nextfree,
            blen = blist->nextfree;
    if (asegment->code < bsegment->code)
    {
      alist = gt_diagbandseed_kmer_iter_next(aiter);
    } else
    {
      if (asegment->code > bsegment->code)
      {
        blist = gt_diagbandseed_kmer_iter_next(biter);
      } else
      {
        GtUword frequency = GT_MAX(alen, blen);

        if (frequency <= maxfreq)
        {
          /* add all equal k-mers */
          frequency = GT_MIN(maxgram, frequency);
          gt_assert(frequency > 0);
          if (count_cartesian)
          {
            histogram[frequency - 1] += alen * blen;
          } else
          {
            const GtDiagbandseedKmerPos *aptr, *bptr;
            for (aptr = asegment; aptr < asegment + alen; aptr++)
            {
              for (bptr = bsegment; bptr < bsegment + blen; bptr++)
              {
                if (!selfcomp || aptr->seqnum < bptr->seqnum ||
                    (aptr->seqnum == bptr->seqnum &&
                     aptr->endpos + seedpairdistance->start <= bptr->endpos &&
                     aptr->endpos + seedpairdistance->end >= bptr->endpos))
                {
                  /* no duplicates from the same dataset */
                  if (histogram == NULL)
                  {
                    /* save SeedPair in seedpairlist */
                    gt_seedpairlist_add(seedpairlist,
                                        knowthesize,
                                        aptr->seqnum,
                                        bptr->seqnum,
                                        bptr->endpos,
                                        aptr->endpos);
                  } else
                  {
                    /* count seed frequency in histogram */
                    histogram[frequency - 1]++;
                  }
                }
              }
            }
          }
        } /* else: ignore all equal elements */
        alist = gt_diagbandseed_kmer_iter_next(aiter);
        blist = gt_diagbandseed_kmer_iter_next(biter);
      }
    }
  }
}

/* Verify seeds in the original sequences */
static int gt_diagbandseed_verify(const GtSeedpairlist *seedpairlist,
                                  const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  unsigned int seedlength,
                                  bool reverse,
                                  bool verbose,
                                  FILE *stream,
                                  GtError *err) {
  const GtUword mlistlen = gt_seedpairlist_length(seedpairlist);
  GtTimer *timer = NULL;
  GtUword idx;
  char *buf1 = gt_malloc(3 * (seedlength + 1) * sizeof *buf1);
  char *buf2 = buf1 + 1 + seedlength;
  char *buf3 = buf2 + 1 + seedlength;
  buf1[seedlength] = buf2[seedlength] = buf3[seedlength] = '\0';

  if (verbose) {
    fprintf(stream, "# start verifying seeds ...\n");
    timer = gt_timer_new();
    gt_timer_start(timer);
  }

  gt_assert(aencseq != NULL && bencseq != NULL);
  for (idx = 0; idx < mlistlen; idx++)
  {
    GtDiagbandseedSeedPair seedpair;
    GtDiagbandseedPosition abs_apos, abs_bpos;

    gt_seedpairlist_at(&seedpair,seedpairlist,idx);
    /* extract decoded k-mers at seed positions */
    abs_apos = GT_DIAGBANDSEED_GETPOS_A(&seedpair) +
               gt_encseq_seqstartpos(aencseq, seedpair.aseqnum);
    gt_encseq_extract_decoded(aencseq, buf1, abs_apos + 1 - seedlength,
                              abs_apos);

    if (!reverse) {
      GtDiagbandseedPosition apos = GT_DIAGBANDSEED_GETPOS_A(&seedpair),
                             bpos = GT_DIAGBANDSEED_GETPOS_B(&seedpair);

      abs_bpos = GT_DIAGBANDSEED_CONV_B(apos,bpos) +
                 gt_encseq_seqstartpos(bencseq,seedpair.bseqnum);
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos + 1 - seedlength,
                                abs_bpos);
      if (strcmp(buf1, buf2) != 0) {
        gt_error_set(err, "Wrong SeedPair (%" PRIu32
                          ",%" PRIu32 ",%" PRIu32 ",%" GT_WUS "): %s != %s\n",
                     seedpair.aseqnum, seedpair.bseqnum,
                     apos,
                     GT_DIAGBANDSEED_CONV_B(apos,bpos),
                     buf1, buf2);
        gt_free(buf1);
        if (verbose)
        {
          gt_timer_delete(timer);
        }
        return -1;
      }
    } else {
      /* get reverse k-mer */
      char *bufptr;
      abs_bpos = gt_encseq_seqstartpos(bencseq, seedpair.bseqnum) +
                 gt_encseq_seqlength(bencseq, seedpair.bseqnum)
                 - GT_DIAGBANDSEED_GETPOS_B(&seedpair) - 1;
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos,
                                abs_bpos + seedlength - 1);

      for (bufptr = buf3; bufptr < buf3 + seedlength; bufptr++) {
        gt_complement(bufptr, buf2[seedlength + buf3 - bufptr - 1], NULL);
      }
      if (strcmp(buf1, buf3) != 0) {
        gt_error_set(err, "Wrong SeedPair (" "%"PRIu32
                          ",%"PRIu32",%"PRIu32",%"PRIu32"): %s != %s\n",
                     seedpair.aseqnum, seedpair.bseqnum,
                     GT_DIAGBANDSEED_GETPOS_A(&seedpair),
                     GT_DIAGBANDSEED_GETPOS_B(&seedpair), buf1, buf3);
        gt_free(buf1);
        if (verbose)
        {
          gt_timer_delete(timer);
        }
        return -1;
      }
    }
  }
  if (verbose) {
    fprintf(stream, "# ...successfully verified all seeds ");
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }
  gt_free(buf1);
  return 0;
}

/* Return estimated length of mlist, and maxfreq w.r.t. the given memlimit */
static int gt_diagbandseed_get_mlistlen_maxfreq(GtUword *mlistlen,
                                            GtUword *maxfreq,
                                            GtDiagbandseedKmerIterator *aiter,
                                            GtDiagbandseedKmerIterator *biter,
                                            GtUword memlimit,
                                            size_t sizeofunit,
                                            const GtRange *seedpairdistance,
                                            GtUword len_used,
                                            bool selfcomp,
                                            bool alist_blist_id,
                                            bool verbose,
                                            FILE *stream,
                                            GtError *err)
{
  const GtUword maxgram = GT_MIN(*maxfreq, 8190) + 1; /* Cap on k-mer count */
  GtUword *histogram = NULL;
  GtTimer *timer = NULL;
  int had_err = 0;

  gt_assert(memlimit < GT_UWORD_MAX);
  if (verbose) {
    fprintf(stream, "# start calculating k-mer frequency histogram...\n");
    timer = gt_timer_new();
    gt_timer_start(timer);
  }

  /* build histogram; histogram[maxgram] := estimation for mlistlen */
  histogram = gt_calloc(maxgram + 1, sizeof *histogram);
  gt_diagbandseed_merge(NULL, /* mlist not needed: just count */
                        histogram,
                        false,
                        aiter,
                        biter,
                        *maxfreq,
                        maxgram,
                        seedpairdistance,
                        selfcomp);
  *maxfreq = gt_diagbandseed_processhistogram(histogram,
                                              *maxfreq,
                                              maxgram,
                                              memlimit,
                                              len_used *
                                                sizeof (GtDiagbandseedKmerPos),
                                              alist_blist_id,
                                              sizeofunit);
  *mlistlen = histogram[maxgram];
  gt_free(histogram);

  if (verbose) {
    gt_timer_show_formatted(timer,"# ... finished histogram "
                            GT_DIAGBANDSEED_FMT,stream);
    gt_timer_delete(timer);
  }

  /* check maxfreq value */
  if (*maxfreq == 0 || (*maxfreq == 1 && alist_blist_id)) {
    gt_error_set(err,
                 "option -memlimit too strict: need at least " GT_WU "MB",
                 (*mlistlen >> 20) + 1);
    *mlistlen = 0;
    had_err = -1;
  } else if (verbose) {
    if (*maxfreq == GT_UWORD_MAX) {
      fprintf(stream, "# disable k-mer maximum frequency, ");
    } else {
      fprintf(stream, "# set k-mer maximum frequency to " GT_WU ", ", *maxfreq);
    }
    fprintf(stream, "expect " GT_WU " seeds.\n", *mlistlen);
  } else if (*maxfreq <= 5) {
    gt_warning("only k-mers occurring <= " GT_WU " times will be considered, "
               "due to small memlimit.", *maxfreq);
  }

  return had_err;
}

/* Return a sorted list of SeedPairs from given Kmer-Iterators.
 * Parameter known_size > 0 can be given to allocate the memory beforehand.
 * The caller is responsible for freeing the result. */
static void gt_diagbandseed_get_seedpairs(GtSeedpairlist *seedpairlist,
                                          GtDiagbandseedKmerIterator *aiter,
                                          GtDiagbandseedKmerIterator *biter,
                                          GtUword maxfreq,
                                          GtUword known_size,
                                          const GtRange *seedpairdistance,
                                          bool selfcomp,
                                          bool debug_seedpair,
                                          bool verbose,
                                          FILE *stream)
{
  GtTimer *timer = NULL;
  GtUword mlistlen;

  if (verbose) {
    timer = gt_timer_new();
    if (known_size > 0) {
      fprintf(stream, "# start collecting " GT_WU " seeds in %.0f MB ...\n",
                      known_size,
                      GT_MEGABYTES(known_size *
                                   gt_seedpairlist_sizeofunit(seedpairlist)));
    } else {
      fprintf(stream, "# start collecting seeds ...\n");
    }
    gt_timer_start(timer);
  }

  /* allocate mlist space according to seed count */
  gt_seedpairlist_init(seedpairlist,known_size);

  /* create mlist */
  (void) gt_diagbandseed_merge(seedpairlist,
                        NULL, /* histogram not needed: save seeds */
                        known_size > 0 ? true : false,
                        aiter,
                        biter,
                        maxfreq,
                        GT_UWORD_MAX, /* maxgram */
                        seedpairdistance,
                        selfcomp);
  mlistlen = gt_seedpairlist_length(seedpairlist);
  if (verbose) {
    fprintf(stream, "# ... collected " GT_WU " seeds ", mlistlen);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
  }

  /* sort mlist */
  if (mlistlen > 0) {
    if (verbose) {
      gt_timer_start(timer);
    }
    gt_diagbandseed_seedpairlist_sort(seedpairlist);
    if (verbose) {
      fprintf(stream, "# ... sorted " GT_WU " seeds ", mlistlen);
      gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    }
    if (debug_seedpair) {
      gt_diagbandseed_seedpairlist_out(stream,seedpairlist);
    }
  }
  if (timer != NULL)
  {
    gt_timer_delete(timer);
  }
}

typedef struct
{
  GtUword totalseeds,
          totalseqpairs,
          extended_seeds,
          selected_seeds,
          countmatches,
          failedmatches,
          seqpairs_with_minsegment,
          maxmatchespersegment,
          total_extension_time_usec,
          total_process_seeds_usec,
          filteredbydiagonalscore;
  bool withtiming;
  GtBittab *used_a_sequences,
           *used_b_sequences;
  GtUword current_gfa2_edge_num;
} GtDiagbandseedState;

static GtDiagbandseedState *gt_diagbandseed_dbs_state_new(bool verbose,
                                                        GtUword a_num_sequences,
                                                        GtUword b_num_sequences)
{
  GtDiagbandseedState *dbs_state = gt_calloc(1,sizeof *dbs_state);
  dbs_state->withtiming = verbose;
  if (a_num_sequences > 0)
  {
    dbs_state->used_a_sequences = gt_bittab_new(a_num_sequences);
  } else
  {
    dbs_state->used_a_sequences = NULL;
  }
  if (b_num_sequences > 0)
  {
    gt_assert(a_num_sequences > 0);
    dbs_state->used_b_sequences = gt_bittab_new(b_num_sequences);
  } else
  {
    dbs_state->used_b_sequences = NULL;
  }
  return dbs_state;
}

static void gt_diagbandseed_dbs_state_delete(GtDiagbandseedState *dbs_state)
{
  if (dbs_state != NULL)
  {
    gt_bittab_delete(dbs_state->used_a_sequences);
    gt_bittab_delete(dbs_state->used_b_sequences);
    gt_free(dbs_state);
  }
}

static void gt_diagbandseed_dbs_state_update(GtDiagbandseedState *dbs_state,
                                          GtUword mlistlen,
                                          GtUword numseqpairs)
{
  dbs_state->totalseeds += mlistlen;
  dbs_state->totalseqpairs += numseqpairs;
}

static void gt_diagbandseed_dbs_state_out(const GtDiagbandseedState *dbs_state,
                                       bool maxmat_show)
{
  printf("# total number of seeds: " GT_WU,dbs_state->totalseeds);
  if (!maxmat_show)
  {
    printf("; " GT_WU ", i.e. %.2f%% of all seeds were filtered by "
                   "diagonal score"
                   "; " GT_WU ", i.e. %.2f%% of all seeds were selected"
                   "; " GT_WU ", i.e. %.2f%% of all seeds were extended\n",
            dbs_state->filteredbydiagonalscore,
            100.0 * (double) dbs_state->filteredbydiagonalscore/
                             dbs_state->totalseeds,
            dbs_state->selected_seeds,
            100.0 * (double) dbs_state->selected_seeds/
                             dbs_state->totalseeds,
            dbs_state->extended_seeds,
            100.0 * (double) dbs_state->extended_seeds/
                             dbs_state->totalseeds);
    printf( "# number of unsuccessful extensions: " GT_WU "\n",
            dbs_state->failedmatches);
    printf( "# sequence pairs with selected seeds: " GT_WU
                    " (%.2f%% of all " GT_WU ")\n",
            dbs_state->seqpairs_with_minsegment,
            100.0 * (double) dbs_state->seqpairs_with_minsegment/
                             dbs_state->totalseqpairs,
            dbs_state->totalseqpairs);
  } else
  {
    printf("\n# maximum number of MEMs per segment: " GT_WU "\n",
            dbs_state->maxmatchespersegment);
  }
  printf("# number of matches output: " GT_WU " (%.4f per seed)\n",
          dbs_state->countmatches,
          (double) dbs_state->countmatches/
                   dbs_state->totalseeds);
}

typedef bool (*GtExtendRelativeCoordsFunc)(void *,
                                           const GtSeqorEncseq *,
                                           GtUword,
                                           GtUword,
                                           const GtSeqorEncseq *,
                                           bool,
                                           GtUword,
                                           GtUword,
                                           GtUword,
                                           GtReadmode);

static bool gt_diagbandseed_has_overlap_with_previous_match(
     const GtArrayGtDiagbandseedRectangle *previous_extensions,
     GtUword previous_match_a_start,
     GtUword previous_match_a_end,
     GtUword previous_match_b_start,
     GtUword previous_match_b_end,
     GtUword apos,
     GtUword bpos,
     GtUword matchlength,
     bool debug)
{
  GtDiagbandseedRectangle maxmatch;

  if (debug)
  {
    printf("# overlap of " GT_WU " " GT_WU " " GT_WU " " GT_WU " "
                           GT_WU " " GT_WU " " GT_WU " " GT_WU "?\n",
           apos + 1 - matchlength,apos,bpos + 1 - matchlength,bpos,
           previous_match_a_start,
           previous_match_a_end,
           previous_match_b_start,
           previous_match_b_end);
  }

  maxmatch.a_start = apos + 1 - matchlength;
  maxmatch.a_end = apos;
  maxmatch.b_start = bpos + 1 - matchlength;
  maxmatch.b_end = bpos;

  if (debug)
  {
    gt_rectangle_store_show(previous_extensions);
  }
  if (gt_rectangle_overlap(previous_extensions,&maxmatch))
  {
    if (debug)
    {
      printf("# overlap with previous_ext (both dim) => return true\n");
    }
    return true;
  } else
  {
    if (debug)
    {
      printf("# not overlap with previous_ext (both dim) => return false\n");
    }
    return false;
  }
}

typedef struct
{
  bool b_differs_from_a, a_haswildcards, b_haswildcards;
  const GtUchar *characters;
  GtUchar wildcardshow;
  GtSeqorEncseq aseqorencseq, bseqorencseq;
  GtUchar *a_byte_sequence, *b_byte_sequence;
  GtUword previous_aseqnum,
          a_first_seqnum,
          a_first_seqstartpos,
          b_first_seqnum,
          b_first_seqstartpos;
  const GtEncseq *a_encseq_for_seq_desc, *b_encseq_for_seq_desc;
} GtDiagbandSeedPlainSequence;

static void gt_diagbandseed_plainsequence_init(GtDiagbandSeedPlainSequence *ps,
                                          bool s_desc_display,
                                          bool q_desc_display,
                                          const GtEncseq *aencseq,
                                          const GtSequencePartsInfo *aseqranges,
                                          GtUword aidx,
                                          bool with_a_bytestring,
                                          const GtEncseq *bencseq,
                                          const GtSequencePartsInfo *bseqranges,
                                          GtUword bidx,
                                          bool with_b_bytestring)
{
  ps->previous_aseqnum = GT_UWORD_MAX;
  if (s_desc_display && aencseq != NULL &&
      gt_encseq_has_description_support(aencseq))
  {
    ps->a_encseq_for_seq_desc = aencseq;
  } else
  {
    ps->a_encseq_for_seq_desc = NULL;
  }
  if (q_desc_display && bencseq != NULL &&
      gt_encseq_has_description_support(bencseq))
  {
    ps->b_encseq_for_seq_desc = bencseq;
  } else
  {
    ps->b_encseq_for_seq_desc = NULL;
  }
  if (with_a_bytestring)
  {
    ps->a_first_seqnum = gt_sequence_parts_info_start_get(aseqranges,aidx);
    ps->a_first_seqstartpos
      = gt_sequence_parts_info_seqstartpos(aseqranges,ps->a_first_seqnum),
    ps->a_byte_sequence = gt_sequence_parts_info_seq_extract(aencseq,aseqranges,
                                                             aidx);
    ps->a_haswildcards = (gt_encseq_wildcards(aencseq) > 0) ? true : false;
  } else
  {
    ps->a_byte_sequence = NULL;
    GT_SEQORENCSEQ_INIT_ENCSEQ(&ps->aseqorencseq,aencseq);
    ps->a_haswildcards = true;
  }
  if (with_b_bytestring)
  {
    ps->b_first_seqnum = gt_sequence_parts_info_start_get(bseqranges,bidx);
    ps->b_first_seqstartpos
      = gt_sequence_parts_info_seqstartpos(bseqranges,ps->b_first_seqnum);
    if (with_a_bytestring && aencseq == bencseq && aidx == bidx)
    {
      ps->b_differs_from_a = false;
      ps->b_byte_sequence = ps->a_byte_sequence;
      ps->b_haswildcards = ps->a_haswildcards;
    } else
    {
      ps->b_differs_from_a = true;
      ps->b_byte_sequence = gt_sequence_parts_info_seq_extract(bencseq,
                                                               bseqranges,bidx);
      ps->b_haswildcards = (gt_encseq_wildcards(bencseq) > 0) ? true : false;
    }
  } else
  {
    ps->b_byte_sequence = NULL;
    GT_SEQORENCSEQ_INIT_ENCSEQ(&ps->bseqorencseq,bencseq);
    ps->b_haswildcards = true;
  }
  if (with_a_bytestring || with_b_bytestring)
  {
    ps->characters = gt_encseq_alphabetcharacters(aencseq);
    ps->wildcardshow = gt_alphabet_wildcard_show(gt_encseq_alphabet(aencseq));
  }
}

static void gt_diagbandseed_set_sequence(GtSeqorEncseq *seqorencseq,
                                         const GtSequencePartsInfo *seqranges,
                                         const GtUchar *bytesequence,
                                         GtUword first_seqstartpos,
                                         GtUword seqnum,
                                         const GtUchar *characters,
                                         GtUchar wildcardshow,
                                         bool haswildcards,
                                         const GtEncseq *encseq_for_seq_desc)
{
  const GtUword
    seqstartpos = gt_sequence_parts_info_seqstartpos(seqranges,seqnum),
    seqendpos = gt_sequence_parts_info_seqendpos(seqranges,seqnum),
    b_off = seqstartpos - first_seqstartpos;
  const char *seqdescptr;
  if (encseq_for_seq_desc != NULL)
  {
    GtUword desclen;
    seqdescptr = gt_encseq_description(encseq_for_seq_desc,
                                       &desclen,
                                       seqnum);
  } else
  {
    seqdescptr = "Unknown";
  }
  GT_SEQORENCSEQ_INIT_SEQ(seqorencseq,bytesequence + b_off,seqdescptr,
                          seqendpos - seqstartpos + 1,characters,
                          wildcardshow,
                          haswildcards);
}

static void gt_diagbandseed_plainsequence_next_segment(
                         GtDiagbandSeedPlainSequence *ps,
                         const GtSequencePartsInfo *aseqranges,
                         GtUword currsegm_aseqnum,
                         const GtSequencePartsInfo *bseqranges,
                         GtUword currsegm_bseqnum)
{
  if (ps->previous_aseqnum == GT_UWORD_MAX ||
      ps->previous_aseqnum < currsegm_aseqnum)
  {
    if (ps->a_byte_sequence != NULL)
    {
      gt_diagbandseed_set_sequence(&ps->aseqorencseq,
                                   aseqranges,
                                   ps->a_byte_sequence,
                                   ps->a_first_seqstartpos,
                                   currsegm_aseqnum,
                                   ps->characters,
                                   ps->wildcardshow,
                                   ps->a_haswildcards,
                                   ps->a_encseq_for_seq_desc);
      ps->previous_aseqnum = currsegm_aseqnum;
    } else
    {
      GtUword seqstartpos = gt_sequence_parts_info_seqstartpos(aseqranges,
                                                             currsegm_aseqnum),
              seqendpos = gt_sequence_parts_info_seqendpos(aseqranges,
                                                           currsegm_aseqnum);

      GT_SEQORENCSEQ_ADD_SEQ_COORDS(&ps->aseqorencseq,
                                    seqstartpos,
                                    seqendpos - seqstartpos + 1);
    }
  } else
  {
    gt_assert(ps->previous_aseqnum == currsegm_aseqnum);
  }
  if (ps->b_byte_sequence != NULL)
  {
    gt_diagbandseed_set_sequence(&ps->bseqorencseq,
                                 bseqranges,
                                 ps->b_byte_sequence,
                                 ps->b_first_seqstartpos,
                                 currsegm_bseqnum,
                                 ps->characters,
                                 ps->wildcardshow,
                                 ps->b_haswildcards,
                                 ps->b_encseq_for_seq_desc);
  } else
  {
    GtUword seqstartpos = gt_sequence_parts_info_seqstartpos(bseqranges,
                                                             currsegm_bseqnum),
            seqendpos = gt_sequence_parts_info_seqendpos(bseqranges,
                                                         currsegm_bseqnum);

    GT_SEQORENCSEQ_ADD_SEQ_COORDS(&ps->bseqorencseq,
                                  seqstartpos,
                                  seqendpos - seqstartpos + 1);
  }
}

static void gt_diagbandseed_plainsequence_delete(
                         GtDiagbandSeedPlainSequence *ps)
{
  if (ps->b_byte_sequence != NULL && ps->b_differs_from_a)
  {
    gt_free(ps->b_byte_sequence);
  }
  if (ps->a_byte_sequence != NULL)
  {
    gt_free(ps->a_byte_sequence);
  }
}

typedef struct
{
  GtUword userdefinedleastlength,
          errorpercentage,
          use_apos,
          mincoverage;
  bool only_selected_seqpairs;
  double evalue_threshold;
  GtDiagbandSeedPlainSequence plainsequence_info;
  GtReadmode query_readmode;
  bool same_encseq, debug;
  GtProcessinfo_and_querymatchspaceptr info_querymatch;
  GtExtendRelativeCoordsFunc extend_relative_coords_function;
  GtSegmentRejectFunc segment_reject_func;
  GtSegmentRejectInfo *segment_reject_info;
  const GtKarlinAltschulStat *karlin_altschul_stat;
  const GtSeedExtendDisplayFlag *out_display_flag;
  bool benchmark;
  GtAniAccumulate *ani_accumulate;
  GtDiagbandseedState *dbs_state;
} GtDiagbandseedExtendSegmentInfo;

static int gt_diagbandseed_possibly_extend(const GtArrayGtDiagbandseedRectangle
                                             *previous_extensions,
                                           bool haspreviousmatch,
                                           GtUword use_apos,
                                           GtUword aseqnum,
                                           GtUword apos,
                                           GtUword bseqnum,
                                           GtUword bpos,
                                           GtUword matchlength,
                                           GtDiagbandseedExtendSegmentInfo *esi)
{
  int ret = 0;

  if (esi->debug)
  {
    printf("# %s haspreviousmatch=%s\n",__func__,
           haspreviousmatch ? "true" : "false");
  }
  if (!haspreviousmatch ||
      (use_apos == 0 &&
       esi->info_querymatch.previous_match_b_end < bpos) || /* no overlap */
      (use_apos > 0 && !gt_diagbandseed_has_overlap_with_previous_match(
                             previous_extensions,
                             esi->info_querymatch.previous_match_a_start,
                             esi->info_querymatch.previous_match_a_end,
                             esi->info_querymatch.previous_match_b_start,
                             esi->info_querymatch.previous_match_b_end,
                             apos,
                             bpos,
                             matchlength,
                             esi->debug)))
  {
    bool success;

    /* relative seed start position in A and B */
    const GtUword bstart = bpos + 1 - matchlength;
    const GtUword astart = apos + 1 - matchlength;
#ifndef _WIN32
    struct timeval tvalBefore = {0};

    if (esi->dbs_state != NULL &&
        esi->dbs_state->withtiming)
    {
      gettimeofday (&tvalBefore, NULL);
    }
#endif
    ret = 1; /* perform extension */
    /* the following function called is either
         gt_greedy_extend_seed_relative or
         gt_xdrop_extend_seed_relative */
    success = esi->extend_relative_coords_function(&esi->info_querymatch,
                                                   &esi->plainsequence_info.
                                                     aseqorencseq,
                                                   aseqnum,
                                                   astart,
                                                   &esi->plainsequence_info.
                                                     bseqorencseq,
                                                   esi->same_encseq,
                                                   bseqnum,
                                                   bstart,
                                                   matchlength,
                                                   esi->query_readmode);
#ifndef _WIN32
    if (esi->dbs_state != NULL &&
        esi->dbs_state->withtiming)
    {
      struct timeval tvalAfter;
      gettimeofday (&tvalAfter, NULL);
      esi->dbs_state->total_extension_time_usec
        += (tvalAfter.tv_sec - tvalBefore.tv_sec) * 1000000L
            + tvalAfter.tv_usec - tvalBefore.tv_usec;
    }
#endif
    if (success)
    {
      double evalue, bit_score;

      if (esi->ani_accumulate != NULL)
      {
        const GtUword query_seqlen
          = esi->plainsequence_info.bseqorencseq.seqlength,
        aligned_len
          = esi->info_querymatch.previous_match_a_end -
            esi->info_querymatch.previous_match_a_start + 1 +
            esi->info_querymatch.previous_match_b_end -
            esi->info_querymatch.previous_match_b_start + 1;
        if (gt_querymatch_check_final_generic(
                               &evalue,
                               &bit_score,
                               esi->karlin_altschul_stat,
                               query_seqlen,
                               aligned_len,
                               esi->info_querymatch.previous_match_distance,
                               esi->info_querymatch.previous_match_mismatches,
                               esi->userdefinedleastlength,
                               esi->errorpercentage,
                               esi->evalue_threshold,
                               stdout))
        {
          esi->ani_accumulate->sum_of_aligned_len += aligned_len;
          esi->ani_accumulate->sum_of_distance
            += esi->info_querymatch.previous_match_distance;
          ret = 3;
        } else
        {
          ret = 2; /* found match, which does not satisfy length or similarity
                      constraints */
        }
      } else
      {
        const GtQuerymatch *querymatch
          = esi->info_querymatch.querymatchspaceptr;

        /* show extension results */
        gt_assert(querymatch != NULL);
        if (gt_querymatch_check_final(&evalue,
                                      &bit_score,
                                      esi->karlin_altschul_stat,
                                      querymatch,
                                      esi->userdefinedleastlength,
                                      esi->errorpercentage,
                                      esi->evalue_threshold))
        {
          if (!esi->benchmark) {
            if (gt_querymatch_gfa2_display(esi->out_display_flag))
            {
              gt_assert(esi->dbs_state != NULL);
              gt_querymatch_gfa2_edge(querymatch,
                                      esi->dbs_state->
                                           current_gfa2_edge_num);
              esi->dbs_state->current_gfa2_edge_num++;
              gt_bittab_set_bit(esi->dbs_state->used_a_sequences,
                                aseqnum);
              gt_bittab_set_bit(
                 esi->dbs_state->used_b_sequences != NULL
                   ? esi->dbs_state->used_b_sequences
                   : esi->dbs_state->used_a_sequences,
                   bseqnum);
            }
            gt_querymatch_prettyprint(evalue,bit_score,esi->out_display_flag,
                                      querymatch);
          }
          ret = 3; /* output match */
        } else
        {
          if (!esi->benchmark) {
            gt_querymatch_show_failed_seed(esi->out_display_flag,querymatch);
          }
          ret = 2; /* found match, which does not satisfy length or similarity
                      constraints */
        }
      }
    }
    /* else reference and query are the same sequence and overlap so that
       no extension was performed */
  }
  return ret;
}

void gt_diagbandseed_printchainelem(FILE *outfp,
                                    GtUword aseqnum,
                                    GtUword astart,
                                    GtUword bseqnum,
                                    GtUword bstart,
                                    GtUword matchlength)
{
  fprintf(outfp,GT_WU " " GT_WU " " GT_WU " " GT_WU " " GT_WU " " GT_WU "\n",
          matchlength,aseqnum,astart,matchlength,bseqnum,bstart);
}

typedef struct
{
  GtDiagbandseedSeqnum aseqnum,
                       bseqnum;
} GtDiagbandseedSequencePair;

static void gt_diagbandseed_chain_out(void *data,
                                      const GtChain2Dimmatchtable *matchtable,
                                      const GtChain2Dim *chain)
{
  GtUword idx, chainlength;
  const GtDiagbandseedSequencePair *sequencepair
    = (const GtDiagbandseedSequencePair *) data;

  gt_assert(chain != NULL);
  chainlength = gt_chain_chainlength(chain);
  printf("# chain of length " GT_WU " with score " GT_WD "\n",
         chainlength,gt_chain_chainscore(chain));

  gt_assert(!gt_chain_storedinreverseorder(chain));
  for (idx = 0; idx < chainlength; idx++)
  {
    GtChain2Dimmatchvalues value;

    gt_chain_extractchainelem(&value, matchtable, chain, idx);
    gt_assert(value.startpos[0] <= value.endpos[0] &&
              value.weight >= 0 &&
              value.endpos[0] - value.startpos[0] + 1
                        == (GtUword) value.weight &&
              value.startpos[1] <= value.endpos[1] &&
              value.endpos[1] - value.startpos[1] + 1
                        == (GtUword) value.weight);

    gt_diagbandseed_printchainelem(stdout,
                                   sequencepair->aseqnum,
                                   value.startpos[0],
                                   sequencepair->bseqnum,
                                   value.startpos[1],
                                   (GtUword) value.weight);
  }
}

static int gt_diagbandseed_bstart_ldesc_compare_mems(const void *vl,
                                                     const void *vr)
{
  GtDiagbandseedMaximalmatch *l = (GtDiagbandseedMaximalmatch *) vl,
                             *r = (GtDiagbandseedMaximalmatch *) vr;
  GtDiagbandseedPosition lapos = l->apos + 1 - l->len,
                         lbpos = l->bpos + 1 - l->len,
                         rapos = r->apos + 1 - r->len,
                         rbpos = r->bpos + 1 - r->len;
  if (lbpos < rbpos)
  {
    return -1;
  }
  if (lbpos > rbpos)
  {
    return 1;
  }
  if (l->len < r->len)
  {
    return 1;
  }
  if (l->len > r->len)
  {
    return -1;
  }
  if (lapos < rapos)
  {
    return -1;
  }
  if (lapos > rapos)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

typedef struct
{
  GtDiagbandseedMaximalmatch *spaceGtDiagbandseedMaximalmatch;
  GtUword allocatedGtDiagbandseedMaximalmatch,
          nextfreeGtDiagbandseedMaximalmatch;
} GtArrayGtDiagbandseedMaximalmatch;

static bool gt_diagbandseed_process_mem(
                                    bool forward,
                                    GtUword aseqnum,
                                    GtUword bseqnum,
                                    GtArrayGtDiagbandseedMaximalmatch *memstore,
                                    GtUword amaxlen,
                                    const GtSeedpairPositions *previous,
                                    GtUword previous_matchlength,
                                    GtUword userdefinedleastlength,
                                    FILE *fpout)
{
  if (memstore == NULL)
  {
    if (previous_matchlength >= userdefinedleastlength)
    {
      fprintf(fpout,"%8" GT_WUS "  %8" GT_WUS "  %8" GT_WUS "  %c  %8" GT_WUS
                    "  %8" GT_WUS "\n",
                  previous_matchlength,
                  aseqnum,
                  previous->apos + 2 - previous_matchlength,
                  forward ? 'F' : 'P',
                  bseqnum,
                  GT_DIAGBANDSEED_DIAGONAL2BPOS(amaxlen,
                                                previous->apos,
                                                previous->bpos)
                                                + 2 - previous_matchlength);
      return true;
    }
  } else
  {
    GtDiagbandseedMaximalmatch *memstore_ptr;

    GT_GETNEXTFREEINARRAY(memstore_ptr,memstore,
                          GtDiagbandseedMaximalmatch,
                          memstore->allocatedGtDiagbandseedMaximalmatch
                          * 0.2 + 256);
    memstore_ptr->apos = previous->apos;
    memstore_ptr->bpos = GT_DIAGBANDSEED_DIAGONAL2BPOS(amaxlen,
                                                       previous->apos,
                                                       previous->bpos);
    memstore_ptr->len = previous_matchlength;
  }
  return false;
}

/* This will later be replaced by an encoding of the three
   values in a GtUword as to speed up sorting */

static void gt_diagbandseed_segment2maxmatches(
             bool forward,
             GtArrayGtDiagbandseedMaximalmatch *memstore,
             GtUword aseqnum,
             GtUword bseqnum,
             unsigned int seedlength,
             GtUword userdefinedleastlength,
             GtUword amaxlen,
             const GtSeedpairPositions *segment_positions,
             GtUword segment_length,
             GtDiagbandseedState *dbs_state,
             GtSegmentRejectFunc segment_reject_func,
             GtSegmentRejectInfo *segment_reject_info,
             const GtChain2Dimmode *chainmode,
             FILE *fpout)
{
  GtUword previous_matchlength = seedlength, localmatchcount = 0;
  const GtSeedpairPositions *current;
  GtSeedpairPositions previous;
  GtDiagbandseedSequencePair sequencepair = {aseqnum,bseqnum};
  const bool anchor_pairs = false;
  bool rejected = false;
#ifndef NDEBUG
  GtUword idx;
#endif

  gt_assert(segment_length > 0 && seedlength <= userdefinedleastlength);
#ifndef NDEBUG
  for (idx = 1; idx < segment_length; idx++)
  {
    gt_assert(segment_positions[idx-1].bpos < segment_positions[idx].bpos ||
              (segment_positions[idx-1].bpos == segment_positions[idx].bpos &&
               segment_positions[idx-1].apos <= segment_positions[idx].apos));
  }
#endif
  previous = segment_positions[0];
  if (memstore != NULL)
  {
    memstore->nextfreeGtDiagbandseedMaximalmatch = 0;
  }
  for (current = segment_positions + 1;
       current < segment_positions + segment_length && !rejected; current++)
  {
    if (previous.bpos == current->bpos && previous.apos + 1 == current->apos)
    {
      previous_matchlength++;
      previous.apos++;
    } else
    {
      if (previous.bpos == current->bpos &&
          previous.apos + seedlength - 1 >= current->apos)
      {
        /* This case can only happen if maxfreq exclude some intermediate
           matches */
        gt_assert(previous.apos <= current->apos);
        previous_matchlength += (current->apos - previous.apos + 1);
        previous.apos = current->apos;
      } else
      {
        /*
        printf("previous.diag = %u, current.diag=%u\n",previous.bpos,
                current->bpos);
        printf("previous.apos=%u,current.apos=%u,seedlength=%u\n",
                previous.apos,current->apos,seedlength);
        gt_assert (previous.bpos != current->bpos ||
                   previous.apos + seedlength - 1 < current->apos);*/
        if (anchor_pairs)
        {
          if (previous.bpos == current->bpos &&
              previous.apos + seedlength - 1 < current->apos)
          {
            printf("A %u\n",previous.apos+1);
          }
        } else
        {
          if (gt_diagbandseed_process_mem(forward,
                                          aseqnum,
                                          bseqnum,
                                          memstore,
                                          amaxlen,
                                          &previous,
                                          previous_matchlength,
                                          userdefinedleastlength,
                                          fpout))
          {
            localmatchcount++;
            if (segment_reject_func != NULL)
            {
              gt_segment_reject_register_match(segment_reject_info,bseqnum);
              rejected = true;
            }
          }
        }
        previous_matchlength = seedlength;
        previous = *current;
      }
    }
  }
  if (!anchor_pairs)
  {
    if (!rejected && gt_diagbandseed_process_mem(forward,
                                                 aseqnum,
                                                 bseqnum,
                                                 memstore,
                                                 amaxlen,
                                                 &previous,
                                                 previous_matchlength,
                                                 userdefinedleastlength,
                                                 fpout))
    {
      localmatchcount++;
      if (segment_reject_func != NULL)
      {
        gt_segment_reject_register_match(segment_reject_info,bseqnum);
      }
    }
  }
  if (dbs_state != NULL)
  {
    dbs_state->countmatches += localmatchcount;
    if (dbs_state->maxmatchespersegment < localmatchcount)
    {
      dbs_state->maxmatchespersegment = localmatchcount;
    }
  }
  if (memstore != NULL)
  {
    if (memstore->nextfreeGtDiagbandseedMaximalmatch >= 2)
    {
      qsort(memstore->spaceGtDiagbandseedMaximalmatch,
            memstore->nextfreeGtDiagbandseedMaximalmatch,
            sizeof *memstore->spaceGtDiagbandseedMaximalmatch,
            gt_diagbandseed_bstart_ldesc_compare_mems);
    }
    if (chainmode != NULL)
    {
      const GtUword presortdim = 1;
#ifndef NDEBUG
      GtUword previous_start_b = GT_UWORD_MAX;
#endif
      GtChain2Dimmatchvalues inmatch;
      GtChain2Dim *localchain = gt_chain_chain_new();
      GtDiagbandseedMaximalmatch *mmptr;
      GtChain2Dimmatchtable *chainmatchtable
        = gt_chain_matchtable_new(memstore->nextfreeGtDiagbandseedMaximalmatch);

      for (mmptr = memstore->spaceGtDiagbandseedMaximalmatch;
           mmptr < memstore->spaceGtDiagbandseedMaximalmatch +
                   memstore->nextfreeGtDiagbandseedMaximalmatch;
           mmptr++)
      {
        inmatch.startpos[0] = mmptr->apos + 1 - mmptr->len;
        inmatch.endpos[0] = mmptr->apos;
        inmatch.startpos[1] = mmptr->bpos + 1 - mmptr->len;
        inmatch.endpos[1] = mmptr->bpos;
        gt_assert(previous_start_b == GT_UWORD_MAX ||
                  previous_start_b <= inmatch.startpos[1]);
        inmatch.weight = mmptr->len;
        gt_chain_matchtable_add(chainmatchtable,&inmatch);
#ifndef NDEBUG
        previous_start_b = inmatch.startpos[1];
#endif
      }
      gt_chain_fastchaining(chainmode,
                            localchain,
                            chainmatchtable,
                            false,
                            presortdim,
                            true,
                            gt_diagbandseed_chain_out,
                            &sequencepair,
                            NULL);
      gt_chain_matchtable_delete(chainmatchtable);
      gt_chain_chain_delete(localchain);
    }
  }
}

/* * * * * SEED EXTENSION * * * * */

static void gt_diagbandseed_segment2matches(
             void *v_segment2matches_info,
             GT_UNUSED const GtEncseq *aencseq,
             GT_UNUSED const GtEncseq *bencseq,
             GtUword aseqnum,
             GtUword bseqnum,
             const GtDiagbandStruct *diagband_struct,
             const GtDiagbandseedMaximalmatch *memstore,
             unsigned int seedlength,
             const GtSeedpairPositions *segment_positions,
             GtUword segment_length)
{
  GtDiagbandseedExtendSegmentInfo *esi
    = (GtDiagbandseedExtendSegmentInfo *) v_segment2matches_info;
  bool found_selected = false, haspreviousmatch = false;
  GtUword idx, matchlength = seedlength;
  GtArrayGtDiagbandseedRectangle *previous_extensions = NULL;

  if (esi->use_apos > 0)
  {
    previous_extensions = gt_rectangle_store_new();
  }
  gt_assert(memstore != NULL || segment_positions != NULL);
  for (idx = 0; idx < segment_length; idx++)
  {
    GtDiagbandseedPosition apos, bpos;
    GtUword coverage;

    if (memstore != NULL)
    {
      apos = memstore[idx].apos;
      bpos = memstore[idx].bpos;
      matchlength = memstore[idx].len;
    } else
    {
      apos = segment_positions[idx].apos;
      bpos = segment_positions[idx].bpos;
    }
    if (esi->debug)
    {
      printf("# apos=%u,bpos=%u,matchlength=" GT_WU "\n",
             apos,bpos,matchlength);
    }

    /* The filter sums the score of the current diagonalband
       as well as the maximum of the score of the previous and next
       diagonal band */
    coverage = gt_diagband_struct_coverage(diagband_struct,apos,bpos);
    if (coverage >= esi->mincoverage)
    {
      int ret;

      if (esi->dbs_state != NULL)
      {
        esi->dbs_state->selected_seeds++;
      }
      if (!found_selected)
      {
        if (esi->dbs_state != NULL)
        {
          esi->dbs_state->seqpairs_with_minsegment++;
        }
        found_selected = true;
      }
      if (esi->only_selected_seqpairs)
      {
        printf("# " GT_WU "%c" GT_WU "\n",aseqnum,
               esi->query_readmode == GT_READMODE_REVCOMPL ? '-' : '+',
               bseqnum);
        break;
      }
      ret = gt_diagbandseed_possibly_extend(
                     previous_extensions,
                     haspreviousmatch,
                     esi->use_apos,
                     aseqnum,
                     apos,
                     bseqnum,
                     bpos,
                     matchlength,
                     esi);
      if (ret >= 1 && esi->dbs_state != NULL)
      {
        esi->dbs_state->extended_seeds++;
      }
      if (ret >= 2)
      {
        haspreviousmatch = true;
        if (esi->use_apos == 2 || /* add all previous matches */
            (esi->use_apos == 1 /* only add successful match */
             && ret == 3))
        {
          GtDiagbandseedRectangle newrectangle;

          newrectangle.a_start
            = esi->info_querymatch.previous_match_a_start;
          newrectangle.a_end = esi->info_querymatch.previous_match_a_end;
          newrectangle.b_start
            = esi->info_querymatch.previous_match_b_start;
          newrectangle.b_end
            = esi->info_querymatch.previous_match_b_end;
          if (esi->debug)
          {
            printf("# add rectangle (%u,%u) (%u,%u) from seed (%u,%u,"
                   GT_WU ")\n",
                    newrectangle.a_start,
                    newrectangle.a_end,
                    newrectangle.b_start,
                    newrectangle.b_end,
                    apos,bpos,matchlength);
          }
          gt_rectangle_store_add(previous_extensions,&newrectangle);
        }
        if (ret == 2 && esi->dbs_state != NULL)
        {
          esi->dbs_state->failedmatches++;
        }
      }
      if (ret == 3)
      {
        if (esi->dbs_state != NULL)
        {
          esi->dbs_state->countmatches++;
        }
        if (esi->segment_reject_func != NULL)
        {
          gt_segment_reject_register_match(esi->segment_reject_info,bseqnum);
          break;
        }
      }
    } else
    {
      if (esi->dbs_state != NULL)
      {
        esi->dbs_state->filteredbydiagonalscore++;
      }
      if (esi->debug)
      {
        printf("# filtered as diagonal score " GT_WU " < " GT_WU "\n",
               coverage,
               esi->mincoverage);
      }
    }
  }
  if (esi->use_apos > 0)
  {
    gt_rectangle_store_delete(previous_extensions);
  }
}

#define GT_DIAGBANDSEED_PROCESS_SEGMENT\
        if (segment_reject_func == NULL ||\
            !segment_reject_func(segment_reject_info,currsegm_bseqnum))\
        {\
          if (seedpairlist->maxmat_compute)\
          {\
            gt_diagbandseed_segment2maxmatches(forward,\
                                               memstore,\
                                               currsegm_aseqnum,\
                                               currsegm_bseqnum,\
                                               seedlength,\
                                               extp->userdefinedleastlength,\
                                               seedpairlist->amaxlen,\
                                               segment_positions,\
                                               segment_length,\
                                               dbs_state,\
                                               segment_reject_func,\
                                               segment_reject_info,\
                                               chainmode,\
                                               stream);\
            if (memstore != NULL)\
            {\
              gt_assert(diagband_struct != NULL &&\
                        gt_diagband_struct_empty(diagband_struct));\
              segment_length = memstore->nextfreeGtDiagbandseedMaximalmatch;\
              gt_diagband_struct_multi_update(\
                      diagband_struct,\
                      memstore->spaceGtDiagbandseedMaximalmatch,\
                      segment_length);\
            } else\
            {\
              gt_assert(diagband_struct == NULL);\
            }\
          }\
          if (!seedpairlist->maxmat_compute || memstore != NULL)\
          {\
            gt_assert(segment_proc_func != NULL && segment_proc_info != NULL);\
            segment_proc_func(segment_proc_info,\
                              aencseq,\
                              bencseq,\
                              currsegm_aseqnum,\
                              currsegm_bseqnum,\
                              diagband_struct,\
                              memstore == NULL\
                                ? NULL \
                                : memstore->spaceGtDiagbandseedMaximalmatch,\
                              seedlength,\
                              segment_positions,\
                              segment_length);\
          }\
        }\
        if (diagband_struct != NULL && \
            !gt_diagband_struct_empty(diagband_struct))\
        {\
          gt_diagband_struct_reset(diagband_struct,\
                                   memstore == NULL ? segment_positions\
                                                        : NULL,\
                                   memstore == NULL\
                                     ? NULL\
                                     : memstore->\
                                           spaceGtDiagbandseedMaximalmatch,\
                                   segment_length);\
        }

static void gt_diagbandseed_match_header(FILE *stream,
                                         const GtDiagbandseedExtendParams *extp,
                                         const void *processinfo,
                                         unsigned int spacedseedweight,
                                         unsigned int seedlength,
                                         GtUword num_diagbands,
                                         GtUword minsegmentlen)
{
  fprintf(stream,"# start processing of seeds ...\n");
  fprintf(stream,"# parameters for selecting seeds: ");
  if (spacedseedweight < seedlength)
  {
    fprintf(stream,"weight=%u, span=%u for spaced seeds,",spacedseedweight,
                                                          seedlength);
  } else
  {
    fprintf(stream,"seedlength=%u,",seedlength);
  }
  fprintf(stream," diagonal bands=" GT_WU ", minimal segmentsize=" GT_WU
                 ", minimal coverage=" GT_WU "\n",
                 num_diagbands,minsegmentlen,extp->mincoverage);
  if (extp->extendgreedy)
  {
    const GtGreedyextendmatchinfo *ggemi
      = (GtGreedyextendmatchinfo *) processinfo;

    fprintf(stream,"# parameters for greedy extension of seeds: history=" GT_WU
                   ", max_aligned_length_difference=" GT_WU
                   ", percent_match_history=" GT_WU "\n",
                    extp->history_size,
                    gt_greedy_extend_maxalignedlendifference(ggemi),
                    gt_greedy_extend_perc_mat_history(ggemi));
  } else
  {
    const GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) processinfo;
    fprintf(stream,"# parameters for xdrop extension of seeds: "
                   "xdrop_below_score=" GT_WU "\n",
                   gt_xdrop_extend_belowscore(xdropmatchinfo));
  }
}

static void gt_diagbandseed_info_qm_set(
                                   GtProcessinfo_and_querymatchspaceptr *ifqm,
                                   const GtDiagbandseedExtendParams *extp,
                                   GtQuerymatchoutoptions *querymoutopt,
                                   GtReadmode query_readmode,
                                   FILE *stream,
                                   const GtKarlinAltschulStat
                                     *karlin_altschul_stat,
                                   void *processinfo)
{
  ifqm->processinfo = processinfo;
  if (extp->ani_accumulate != NULL)
  {
    ifqm->querymatchspaceptr = NULL;
  } else
  {
    ifqm->querymatchspaceptr = gt_querymatch_new();
    if (extp->verify_alignment)
    {
      gt_querymatch_verify_alignment_set(ifqm->querymatchspaceptr);
    }
    if (querymoutopt != NULL) {
      gt_querymatch_outoptions_set(ifqm->querymatchspaceptr,querymoutopt);
    }
    gt_querymatch_query_readmode_set(ifqm->querymatchspaceptr,query_readmode);
    gt_querymatch_file_set(ifqm->querymatchspaceptr, stream);
  }
  ifqm->karlin_altschul_stat = karlin_altschul_stat;
  ifqm->out_display_flag = extp->out_display_flag;;
  ifqm->previous_match_a_start = 0;
  ifqm->previous_match_a_end = 0;
  ifqm->previous_match_b_start = 0;
  ifqm->previous_match_b_end = 0;
  ifqm->previous_match_distance = 0;
  ifqm->previous_match_mismatches = 0;
}

#define GT_USEC2SEC(TIME_IN_USEC)\
         ((GtUword) (TIME_IN_USEC)/1000000)

#define GT_USECREMAIN(TIME_IN_USEC) ((TIME_IN_USEC) -\
                                     GT_USEC2SEC(TIME_IN_USEC) * 1000000UL)

#ifndef _WIN32
static void gt_diagbandseed_process_seeds_times(
                 bool extendgreedy,
                 bool maxmat_show,
                 GtUword mlistlen,
                 GtUword extended_seeds,
                 GtUword total_process_seeds_usec,
                 GtUword total_extension_time_usec)
{
  GtUword process_seeds_usec;

  gt_assert(total_process_seeds_usec >= total_extension_time_usec);
  process_seeds_usec = total_process_seeds_usec - total_extension_time_usec;
  if (!maxmat_show)
  {
    printf("# ... %s extension of " GT_WU " seeds in "
                   GT_WD ".%.06ld seconds.\n",
            extendgreedy ? "greedy" : "xdrop",
            extended_seeds,
            GT_USEC2SEC(total_extension_time_usec),
            GT_USECREMAIN(total_extension_time_usec));
  }
  printf( "# ... processed " GT_WU " seeds %sin "
                  GT_WD ".%06ld seconds.\n",mlistlen,
          maxmat_show ? "" : "(excluding extensions) ",
          GT_USEC2SEC(process_seeds_usec),
          GT_USECREMAIN(process_seeds_usec));
}
#endif

static GtDiagbandseedExtendSegmentInfo *gt_diagbandseed_extendSI_new(
                                         const GtDiagbandseedExtendParams *extp,
                                         void *processinfo,
                                         GtQuerymatchoutoptions *querymoutopt,
                                         const GtEncseq *aencseq,
                                         const GtSequencePartsInfo *aseqranges,
                                         GtUword aidx,
                                         const GtEncseq *bencseq,
                                         const GtSequencePartsInfo *bseqranges,
                                         GtUword bidx,
                                         const GtKarlinAltschulStat
                                           *karlin_altschul_stat,
                                         GtReadmode query_readmode,
                                         FILE *stream,
                                         GtDiagbandseedState
                                           *dbs_state,
                                         GtSegmentRejectFunc
                                           segment_reject_func,
                                         GtSegmentRejectInfo
                                           *segment_reject_info)
{
  GtDiagbandseedExtendSegmentInfo *esi = gt_malloc(sizeof *esi);

  esi->extend_relative_coords_function = extp->extendgreedy
                                            ? gt_greedy_extend_seed_relative
                                            : gt_xdrop_extend_seed_relative;
  gt_diagbandseed_plainsequence_init(&esi->plainsequence_info,
                                     gt_querymatch_subjectid_display(
                                             extp->out_display_flag),
                                     gt_querymatch_queryid_display(
                                             extp->out_display_flag),
                                     aencseq,
                                     aseqranges,
                                     aidx,
                                     extp->a_extend_char_access ==
                                         GT_EXTEND_CHAR_ACCESS_DIRECT ? true
                                                                      : false,
                                     bencseq,
                                     bseqranges,
                                     bidx,
                                     extp->b_extend_char_access ==
                                        GT_EXTEND_CHAR_ACCESS_DIRECT ? true
                                                                     : false);
  gt_diagbandseed_info_qm_set(&esi->info_querymatch,
                              extp,
                              querymoutopt,
                              query_readmode,
                              stream,
                              karlin_altschul_stat,
                              processinfo);
  esi->dbs_state = dbs_state;
  /* the following are constant and depends only on the given parameters */
  esi->userdefinedleastlength = extp->userdefinedleastlength;
  esi->errorpercentage = extp->errorpercentage;
  esi->use_apos = extp->use_apos;
  esi->mincoverage = extp->mincoverage;
  esi->only_selected_seqpairs = extp->only_selected_seqpairs;
  esi->evalue_threshold = extp->evalue_threshold;
  esi->query_readmode = query_readmode;
  esi->same_encseq = (aencseq == bencseq) ? true : false;
  esi->debug = gt_log_enabled() ? true : false;
  esi->segment_reject_func = segment_reject_func;
  esi->segment_reject_info = segment_reject_info;
  esi->karlin_altschul_stat = karlin_altschul_stat;
  esi->out_display_flag = extp->out_display_flag;
  esi->benchmark = extp->benchmark;
  if (extp->ani_accumulate != NULL)
  {
    if (GT_ISDIRREVERSE(query_readmode))
    {
      esi->ani_accumulate = extp->ani_accumulate + 1;
    } else
    {
      esi->ani_accumulate = extp->ani_accumulate;
    }
  } else
  {
    esi->ani_accumulate = NULL;
  }
  return esi;
}

static void gt_diagbandseed_extendSI_delete(
                                        GtDiagbandseedExtendSegmentInfo * esi)
{
  if (esi != NULL)
  {
    gt_querymatch_delete(esi->info_querymatch.querymatchspaceptr);
    gt_diagbandseed_plainsequence_delete(&esi->plainsequence_info);
    gt_free(esi);
  }
}

typedef void (*GtDiagbandseedProcessSegmentFunc)(
                        void *v_process_segment_info,
                        const GtEncseq *aencseq,
                        const GtEncseq *bencseq,
                        GtUword aseqnum,
                        GtUword bseqnum,
                        const GtDiagbandStruct *diagband_struct,
                        const GtDiagbandseedMaximalmatch *memstore,
                        unsigned int seedlength,
                        const GtSeedpairPositions *segment_positions,
                        GtUword segment_length);

/* start seed extension for seeds in mlist */
static void gt_diagbandseed_process_seeds(GtSeedpairlist *seedpairlist,
                                         const GtDiagbandseedExtendParams *extp,
                                          void *processinfo,
                                          GtQuerymatchoutoptions *querymoutopt,
                                          const GtEncseq *aencseq,
                                          const GtSequencePartsInfo *aseqranges,
                                          GtUword aidx,
                                          const GtEncseq *bencseq,
                                          const GtSequencePartsInfo *bseqranges,
                                          GtUword bidx,
                                          const GtKarlinAltschulStat
                                            *karlin_altschul_stat,
                                          GtArrayGtDiagbandseedMaximalmatch
                                            *memstore,
                                          const GtChain2Dimmode *chainmode,
                                          unsigned int spacedseedweight,
                                          unsigned int seedlength,
                                          GtReadmode query_readmode,
                                          bool verbose,
                                          FILE *stream,
                                          const GtStr *diagband_statistics_arg,
                                          GtDiagbandseedState
                                            *dbs_state,
                                          GtSegmentRejectFunc
                                            segment_reject_func,
                                          GtSegmentRejectInfo
                                            *segment_reject_info)
{
  const bool forward = query_readmode == GT_READMODE_REVCOMPL ? false : true;
  /* Although the sequences of the parts processed are shorter, we need to
     set amaxlen and bmaxlen to the maximum size of all sequences
     to get the same division into diagonal bands for all parts and thus
     obtain results independent of the number of parts chosen. */
  const GtUword mlistlen = gt_seedpairlist_length(seedpairlist),
                minsegmentlen = (extp->mincoverage - 1) / seedlength + 1;
  GtTimer *timer = NULL;
  GtDiagbandStruct *diagband_struct = NULL;
  GtDiagbandseedExtendSegmentInfo *esi = NULL;
  GtDiagbandStatistics *diagband_statistics = NULL;
  GtDiagbandseedProcessSegmentFunc segment_proc_func = NULL;
  void *segment_proc_info = NULL;

  gt_assert(extp->mincoverage >= seedlength && minsegmentlen >= 1);
  if (mlistlen == 0 || mlistlen < minsegmentlen ||
      (!extp->extendgreedy && !extp->extendxdrop))
  {
    return;
  }
  if (verbose)
  {
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  if (seedpairlist->maxmat_show)
  {
    fprintf(stream,"# Fields: s.len, s.seqnum, s.start, strand, q.seqnum, "
                   "q.start\n");
  } else
  {
    const GtUword bmaxlen = gt_encseq_max_seq_length(bencseq);
    if (verbose)
    {
      gt_diagbandseed_match_header(stream,extp,processinfo,
                                   spacedseedweight,
                                   seedlength,
                                   gt_diagband_struct_num_diagbands(
                                              seedpairlist->amaxlen,bmaxlen,
                                              extp->logdiagbandwidth),
                                   minsegmentlen);
    }
    diagband_struct = gt_diagband_struct_new(seedpairlist->amaxlen,bmaxlen,
                                             extp->logdiagbandwidth);
    if (gt_str_length(diagband_statistics_arg) == 0)
    {
      esi = gt_diagbandseed_extendSI_new(extp,
                                         processinfo,
                                         querymoutopt,
                                         aencseq,
                                         aseqranges,
                                         aidx,
                                         bencseq,
                                         bseqranges,
                                         bidx,
                                         karlin_altschul_stat,
                                         query_readmode,
                                         stream,
                                         dbs_state,
                                         segment_reject_func,
                                         segment_reject_info);
      if (verbose)
      {
        if (esi->plainsequence_info.a_byte_sequence != NULL ||
            esi->plainsequence_info.b_byte_sequence != NULL)
        {
          fprintf(stream, "# ... extracted sequences ");
          gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
          gt_timer_start(timer);
        }
      }
      segment_proc_func = gt_diagbandseed_segment2matches;
      segment_proc_info = esi;
    } else
    {
      diagband_statistics = gt_diagband_statistics_new(diagband_statistics_arg,
                                                       forward);
      segment_proc_func = gt_diagband_statistics_add;
      segment_proc_info = diagband_statistics;
    }
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_STRUCT)
  {
    const GtDiagbandseedSeedPair
      *mlist = gt_seedpairlist_mlist_struct(seedpairlist),
      *mlistend = mlist + mlistlen,
      *last_segment_start = mlistend - minsegmentlen,
      *nextsegm = mlist;

    gt_assert(nextsegm <= last_segment_start);
    /* iterate through segments of equal k-mers */
    while (nextsegm <= last_segment_start)
    {
      GtSeedpairPositions *spp_ptr, *segment_positions;
      GtUword segment_length;
      const GtDiagbandseedSeedPair *currsegm = nextsegm;
      const GtDiagbandseedSeqnum currsegm_aseqnum = currsegm->aseqnum;
      const GtDiagbandseedSeqnum currsegm_bseqnum = currsegm->bseqnum;

      /* if insuffienct number of kmers in segment: skip whole segment */
      if (currsegm_aseqnum != currsegm[minsegmentlen - 1].aseqnum ||
          currsegm_bseqnum != currsegm[minsegmentlen - 1].bseqnum)
      {
        do
        {
          nextsegm++;
        } while (nextsegm < mlistend &&
                 currsegm_aseqnum == nextsegm->aseqnum &&
                 currsegm_bseqnum == nextsegm->bseqnum);
        continue; /* process next segment */
      }

      /* this segment begining with nextsegm possibly has enough seeds */
      spp_ptr = segment_positions = (GtSeedpairPositions *) currsegm;
      do
      {
        if (!seedpairlist->maxmat_compute)
        {
          gt_diagband_struct_single_update(diagband_struct,
                                           GT_DIAGBANDSEED_GETPOS_A(nextsegm),
                                           GT_DIAGBANDSEED_GETPOS_B(nextsegm),
                                           (GtDiagbandseedPosition) seedlength);
        }
        spp_ptr->apos = GT_DIAGBANDSEED_GETPOS_A(nextsegm);
        spp_ptr->bpos = GT_DIAGBANDSEED_GETPOS_B(nextsegm);
        spp_ptr++;
        nextsegm++;
      } while (nextsegm < mlistend &&
               currsegm_aseqnum == nextsegm->aseqnum &&
               currsegm_bseqnum == nextsegm->bseqnum);

      if (esi != NULL)
      {
        gt_diagbandseed_plainsequence_next_segment(&esi->plainsequence_info,
                                                   aseqranges,
                                                   currsegm_aseqnum,
                                                   bseqranges,
                                                   currsegm_bseqnum);
      }

      /* from here on we only need the apos and bpos values of the segment, as
         the segment boundaries have been identified.
         second scan: test for mincoverage and overlap to previous extension,
         based on apos and bpos values. */
      segment_length = (GtUword) (nextsegm - currsegm);
      GT_DIAGBANDSEED_PROCESS_SEGMENT;
    }
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
    {
      const GtUword *mlist = gt_seedpairlist_mlist_ulong(seedpairlist),
                    *mlistend = mlist + mlistlen,
                    *last_segment_start = mlistend - minsegmentlen,
                    *nextsegm = mlist;
     GtUword nextsegm_a_bseqnum, segment_length;

      gt_assert(nextsegm <= last_segment_start);
      /* iterate through segments of equal k-mers */
      nextsegm_a_bseqnum
        = gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*mlist);
      while (nextsegm <= last_segment_start)
      {
        GtSeedpairPositions *spp_ptr, *segment_positions;
        GtUword currsegm_aseqnum, currsegm_bseqnum;
        const GtUword *currsegm = nextsegm;
        const GtUword currsegm_a_bseqnum = nextsegm_a_bseqnum;

        /* if insuffienct number of kmers in segment: skip whole segment */
        if (currsegm_a_bseqnum !=
            gt_seedpairlist_a_bseqnum_ulong (seedpairlist,
                                             currsegm[minsegmentlen-1]))
        {
          do
          {
            nextsegm++;
          } while (nextsegm < mlistend &&
                   currsegm_a_bseqnum ==
                   (nextsegm_a_bseqnum =
                   gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*nextsegm)));
          continue; /* process next segment */
        }

        /* this segment begining with nextsegm pssibly has enough seeds */
        currsegm_aseqnum = gt_seedpairlist_extract_ulong(seedpairlist,*currsegm,
                                                         idx_aseqnum) +
                           seedpairlist->aseqrange_start;
        currsegm_bseqnum = gt_seedpairlist_extract_ulong(seedpairlist,*currsegm,
                                                         idx_bseqnum) +
                           seedpairlist->bseqrange_start;
        spp_ptr = segment_positions = (GtSeedpairPositions *) currsegm;
        do
        {
          GtUword apos = gt_seedpairlist_extract_ulong(seedpairlist,*nextsegm,
                                                       idx_apos);

          spp_ptr->bpos
            = gt_seedpairlist_extract_ulong(seedpairlist,*nextsegm,idx_bpos);
          spp_ptr->apos = apos;
          if (!seedpairlist->maxmat_compute)
          {
            gt_diagband_struct_single_update(diagband_struct,
                                             spp_ptr->apos,
                                             spp_ptr->bpos,
                                             (GtDiagbandseedPosition)
                                                 seedlength);
          }
          spp_ptr++;
          nextsegm++;
        } while (nextsegm < mlistend &&
                 currsegm_a_bseqnum ==
                 (nextsegm_a_bseqnum =
                 gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*nextsegm)));

        if (esi != NULL)
        {
          gt_diagbandseed_plainsequence_next_segment(&esi->plainsequence_info,
                                                     aseqranges,
                                                     currsegm_aseqnum,
                                                     bseqranges,
                                                     currsegm_bseqnum);
        }

        /* from here on we only need the apos and bpos values of the segment, as
           the segment boundaries have been identified.
           second scan: test for mincoverage and overlap to previous extension,
           based on apos and bpos values. */
        segment_length = (GtUword) (nextsegm - currsegm);
        GT_DIAGBANDSEED_PROCESS_SEGMENT;
      }
    } else
    {
      GtDiagbandseedSeedPair nextsegment;
      const GtUword minsegmentlen_offset = (minsegmentlen - 1) *
                                           seedpairlist->bytes_seedpair,
                    last_segment_offset = (mlistlen - minsegmentlen) *
                                          seedpairlist->bytes_seedpair,
                    mlistlen_offset = mlistlen * seedpairlist->bytes_seedpair;
      GtUword nextsegment_offset = 0, segment_length;

      gt_assert(seedpairlist->splt == GT_DIAGBANDSEED_BASE_LIST_BYTESTRING);
      /* iterate through segments of equal k-mers, segment has length > 0 */
      gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,0);
      while (nextsegment_offset <= last_segment_offset)
      {
        GtSeedpairPositions *spp_ptr, *segment_positions;
        GtDiagbandseedSeedPair endminsegment;
        GtDiagbandseedSeqnum currsegm_aseqnum = nextsegment.aseqnum;
        GtDiagbandseedSeqnum currsegm_bseqnum = nextsegment.bseqnum;

        gt_diagbandseed_decode_seedpair(&endminsegment,seedpairlist,
                                        nextsegment_offset +
                                        minsegmentlen_offset);
        if (currsegm_aseqnum != endminsegment.aseqnum ||
            currsegm_bseqnum != endminsegment.bseqnum)
        {
          /* insuffienct number of kmers in segment: skip whole segment */
          while (true)
          {
            nextsegment_offset += seedpairlist->bytes_seedpair;
            if (nextsegment_offset >= mlistlen_offset)
            {
              break;
            }
            gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,
                                            nextsegment_offset);
            if (currsegm_aseqnum != nextsegment.aseqnum ||
                currsegm_bseqnum != nextsegment.bseqnum)
            {
              break;
            }
          }
          continue; /* process next segment */
        }

        spp_ptr = segment_positions
                = (GtSeedpairPositions *)
                  (gt_seedpairlist_mlist_bytestring(seedpairlist) +
                   nextsegment_offset);
        do
        {
          if (!seedpairlist->maxmat_compute)
          {
            gt_diagband_struct_single_update(
                                         diagband_struct,
                                         GT_DIAGBANDSEED_GETPOS_A(&nextsegment),
                                         GT_DIAGBANDSEED_GETPOS_B(&nextsegment),
                                         (GtDiagbandseedPosition) seedlength);
          }
          spp_ptr->apos = GT_DIAGBANDSEED_GETPOS_A(&nextsegment);
          spp_ptr->bpos = GT_DIAGBANDSEED_GETPOS_B(&nextsegment);
          spp_ptr++;
          nextsegment_offset += seedpairlist->bytes_seedpair;
          if (nextsegment_offset >= mlistlen_offset)
          {
            break;
          }
          gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,
                                          nextsegment_offset);
        } while (currsegm_aseqnum == nextsegment.aseqnum &&
                 currsegm_bseqnum == nextsegment.bseqnum);

        /* from here on we only need the apos and bpos values of the segment, as
           the segment boundaries have been identified.
           second scan: test for mincoverage and overlap to previous extension,
           based on apos and bpos values. */
        currsegm_aseqnum += seedpairlist->aseqrange_start;
        currsegm_bseqnum += seedpairlist->bseqrange_start;
        if (esi != NULL)
        {
          gt_diagbandseed_plainsequence_next_segment(&esi->plainsequence_info,
                                                     aseqranges,
                                                     currsegm_aseqnum,
                                                     bseqranges,
                                                     currsegm_bseqnum);
        }
        segment_length = (GtUword) (spp_ptr - segment_positions);
        GT_DIAGBANDSEED_PROCESS_SEGMENT;
      }
    }
  }
  if (diagband_struct != NULL)
  {
    if (verbose)
    {
      gt_diagband_struct_reset_counts(diagband_struct,stream);
    }
    gt_diagband_struct_delete(diagband_struct);
  }
  if (verbose)
  {
    if (dbs_state != NULL)
    {
      const GtUword numseqpairs = (seedpairlist->aseqrange_end -
                                   seedpairlist->aseqrange_start + 1) *
                                  (seedpairlist->bseqrange_end -
                                   seedpairlist->bseqrange_start + 1);
      gt_diagbandseed_dbs_state_update(dbs_state,mlistlen,numseqpairs);
#ifndef _WIN32
      dbs_state->total_process_seeds_usec
        += gt_timer_elapsed_usec(timer);
#endif
    }
    gt_timer_delete(timer);
  }
  gt_diagbandseed_extendSI_delete(esi);
  if (diagband_statistics != NULL)
  {
    gt_diagband_statistics_display(diagband_statistics);
    gt_diagband_statistics_delete(diagband_statistics);
  }
}

/* * * * * ALGORITHM STEPS * * * * */

static char *gt_diagbandseed_kmer_filename(const GtEncseq *encseq,
                                           unsigned int spacedseedweight,
                                           unsigned int seedlength,
                                           bool forward,
                                           unsigned int numparts,
                                           unsigned int partindex,
                                           GtDiagbandseedBaseListType kmplt)
{
  char *filename;
  GtStr *str = gt_str_new_cstr(gt_encseq_indexname(encseq));
  if (spacedseedweight < seedlength)
  {
    gt_str_append_char(str, '.');
    gt_str_append_uint(str, spacedseedweight);
  }
  gt_str_append_char(str, '.');
  gt_str_append_uint(str, seedlength);
  gt_str_append_char(str, forward ? 'f' : 'r');
  gt_str_append_uint(str, numparts);
  gt_str_append_char(str, '-');
  gt_str_append_uint(str, partindex + 1);
  if (kmplt == GT_DIAGBANDSEED_BASE_LIST_ULONG)
  {
    gt_str_append_char(str, 'U');
  } else
  {
    gt_assert(kmplt == GT_DIAGBANDSEED_BASE_LIST_STRUCT);
  }
  gt_str_append_cstr(str, ".kmer");
  filename = gt_cstr_dup(gt_str_get(str));
  gt_str_delete(str);
  return filename;
}

static GtUword gt_numof_kmerpos_in_file(const char *filename,
                                        GtDiagbandseedBaseListType kmplt)
{
  off_t file_size = gt_file_size(filename);
  size_t base_type_size = kmplt == GT_DIAGBANDSEED_BASE_LIST_STRUCT
                            ? sizeof (GtDiagbandseedKmerPos)
                            : sizeof (GtUword);

  gt_assert(file_size >= sizeof (GtLongestCodeRunType) &&
            (file_size - sizeof (GtLongestCodeRunType)) % base_type_size == 0);
  return (GtUword) ((file_size - sizeof (GtLongestCodeRunType))/base_type_size);
}

static GtDiagbandseedBaseListType gt_diagbandseed_kmplt(
            const GtKmerPosListEncodeInfo *encode_info)
{
  return encode_info == NULL ? GT_DIAGBANDSEED_BASE_LIST_STRUCT
                             : GT_DIAGBANDSEED_BASE_LIST_ULONG;
}

/* Go through the different steps of the seed and extend algorithm. */
static int gt_diagbandseed_algorithm(const GtDiagbandseedInfo *arg,
                                     const GtKmerPosList *alist,
                                     FILE *stream,
                                     const GtEncseq *aencseq,
                                     const GtSequencePartsInfo *aseqranges,
                                     GtUword aidx,
                                     const GtEncseq *bencseq,
                                     const GtSequencePartsInfo *bseqranges,
                                     GtUword bidx,
                                     const GtKarlinAltschulStat
                                       *karlin_altschul_stat,
                                     GtDiagbandseedState
                                       *dbs_state,
                                     GtFtTrimstat *trimstat,
                                     GtError *err)
{
  GtKmerPosList *blist = NULL;
  GtSeedpairlist *seedpairlist = NULL;
  GtDiagbandseedKmerIterator *aiter = NULL, *biter = NULL;
  GtUword alen = 0, blen = 0, mlistlen = 0, maxfreq, len_used;
  GtRange seedpairdistance;
  char *blist_file = NULL;
  int had_err = 0;
  bool alist_blist_id, both_strands, selfcomp, equalranges, use_blist = false;
  size_t sizeofunit;
  const GtDiagbandseedExtendParams *extp = NULL;
  GtFtPolishing_info *pol_info = NULL;
  void *processinfo = NULL;
  GtQuerymatchoutoptions *querymoutopt = NULL;
  GtSegmentRejectInfo *segment_reject_info = NULL;
  GtSegmentRejectFunc segment_reject_func = NULL;
  const GtUword anumseqranges = gt_sequence_parts_info_number(aseqranges),
                bnumseqranges = gt_sequence_parts_info_number(bseqranges),
                amaxlen = gt_encseq_max_seq_length(aencseq);
  GtArrayGtDiagbandseedMaximalmatch *memstore = NULL;
  GtChain2Dimmode *chainmode = NULL;
  GtKmerPosListEncodeInfo *aencode_info, *bencode_info;

  gt_assert(arg != NULL);
  aencode_info = gt_kmerpos_encode_info_new(arg->kmplt,
                                            arg->aencseq,
                                            arg->spacedseedweight,
                                            aseqranges,
                                            aidx);
  if (arg->aencseq != arg->bencseq || aseqranges != bseqranges || aidx != bidx)
  {
    bencode_info = gt_kmerpos_encode_info_new(arg->kmplt,
                                              arg->bencseq,
                                              arg->spacedseedweight,
                                              bseqranges,
                                              bidx);
  } else
  {
    bencode_info = aencode_info;
  }
  seedpairdistance = *arg->seedpairdistance;
  maxfreq = arg->maxfreq;
  selfcomp = (arg->bencseq == arg->aencseq &&
              gt_sequence_parts_info_overlap(aseqranges,aidx,bseqranges,bidx))
              ? true : false;
  equalranges = gt_sequence_parts_info_equal(aseqranges,aidx,bseqranges,bidx);
  alist_blist_id = (selfcomp && !arg->nofwd && equalranges) ? true : false;
  both_strands = (arg->norev || arg->nofwd) ? false : true;
  if (!alist_blist_id) {
    seedpairdistance.start = 0UL;
  }
  if (arg->verbose && (anumseqranges > 1 || bnumseqranges > 1))
  {
    fprintf(stream, "# process part " GT_WU " (sequences " GT_WU "..." GT_WU
                    ") vs part " GT_WU " (sequences " GT_WU "..." GT_WU ")\n",
            aidx + 1,
            gt_sequence_parts_info_start_get(aseqranges,aidx),
            gt_sequence_parts_info_end_get(aseqranges,aidx),
            bidx + 1,
            gt_sequence_parts_info_start_get(bseqranges,bidx),
            gt_sequence_parts_info_end_get(bseqranges,bidx));
  }
  extp = arg->extp;
  if (gt_querymatch_fstperquery_display(extp->out_display_flag))
  {
    segment_reject_func = gt_segment_reject_check;
    segment_reject_info
      = gt_segment_reject_info_new(
               gt_sequence_parts_info_start_get(bseqranges,bidx),
               gt_sequence_parts_info_numofsequences_get(bseqranges,bidx));
  }

  /* Create k-mer iterator for alist */
  if (alist == NULL) {
    char *alist_file
      = gt_diagbandseed_kmer_filename(arg->aencseq,
                                      arg->spacedseedweight,
                                      arg->seedlength,
                                      true,
                                      anumseqranges,
                                      aidx,
                                      gt_diagbandseed_kmplt(aencode_info));
    FILE *alist_fp = gt_fa_fopen(alist_file, "rb", err);
    if (alist_fp == NULL) {
      gt_kmerpos_encode_info_delete(aencode_info);
      gt_kmerpos_encode_info_delete(bencode_info);
      return -1;
    }
    alen = gt_numof_kmerpos_in_file(alist_file,
                                    gt_diagbandseed_kmplt(aencode_info));
    aiter = gt_diagbandseed_kmer_iter_new_file(alist_fp,aencode_info);
    gt_free(alist_file);
    alist_file = NULL;
  } else {
    gt_assert(alist != NULL);
    alen = alist->nextfree;
    aiter = gt_diagbandseed_kmer_iter_new_list(alist);
  }

  /* Second k-mer list */
  if (alist_blist_id && alist != NULL) {
    biter = gt_diagbandseed_kmer_iter_new_list(alist);
    blen = alen;
  } else if (arg->use_kmerfile) {
    blist_file
      = gt_diagbandseed_kmer_filename(arg->bencseq,
                                      arg->spacedseedweight,
                                      arg->seedlength,
                                      !arg->nofwd,
                                      bnumseqranges,
                                      bidx,
                                      gt_diagbandseed_kmplt(bencode_info));
    if (!gt_file_exists(blist_file)) {
      gt_free(blist_file);
      blist_file = NULL;
    }
  }
  if (blist_file != NULL) {
    FILE *blist_fp = gt_fa_fopen(blist_file, "rb", err);
    if (blist_fp == NULL) {
      gt_free(blist_file);
      gt_diagbandseed_kmer_iter_delete(aiter);
      aiter = NULL;
      gt_kmerpos_encode_info_delete(aencode_info);
      gt_kmerpos_encode_info_delete(bencode_info);
      return -1;
    }
    blen = gt_numof_kmerpos_in_file(blist_file,
                                    gt_diagbandseed_kmplt(bencode_info));
    gt_assert(biter == NULL);
    biter = gt_diagbandseed_kmer_iter_new_file(blist_fp,bencode_info);
    gt_free(blist_file);
    blist_file = NULL;
  } else if (!alist_blist_id) {
    const GtReadmode readmode_kmerscan = arg->nofwd ? GT_READMODE_COMPL
                                                    : GT_READMODE_FORWARD;
    const GtUword known_size = (selfcomp && equalranges) ? alen : 0;
    blist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->spacedseedweight,
                              arg->seedlength,
                              arg->spaced_seed_spec,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
                              bencode_info,
                              arg->debug_kmer,
                              arg->verbose,
                              known_size,
                              stream);
    blen = blist->nextfree;
    biter = gt_diagbandseed_kmer_iter_new_list(blist);
    use_blist = true;
  }

  len_used = alen;
  if (!selfcomp || !arg->norev) {
    len_used += blen;
  }
  seedpairlist = gt_seedpairlist_new(arg->splt,aseqranges,aidx,bseqranges,bidx,
                                     arg->maxmat,amaxlen);
  sizeofunit = gt_seedpairlist_sizeofunit(seedpairlist);
  if (seedpairlist->maxmat_compute && !seedpairlist->maxmat_show)
  {
    memstore = gt_malloc(sizeof *memstore);
    GT_INITARRAY(memstore,GtDiagbandseedMaximalmatch);
    if (gt_str_length(arg->chainarguments) > 0)
    {
      chainmode = gt_chain_chainmode_new(GT_UWORD_MAX,
                                         false,
                                         NULL,
                                         true,
                                         gt_str_get(arg->chainarguments),
                                         err);
      if (chainmode == NULL)
      {
        had_err = -1;
      }
    }
  }
  if (!had_err && !seedpairlist->maxmat_show && arg->memlimit < GT_UWORD_MAX)
  {
    had_err = gt_diagbandseed_get_mlistlen_maxfreq(&mlistlen,
                                                   &maxfreq,
                                                   aiter,
                                                   biter,
                                                   arg->memlimit,
                                                   sizeofunit,
                                                   &seedpairdistance,
                                                   len_used,
                                                   selfcomp,
                                                   alist_blist_id,
                                                   arg->verbose,
                                                   stream,
                                                   err);
  }
  if (!had_err) {
    gt_diagbandseed_kmer_iter_reset(aiter);
    gt_diagbandseed_kmer_iter_reset(biter);
    if (arg->verbose)
    {
      gt_seedpairlist_show_bits(stream,seedpairlist);
    }
    gt_diagbandseed_get_seedpairs(seedpairlist,
                                  aiter,
                                  biter,
                                  maxfreq,
                                  mlistlen,
                                  &seedpairdistance,
                                  selfcomp,
                                  arg->debug_seedpair,
                                  arg->verbose,
                                  stream);
    mlistlen = gt_seedpairlist_length(seedpairlist);
    if (arg->verify && mlistlen > 0) {
      had_err = gt_diagbandseed_verify(seedpairlist,
                                       arg->aencseq,
                                       arg->bencseq,
                                       arg->seedlength,
                                       arg->nofwd,
                                       arg->verbose,
                                       stream,
                                       err);
      if (had_err) {
        gt_seedpairlist_delete(seedpairlist);
        seedpairlist = NULL;
      }
    }
  } else
  {
    gt_seedpairlist_delete(seedpairlist);
    seedpairlist = NULL;
  }
  if (use_blist) {
    gt_kmerpos_list_delete(blist);
  }
  use_blist = false;
  gt_diagbandseed_kmer_iter_delete(biter);
  biter = NULL;

  /* Create extension info objects */
  if (!had_err)
  {
    if (extp->extendgreedy) {
      GtGreedyextendmatchinfo *grextinfo = NULL;
      const double weak_errorperc = (double)(extp->weakends
                                             ? GT_MAX(extp->errorpercentage, 20)
                                             : extp->errorpercentage);

      pol_info = polishing_info_new_with_bias(weak_errorperc,
                                              extp->matchscore_bias,
                                              extp->history_size);
      grextinfo = gt_greedy_extend_matchinfo_new(extp->maxalignedlendifference,
                                                 extp->history_size,
                                                 extp->perc_mat_history,
                                                 extp->userdefinedleastlength,
                                                 extp->errorpercentage,
                                                 extp->evalue_threshold,
                                                 extp->a_extend_char_access,
                                                 extp->b_extend_char_access,
                                                 extp->cam_generic,
                                                 extp->sensitivity,
                                                 pol_info);
      if (trimstat != NULL)
      {
        gt_greedy_extend_matchinfo_trimstat_set(grextinfo,trimstat);
      }
      processinfo = (void *) grextinfo;
    } else if (extp->extendxdrop) {
      GtXdropmatchinfo *xdropinfo = NULL;
      gt_assert(extp->extendgreedy == false);
      xdropinfo = gt_xdrop_matchinfo_new(extp->userdefinedleastlength,
                                         extp->errorpercentage,
                                         extp->evalue_threshold,
                                         extp->xdropbelowscore,
                                         extp->sensitivity);
      processinfo = (void *) xdropinfo;
    }
    if (extp->extendxdrop || extp->verify_alignment ||
        gt_querymatch_alignment_display(extp->out_display_flag) ||
        gt_querymatch_trace_display(extp->out_display_flag) ||
        gt_querymatch_dtrace_display(extp->out_display_flag) ||
        gt_querymatch_cigar_display(extp->out_display_flag) ||
        gt_querymatch_cigarX_display(extp->out_display_flag))
    {
      querymoutopt = gt_querymatchoutoptions_new(extp->out_display_flag,
                                                 NULL,
                                                 NULL);
      gt_assert(querymoutopt != NULL);
      if (extp->extendxdrop || extp->extendgreedy) {
        const GtUword sensitivity = extp->extendxdrop ? 100UL
                                                      : extp->sensitivity;
        gt_querymatchoutoptions_extend(querymoutopt,
                                       extp->errorpercentage,
                                       extp->evalue_threshold,
                                       extp->maxalignedlendifference,
                                       extp->history_size,
                                       extp->perc_mat_history,
                                       extp->a_extend_char_access,
                                       extp->b_extend_char_access,
                                       extp->cam_generic,
                                       extp->weakends,
                                       sensitivity,
                                       extp->matchscore_bias,
                                       extp->always_polished_ends,
                                       extp->out_display_flag);
      }
    }
    /* process first mlist */
    gt_assert(seedpairlist != NULL);
    gt_diagbandseed_process_seeds(seedpairlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  aencseq,aseqranges,aidx,
                                  bencseq,bseqranges,bidx,
                                  karlin_altschul_stat,
                                  memstore,
                                  chainmode,
                                  arg->spacedseedweight,
                                  arg->seedlength,
                                  arg->nofwd ? GT_READMODE_REVCOMPL
                                             : GT_READMODE_FORWARD,
                                  arg->verbose,
                                  stream,
                                  arg->diagband_statistics_arg,
                                  dbs_state,
                                  segment_reject_func,
                                  segment_reject_info);
    gt_seedpairlist_reset(seedpairlist);
    gt_querymatchoutoptions_reset(querymoutopt);

    /* Third (reverse) k-mer list */
    if (both_strands) {
      GtUword mrevlen = 0;
      GtKmerPosList *clist = NULL;

      gt_assert(blist_file == NULL && !use_blist);
      seedpairdistance.start = 0UL;
      if (arg->use_kmerfile) {
        blist_file
          = gt_diagbandseed_kmer_filename(arg->bencseq,
                                          arg->spacedseedweight,
                                          arg->seedlength,
                                          false,
                                          bnumseqranges,
                                          bidx,
                                          gt_diagbandseed_kmplt(bencode_info));
        if (!gt_file_exists(blist_file)) {
          gt_free(blist_file);
          blist_file = NULL;
        }
      }
      if (blist_file != NULL) {
        FILE *blist_fp = gt_fa_fopen(blist_file, "rb", err);
        if (blist_fp == NULL) {
          had_err = -1;
        } else {
          biter = gt_diagbandseed_kmer_iter_new_file(blist_fp,bencode_info);
        }
        gt_free(blist_file);
      } else {
        const GtReadmode readmode_kmerscan = GT_READMODE_COMPL;
        clist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->spacedseedweight,
                              arg->seedlength,
                              arg->spaced_seed_spec,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
                              bencode_info,
                              arg->debug_kmer,
                              arg->verbose,
                              blen,
                              stream);
        biter = gt_diagbandseed_kmer_iter_new_list(clist);
        use_blist = true;
      }

      if (!had_err) {
        gt_diagbandseed_kmer_iter_reset(aiter);
        if (!seedpairlist->maxmat_show && arg->memlimit < GT_UWORD_MAX)
        {
          had_err = gt_diagbandseed_get_mlistlen_maxfreq(&mrevlen,
                                                         &maxfreq,
                                                         aiter,
                                                         biter,
                                                         arg->memlimit,
                                                         sizeofunit,
                                                         &seedpairdistance,
                                                         len_used,
                                                         selfcomp,
                                                         alist_blist_id,
                                                         arg->verbose,
                                                         stream,
                                                         err);
        }
      }

      if (!had_err) {
        gt_diagbandseed_kmer_iter_reset(aiter);
        gt_diagbandseed_kmer_iter_reset(biter);
        gt_diagbandseed_get_seedpairs(seedpairlist,
                                      aiter,
                                      biter,
                                      maxfreq,
                                      mrevlen,
                                      &seedpairdistance,
                                      selfcomp,
                                      arg->debug_seedpair,
                                      arg->verbose,
                                      stream);
        mrevlen = gt_seedpairlist_length(seedpairlist);
        if (arg->verify && mrevlen > 0) {
          had_err = gt_diagbandseed_verify(seedpairlist,
                                           arg->aencseq,
                                           arg->bencseq,
                                           arg->seedlength,
                                           true,
                                           arg->verbose,
                                           stream,
                                           err);
          if (had_err) {
            gt_seedpairlist_delete(seedpairlist);
          }
        }
      }
      if (use_blist) {
        gt_kmerpos_list_delete(clist);
      }
      gt_diagbandseed_kmer_iter_delete(biter);
      biter = NULL;
    }
  }
  gt_diagbandseed_kmer_iter_delete(aiter);
  aiter = NULL;

  /* Process second (reverse) mlist */
  if (!had_err && both_strands) {
    gt_diagbandseed_process_seeds(seedpairlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  aencseq,aseqranges,aidx,
                                  bencseq,bseqranges,bidx,
                                  karlin_altschul_stat,
                                  memstore,
                                  chainmode,
                                  arg->spacedseedweight,
                                  arg->seedlength,
                                  GT_READMODE_REVCOMPL,
                                  arg->verbose,
                                  stream,
                                  arg->diagband_statistics_arg,
                                  dbs_state,
                                  segment_reject_func,
                                  segment_reject_info);
  }
  /* Clean up */
  gt_seedpairlist_delete(seedpairlist);
  if (memstore != NULL)
  {
    GT_FREEARRAY(memstore,GtDiagbandseedMaximalmatch);
    if (chainmode != NULL)
    {
      gt_chain_chainmode_delete(chainmode);
    }
    gt_free(memstore);
  }
  if (extp->extendgreedy)
  {
    polishing_info_delete(pol_info);
    gt_greedy_extend_matchinfo_delete((GtGreedyextendmatchinfo *) processinfo);
  } else
  {
    if (extp->extendxdrop)
    {
      gt_xdrop_matchinfo_delete((GtXdropmatchinfo *) processinfo);
    }
  }
  gt_querymatchoutoptions_delete(querymoutopt);
  if (segment_reject_info != NULL)
  {
    gt_segment_reject_info_delete(segment_reject_info);
  }
  gt_kmerpos_encode_info_delete(aencode_info);
  if (aencode_info != bencode_info)
  {
    gt_kmerpos_encode_info_delete(bencode_info);
  }
  return had_err;
}

#ifdef GT_THREADS_ENABLED
typedef struct
{
  const GtDiagbandseedInfo *arg;
  const GtKmerPosList *alist;
  FILE *stream;
  const GtEncseq *aencseq, *bencseq;
  const GtSequencePartsInfo *aseqranges,
                            *bseqranges;
  GtSegmentRejectFunc segment_reject_func;
  GtArray *combinations;
  int had_err;
  GtError *err;
  const GtKarlinAltschulStat *karlin_altschul_stat;
} GtDiagbandseedThreadInfo;

static void gt_diagbandseed_thread_info_set(GtDiagbandseedThreadInfo *ti,
                                     const GtDiagbandseedInfo *arg,
                                     const GtKmerPosList *alist,
                                     FILE *stream,
                                     const GtEncseq *aencseq,
                                     const GtSequencePartsInfo *aseqranges,
                                     const GtEncseq *bencseq,
                                     const GtSequencePartsInfo *bseqranges,
                                     const GtKarlinAltschulStat
                                       *karlin_altschul_stat,
                                     GtArray *combinations,
                                     GtError *err)
{
  gt_assert(ti != NULL);
  ti->arg = arg;
  ti->alist = alist;
  ti->stream = stream;
  ti->aencseq = aencseq;
  ti->aseqranges = aseqranges;
  ti->bencseq = bencseq;
  ti->bseqranges = bseqranges;
  ti->karlin_altschul_stat = karlin_altschul_stat;
  ti->combinations = gt_array_clone(combinations);
  ti->had_err = 0;
  ti->err = err;
}

static void *gt_diagbandseed_thread_algorithm(void *thread_info)
{
  GtDiagbandseedThreadInfo *info = (GtDiagbandseedThreadInfo *)thread_info;
  if (gt_array_size(info->combinations) > 0) {
    const GtUwordPair *last = gt_array_get_last(info->combinations),
                      *comb;

    for (comb = gt_array_get_first(info->combinations); comb <= last; comb++) {
      info->had_err = gt_diagbandseed_algorithm(
                           info->arg,
                           info->alist,
                           info->stream,
                           info->aencseq,
                           info->aseqranges,
                           comb->a,
                           info->bencseq,
                           info->bseqranges,
                           comb->b,
                           info->karlin_altschul_stat,
                           NULL,
                           NULL,
                           info->err);
      if (info->had_err) break;
    }
  }
  gt_array_delete(info->combinations);
  return NULL;
}
#endif

static int gt_diagbandseed_write_kmers(const GtKmerPosList *kmerpos_list,
                                       const char *path,
                                       unsigned int spacedseedweight,
                                       unsigned int seedlength,
                                       bool verbose,
                                       GtError *err)
{
  FILE *stream;
  GtUword longest_code_run;

  if (verbose) {
    gt_assert(kmerpos_list != NULL);
    printf("# write " GT_WU " %u-mers ",kmerpos_list->nextfree, seedlength);
    if (spacedseedweight < seedlength)
    {
      printf("with weight %u ",spacedseedweight);
    }
    printf("to file %s\n",path);
  }

  stream = gt_fa_fopen(path, "wb", err);
  if (stream == NULL)
  {
    return -1;
  }
  gt_xfwrite(&kmerpos_list->longest_code_run,sizeof longest_code_run,1,stream);
  if (kmerpos_list->encode_info != NULL)
  {
    gt_xfwrite(kmerpos_list->spaceGtUword,
               sizeof *kmerpos_list->spaceGtUword,
               kmerpos_list->nextfree, stream);
  } else
  {
    gt_xfwrite(kmerpos_list->spaceGtDiagbandseedKmerPos,
               sizeof *kmerpos_list->spaceGtDiagbandseedKmerPos,
               kmerpos_list->nextfree, stream);
  }
  gt_fa_fclose(stream);
  return 0;
}

static bool gt_create_or_update_file(const char *path,const GtEncseq *encseq)
{
  if (gt_file_exists(path))
  {
    GtStr *esqfile = gt_str_new_cstr(gt_encseq_indexname(encseq));

    gt_str_append_cstr(esqfile, ".esq");
    if (gt_file_is_newer(path,gt_str_get(esqfile)))
    {
      gt_str_delete(esqfile);
      return false;
    }
    gt_str_delete(esqfile);
  }
  return true;
}

static void gt_diagbandseed_out_sequences_with_matches(
                char seqtype,
                const GtEncseq *encseq,
                const GtBittab *used_sequences)
{
  GtUword seqnum, idx, numbits, max_seq_length;
  char *string_buffer;

  gt_assert(encseq != NULL && used_sequences != NULL);
  numbits = gt_bittab_count_set_bits(used_sequences);
  max_seq_length = gt_encseq_max_seq_length(encseq);
  string_buffer = gt_malloc(sizeof *string_buffer * max_seq_length);
  for (idx = 0, seqnum = gt_bittab_get_first_bitnum(used_sequences);
       idx < numbits;
       idx++, seqnum = gt_bittab_get_next_bitnum(used_sequences,seqnum))
  {
    const GtUword seqstartpos = gt_encseq_seqstartpos(encseq,seqnum),
                  seqlength = gt_encseq_seqlength(encseq,seqnum);
    gt_encseq_extract_decoded(encseq,string_buffer,seqstartpos,
                              seqstartpos + seqlength - 1);
    printf("S\t%c" GT_WU "\t" GT_WU "\t",seqtype,seqnum,seqlength);
    fwrite(string_buffer,sizeof *string_buffer,seqlength,stdout);
    fputc('\n',stdout);
  }
  gt_free(string_buffer);
}

/* Run the algorithm by iterating over all combinations of sequence ranges. */
int gt_diagbandseed_run(const GtDiagbandseedInfo *arg,
                        const GtSequencePartsInfo *aseqranges,
                        const GtSequencePartsInfo *bseqranges,
                        const GtUwordPair *pick,
                        GtError *err)
{
  const bool self = arg->aencseq == arg->bencseq ? true : false;
  const bool apick = pick->a != GT_UWORD_MAX ? true : false;
  const bool bpick = pick->b != GT_UWORD_MAX ? true : false;
  GtKmerPosList *alist = NULL;
  GtUword aidx, bidx;
  int had_err = 0;
  const GtUword anumseqranges = gt_sequence_parts_info_number(aseqranges),
                bnumseqranges = gt_sequence_parts_info_number(bseqranges);
  GtFtTrimstat *trimstat = NULL;
  GtKarlinAltschulStat *karlin_altschul_stat = NULL;
  GtDiagbandseedState *dbs_state = NULL;
#ifdef GT_THREADS_ENABLED
  GtDiagbandseedThreadInfo *tinfo = gt_malloc(gt_jobs * sizeof *tinfo);
  FILE **stream_tab;
  unsigned int tidx;

  /* create output streams */
  stream_tab = gt_malloc(gt_jobs * sizeof *stream_tab);
  stream_tab[0] = stdout;
  for (tidx = 1; !had_err && tidx < gt_jobs; tidx++) {
    stream_tab[tidx]
      = gt_xtmpfp_generic(NULL, GT_TMPFP_OPENBINARY | GT_TMPFP_AUTOREMOVE);
  }
#endif
  if (arg->verbose || gt_querymatch_gfa2_display(arg->extp->out_display_flag))
  {
    GtUword a_num_sequences = 0, b_num_sequences = 0;
    if (gt_querymatch_gfa2_display(arg->extp->out_display_flag))
    {
      a_num_sequences = gt_encseq_num_of_sequences(arg->aencseq);
      if (!self)
      {
        b_num_sequences = gt_encseq_num_of_sequences(arg->bencseq);
      }
    }
    dbs_state = gt_diagbandseed_dbs_state_new(arg->verbose,
                                              a_num_sequences,
                                              b_num_sequences);
  }

  /* create all missing k-mer lists for bencseq */
  if (arg->use_kmerfile) {
    unsigned int count;
    for (count = 0; count < 2; count++) {
      const bool fwd = count == 0 ? true : false;

      if ((fwd && (self || arg->nofwd)) || (!fwd && arg->norev))
      {
        continue;
      }
      for (bidx = 0; !had_err && bidx < bnumseqranges; bidx++)
      {
        char *path;
        GtKmerPosListEncodeInfo *bencode_info;

        if (bpick && pick->b != bidx)
        {
          continue;
        }
        bencode_info = gt_kmerpos_encode_info_new(arg->kmplt,
                                                  arg->bencseq,
                                                  arg->spacedseedweight,
                                                  bseqranges,
                                                  bidx);
        path = gt_diagbandseed_kmer_filename(arg->bencseq,
                                             arg->spacedseedweight,
                                             arg->seedlength,
                                             fwd,
                                             bnumseqranges,
                                             bidx,
                                             gt_diagbandseed_kmplt(
                                               bencode_info));
        if (gt_create_or_update_file(path,arg->bencseq))
        {
          GtKmerPosList *blist;
          GtReadmode readmode_kmerscan = fwd ? GT_READMODE_FORWARD
                                             : GT_READMODE_COMPL;
          blist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->spacedseedweight,
                              arg->seedlength,
                              arg->spaced_seed_spec,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
                              bencode_info,
                              arg->debug_kmer,
                              arg->verbose,
                              0,
                              stdout);
          had_err = gt_diagbandseed_write_kmers(blist, path,
                                                arg->spacedseedweight,
                                                arg->seedlength,
                                                arg->verbose, err);
          gt_kmerpos_list_delete(blist);
        }
        gt_free(path);
        gt_kmerpos_encode_info_delete(bencode_info);
      }
    }
  }
  if (!had_err && arg->trimstat_on)
  {
    trimstat = gt_ft_trimstat_new();
  }
  if (gt_querymatch_evalue_display(arg->extp->out_display_flag) ||
      gt_querymatch_bitscore_display(arg->extp->out_display_flag) ||
      arg->extp->evalue_threshold != DBL_MAX)
  {
    GtTimer *timer = NULL;

    if (arg->verbose)
    {
      timer = gt_timer_new();
      gt_timer_start(timer);
    }
    karlin_altschul_stat = gt_karlin_altschul_stat_new_gapped(
                                      gt_encseq_total_length(arg->aencseq),
                                      gt_encseq_num_of_sequences(arg->aencseq),
                                      arg->bencseq);
    if (arg->verbose)
    {
      gt_timer_show_formatted(timer,
                              "# ... computed lookup table for E-values "
                              GT_DIAGBANDSEED_FMT,stdout);
      gt_timer_delete(timer);
    }
  }
  for (aidx = 0; !had_err && aidx < anumseqranges; aidx++) {
    /* create alist here to prevent redundant calculations */
    char *path = NULL;
    bool use_alist = false;
    GtKmerPosListEncodeInfo *aencode_info;

    if (apick && pick->a != aidx)
    {
      continue;
    }
    aencode_info = gt_kmerpos_encode_info_new(arg->kmplt,
                                              arg->aencseq,
                                              arg->spacedseedweight,
                                              aseqranges,
                                              aidx);
    if (arg->use_kmerfile) {
      path = gt_diagbandseed_kmer_filename(arg->aencseq,
                                           arg->spacedseedweight,
                                           arg->seedlength,
                                           true,
                                           anumseqranges,
                                           aidx,
                                           gt_diagbandseed_kmplt(
                                              aencode_info));
    }
    if (!arg->use_kmerfile || gt_create_or_update_file(path,arg->aencseq))
    {
      use_alist = true;
      alist = gt_diagbandseed_get_kmers(
                              arg->aencseq,
                              arg->spacedseedweight,
                              arg->seedlength,
                              arg->spaced_seed_spec,
                              GT_READMODE_FORWARD,
                              gt_sequence_parts_info_start_get(aseqranges,aidx),
                              gt_sequence_parts_info_end_get(aseqranges,aidx),
                              aencode_info,
                              arg->debug_kmer,
                              arg->verbose,
                              0,
                              stdout);
      if (arg->use_kmerfile)
      {
        had_err = gt_diagbandseed_write_kmers(alist, path,
                                              arg->spacedseedweight,
                                              arg->seedlength,
                                              arg->verbose, err);
      }
    }
    if (arg->use_kmerfile) {
      gt_free(path);
    }
    bidx = self ? aidx : 0;

#ifdef GT_THREADS_ENABLED
    if (gt_jobs <= 1) {
#endif
      while (!had_err && bidx < bnumseqranges) {
        if (!bpick || pick->b == bidx) {
          /* start algorithm with chosen sequence ranges */
          had_err = gt_diagbandseed_algorithm(
                           arg,
                           use_alist ? alist : NULL,
                           stdout,
                           arg->aencseq,aseqranges,aidx,
                           arg->bencseq,bseqranges,bidx,
                           karlin_altschul_stat,
                           dbs_state,
                           trimstat,
                           err);
        }
        bidx++;
      }
#ifdef GT_THREADS_ENABLED
    } else if (!arg->use_kmerfile) {
      const GtUword num_runs = bpick ? 1 : bnumseqranges - bidx;
      const GtUword num_runs_per_thread = (num_runs - 1) / gt_jobs + 1;
      const GtUword num_threads = (num_runs - 1) / num_runs_per_thread + 1;
      GtArray *combinations = gt_array_new(sizeof (GtUwordPair));
      GtArray *threads = gt_array_new(sizeof (GtThread *));

      gt_assert(bidx < bnumseqranges);
      gt_assert(num_threads <= gt_jobs);
      gt_assert(!bpick || num_threads == 1);

      /* start additional threads */
      for (tidx = 1; !had_err && tidx < num_threads; tidx++) {
        GtThread *thread;
        GtUword idx;
        bidx += num_runs_per_thread;
        const GtUword end = GT_MIN(bidx + num_runs_per_thread, bnumseqranges);

        for (idx = bidx; idx < end; idx++) {
          GtUwordPair comb = {aidx, idx};
          gt_array_add(combinations, comb);
        }
        gt_diagbandseed_thread_info_set(tinfo + tidx,
                                        arg,
                                        use_alist ? alist : NULL,
                                        stream_tab[tidx],
                                        arg->aencseq,
                                        aseqranges,
                                        arg->bencseq,
                                        bseqranges,
                                        karlin_altschul_stat,
                                        combinations,
                                        err);
        gt_array_reset(combinations);
        if ((thread = gt_thread_new(gt_diagbandseed_thread_algorithm,
                                    tinfo + tidx, err)) != NULL) {
          gt_array_add(threads, thread);
        } else {
          had_err = -1;
        }
      }

      /* start main thread */
      if (!had_err) {
        GtUword idx;
        bidx = self ? aidx : 0;
        for (idx = bidx;
             idx < GT_MIN(bidx + num_runs_per_thread, bnumseqranges);
             ++idx) {
          if (!bpick || pick->b == idx) {
            GtUwordPair comb = {aidx, idx};
            gt_array_add(combinations, comb);
          }
        }
        gt_diagbandseed_thread_info_set(tinfo,
                                        arg,
                                        use_alist ? alist : NULL,
                                        stream_tab[0],
                                        arg->aencseq,
                                        aseqranges,
                                        arg->bencseq,
                                        bseqranges,
                                        karlin_altschul_stat,
                                        combinations,
                                        err);
        gt_diagbandseed_thread_algorithm(tinfo);
      }

      /* clean up */
      for (tidx = 0; tidx < gt_array_size(threads); tidx++) {
        GtThread *thread = *(GtThread**) gt_array_get(threads, tidx);
        if (!had_err) {
          gt_thread_join(thread);
        }
        gt_thread_delete(thread);
      }
      gt_array_delete(threads);
      for (tidx = 0; tidx < num_threads && !had_err; tidx++) {
        had_err = tinfo[tidx].had_err;
      }
      gt_array_delete(combinations);
    }
#endif
    if (use_alist) {
      gt_kmerpos_list_delete(alist);
    }
    gt_kmerpos_encode_info_delete(aencode_info);
  }
#ifdef GT_THREADS_ENABLED
  if (gt_jobs > 1 && arg->use_kmerfile) {
    GtArray *combinations[gt_jobs];
    GtArray *threads = gt_array_new(sizeof (GtThread *));
    GtUword counter = 0;
    for (tidx = 0; tidx < gt_jobs; tidx++) {
      combinations[tidx] = gt_array_new(sizeof (GtUwordPair));
    }
    for (aidx = 0; aidx < anumseqranges; aidx++) {
      if (apick && pick->a != aidx)
      {
        continue;
      }
      for (bidx = self ? aidx : 0; bidx < bnumseqranges; bidx++) {
        if (!bpick || pick->b == bidx) {
          GtUwordPair comb = {aidx, bidx};
          gt_array_add(combinations[counter++ % gt_jobs], comb);
        }
      }
    }

    for (tidx = 1; !had_err && tidx < gt_jobs; tidx++) {
      GtThread *thread;
      gt_diagbandseed_thread_info_set(tinfo + tidx,
                                      arg,
                                      NULL,
                                      stream_tab[tidx],
                                      arg->aencseq,
                                      aseqranges,
                                      arg->bencseq,
                                      bseqranges,
                                      karlin_altschul_stat,
                                      combinations[tidx],
                                      err);
      if ((thread = gt_thread_new(gt_diagbandseed_thread_algorithm,
                                  tinfo + tidx, err)) != NULL) {
        gt_array_add(threads, thread);
      } else {
        had_err = -1;
      }
    }
    /* start main thread */
    if (!had_err) {
      gt_diagbandseed_thread_info_set(tinfo,
                                      arg,
                                      NULL,
                                      stream_tab[0],
                                      arg->aencseq,
                                      aseqranges,
                                      arg->bencseq,
                                      bseqranges,
                                      karlin_altschul_stat,
                                      combinations[0],
                                      err);
      gt_diagbandseed_thread_algorithm(tinfo);
    }

    /* clean up */
    for (tidx = 0; tidx < gt_array_size(threads); tidx++) {
      GtThread *thread = *(GtThread**) gt_array_get(threads, tidx);
      if (!had_err) {
        gt_thread_join(thread);
      }
      gt_thread_delete(thread);
    }
    for (tidx = 0; tidx < gt_array_size(threads) && !had_err; tidx++) {
      had_err = tinfo[tidx].had_err;
    }
    gt_array_delete(threads);
    for (tidx = 0; tidx < gt_jobs; tidx++) {
      gt_array_delete(combinations[tidx]);
    }
  }
  gt_free(tinfo);

  /* print the threads' output to stdout */
  for (tidx = 1; tidx < gt_jobs; tidx++) {
    int cc;
    rewind(stream_tab[tidx]);
    while ((cc = fgetc(stream_tab[tidx])) != EOF) {
      putchar(cc);
    }
    gt_fa_xfclose(stream_tab[tidx]);
  }
  gt_free(stream_tab);

#endif
  if (arg->verbose)
  {
    const bool maxmat_show = gt_diagbandseed_derive_maxmat_show(arg->maxmat);

    gt_assert(dbs_state != NULL);
    gt_diagbandseed_dbs_state_out(dbs_state,maxmat_show);
#ifndef _WIN32
    gt_diagbandseed_process_seeds_times(
                 arg->extp->extendgreedy,
                 maxmat_show,
                 dbs_state->totalseeds,
                 dbs_state->extended_seeds,
                 dbs_state->total_process_seeds_usec,
                 dbs_state->total_extension_time_usec);
#endif
  }
  if (dbs_state != NULL)
  {
    if (dbs_state->used_a_sequences != NULL)
    {
      gt_diagbandseed_out_sequences_with_matches('S',arg->aencseq,
                  dbs_state->used_a_sequences);
    }
    if (dbs_state->used_b_sequences != NULL)
    {
      gt_diagbandseed_out_sequences_with_matches('Q',arg->bencseq,
                  dbs_state->used_b_sequences);
    }
  }
  gt_diagbandseed_dbs_state_delete(dbs_state);
  gt_karlin_altschul_stat_delete(karlin_altschul_stat);
  gt_ft_trimstat_delete(trimstat);
  return had_err;
}
