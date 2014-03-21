/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>

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

#ifndef S_SPLINT_S
#include <errno.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <sam.h>

#include "core/fa.h"
#include "core/bittab_api.h"
#include "core/chardef.h"
#include "core/compat.h"
#include "core/disc_distri_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/parseutils_api.h"
#include "core/queue_api.h"
#include "core/safearith.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/bitinstream.h"
#include "extended/bitoutstream.h"
#include "extended/elias_gamma.h"
#include "extended/encdesc.h"
#include "extended/golomb.h"
#include "extended/huffcode.h"
#include "extended/sam_alignment.h"
#include "extended/sam_query_name_iterator.h"
#include "extended/samfile_iterator.h"
#include "extended/rcr.h"

#define BAMBASEA 1
#define BAMBASEC 2
#define BAMBASEG 4
#define BAMBASET 8
#define BAMBASEN 15

#define PHREDOFFSET 33

#ifndef PAGESIZE
#define PAGESIZE gt_pagesize()
#endif
#define NUMOFPAGES2MAP 100

#define RCR_NUMOFSEPS 8UL
#define MAX_NUM_VAL_HUF 512UL

#define ENDOFRECORD 9

#define RCR_LINEWIDTH 80UL
#define DESCSEPSEQ '@'
#define DESCSEPQUAL '+'

#define DEFAULTMQUAL 0
#define DEFAULTQUAL '-'

/* TODO: use ONE struct for both, this is duplicating code and stupid */
struct GtRcrEncoder {
  FILE               *output,
                     *unmapped_reads_ptr;
  GtBitOutStream     *bitstream;
  GtCstrIterator     *cstr_iterator;
  GtDiscDistri       *readlength_distr,
                     *readpos_distr,
                     *varpos_distr,
                     *qual_distr,
                     *qual_mapping_distr;
  GtEncdescEncoder   *encdesc_enc;
  GtGolomb           *readpos_golomb,
                     *varpos_golomb;
  GtHuffman          *readlenghts_huff,
                     *qual_huff,
                     *qual_mapping_huff,
                     *cigar_ops_huff,
                     *bases_huff;
  GtSamfileIterator  *sam_iter;
  bam1_t             *sam_align;
  const GtEncseq     *encseq;
  const char         *samfilename;
  GtUint64 *ins_bases;
  GtQueue            *not_exact_matches;
  GtUint64  all_bits,
                      dellen_bits,
                      encodedbases,
                      exact_match_flag_bits,
                      ins_bases_bits,
                      mapqual_bits,
                      pos_bits,
                      present_cigar_ops[ENDOFRECORD + 1],
                      qual_bits,
                      readlen_bits,
                      sclip_bits,
                      skiplen_bits,
                      strand_bits,
                      subs_bits,
                      varpos_bits,
                      vartype_bits;
  GtUword       cur_read,
                      cur_seq_startpos,
                      max_read_length,
                      numofreads,
                      numofunmappedreads,
                      prev_readpos,
                      readlength;
  bool                cons_readlength,
                      is_num_fields_cons,
                      is_verbose,
                      store_all_qual,
                      store_mapping_qual,
                      store_unmmaped_reads,
                      store_var_qual;
};

struct GtRcrDecoder {
  FILE               *fp;
  GtEncdesc          *encdesc;
  GtGolomb           *readpos_golomb,
                     *varpos_golomb;
  GtHuffman          *readlenghts_huff,
                     *qual_huff,
                     *qual_mapping_huff,
                     *cigar_ops_huff,
                     *bases_huff;
  GtStr              *inputname;
  const GtEncseq     *encseq;
  const char         *basename;
  GtUint64 *ins_bases;
  GtUint64  present_cigar_ops[ENDOFRECORD + 1];
  GtUword       numofreads,
                      cur_bit,
                      cur_bitseq,
                      readlength,
                      startofencoding;
  bool                cons_readlength,
                      store_all_qual,
                      store_var_qual,
                      store_mapping_qual,
                      is_num_fields_cons;
};

typedef struct MedianData {
  GtUint64 n,
                     x;
  GtUword      median;
} MedianData;

static inline char rcr_bambase2char(uint8_t base)
{
  switch (base) {
    case BAMBASEA:
      return 'A';
    case BAMBASEC:
      return 'C';
    case BAMBASEG:
      return 'G';
    case BAMBASET:
      return 'T';
    case BAMBASEN:
    default:
      return 'N';
  }
}

static inline GtUchar rcr_bambase2gtbase(uint8_t base, GtAlphabet *alpha)
{
  switch (base) {
    case BAMBASEN:
      return (GtUchar) WILDCARD;
    default:
      return gt_alphabet_encode(alpha, rcr_bambase2char(base));
  }
}

static GtBitsequence rcr_transencode(GtUchar ref, GtUchar base,
                                     GtAlphabet *alpha)
{
  GtUchar code_a = gt_alphabet_encode(alpha, 'A'),
          code_c = gt_alphabet_encode(alpha, 'C'),
          code_g = gt_alphabet_encode(alpha, 'G'),
          code_t = gt_alphabet_encode(alpha, 'T');

  if (base == (GtUchar) WILDCARD)
    return (GtBitsequence) 3;

  if (ref == code_a) {
    if (base == code_c)
      return (GtBitsequence) 0;
    if (base == code_g)
      return (GtBitsequence) 1;
    if (base == code_t)
      return (GtBitsequence) 2;
   }
  if (ref == code_c) {
    if (base == code_a)
      return (GtBitsequence) 0;
    if (base == code_g)
      return (GtBitsequence) 1;
    if (base == code_t)
      return (GtBitsequence) 2;
  }
 if (ref == code_g) {
    if (base == code_a)
      return (GtBitsequence) 0;
    if (base == code_c)
      return (GtBitsequence) 1;
    if (base == code_t)
      return (GtBitsequence) 2;
  }
  if (ref == code_t) {
    if (base == code_a)
      return (GtBitsequence) 0;
    if (base == code_c)
      return (GtBitsequence) 1;
    if (base == code_g)
      return (GtBitsequence) 2;
  }
  if (ref == (GtUchar) WILDCARD) {
    if (base == code_a)
      return (GtBitsequence) 0;
    if (base == code_c)
      return (GtBitsequence) 1;
    if (base == code_g)
      return (GtBitsequence) 2;
    if (base == code_t)
      return (GtBitsequence) 3;
  }
  return (GtBitsequence) GT_UNDEF_UWORD;
}

static GtUchar rcr_transdecode(GtUchar ref, GtBitsequence transcode,
                               GtAlphabet *alpha)
{
  GtUchar code_a = gt_alphabet_encode(alpha, 'A'),
          code_c = gt_alphabet_encode(alpha, 'C'),
          code_g = gt_alphabet_encode(alpha, 'G'),
          code_t = gt_alphabet_encode(alpha, 'T');

  if (ref == (GtUchar) WILDCARD) {
    switch (transcode) {
       case 0:
         return code_a;
       case 1:
         return code_c;
       case 2:
         return code_g;
       case 3:
         return code_t;
    }
  }
  else if (transcode == (GtBitsequence) 3)
    return (GtUchar) WILDCARD;
  else {
    if (ref == code_a) {
      switch (transcode) {
        case 0:
          return code_c;
        case 1:
          return code_g;
        case 2:
          return code_t;
      }
    }
    if (ref == code_c) {
       switch (transcode) {
         case 0:
           return code_a;
         case 1:
           return code_g;
         case 2:
           return code_t;
         }
    }
    if (ref == code_g) {
      switch (transcode) {
        case 0:
          return code_a;
        case 1:
          return code_c;
        case 2:
          return code_t;
      }
    }
    if (ref == code_t) {
      switch (transcode) {
        case 0:
          return code_a;
        case 1:
          return code_c;
        case 2:
          return code_g;
      }
    }
  }
  return (GtUchar) GT_UNDEF_UCHAR;
}

static void rcr_convert_cigar_string(GtStr *cigar_str)
{
  GtUword i,
                cur_cigar_len = 1UL;
  char cur_cigar_op = gt_str_get(cigar_str)[0];
  GtStr *new_cigar_str = gt_str_new();

  for (i = 1UL; i < gt_str_length(cigar_str); i++) {
    if ((cur_cigar_op != gt_str_get(cigar_str)[i])) {
      gt_str_append_ulong(new_cigar_str, cur_cigar_len);
      gt_str_append_char(new_cigar_str, cur_cigar_op);
      cur_cigar_op = gt_str_get(cigar_str)[i];
      cur_cigar_len = 1UL;
      if (i == gt_str_length(cigar_str) - 1) {
        gt_str_append_ulong(new_cigar_str, cur_cigar_len);
        gt_str_append_char(new_cigar_str, cur_cigar_op);
      }
    }
    else if (i == gt_str_length(cigar_str) - 1) {
      gt_str_append_ulong(new_cigar_str, ++cur_cigar_len);
      gt_str_append_char(new_cigar_str, cur_cigar_op);
    }
    else
      cur_cigar_len++;
  }
  gt_str_reset(cigar_str);
  gt_str_append_str(cigar_str, new_cigar_str);
  gt_str_delete(new_cigar_str);
}

static void rcr_write_read_to_file(FILE *fp, uint8_t *seq, uint8_t *qual,
                                  const char *desc, GtUword seq_l)
{
  GtUword i,
                cur_width;
  gt_xfputc(DESCSEPSEQ, fp);
  gt_xfputs(desc, fp);
  gt_xfputc('\n', fp);

  for (i = 0, cur_width = 0; i < seq_l; i++, cur_width++) {
    if (cur_width == RCR_LINEWIDTH) {
      cur_width = 0;
      gt_xfputc('\n', fp);
    }
    gt_xfputc(rcr_bambase2char((uint8_t) bam1_seqi(seq, i)), fp);
  }
  gt_xfputc('\n', fp);
  gt_xfputc(DESCSEPQUAL, fp);
  gt_xfputc('\n', fp);

  for (i = 0, cur_width = 0; i < seq_l; i++, cur_width++) {
    if (cur_width == RCR_LINEWIDTH) {
      cur_width = 0;
      gt_xfputc('\n', fp);
    }
    gt_xfputc((int) qual[i] + PHREDOFFSET, fp);
  }
  gt_xfputc('\n', fp);
}

static void rcr_huff_encode_write(GtRcrEncoder *rcr_enc,
                                  GtHuffman *huff,
                                  GtUword val)
{
  GtBitsequence code;
  unsigned bits_to_write;

  gt_huffman_encode(huff, val, &code, &bits_to_write);
  rcr_enc->all_bits += bits_to_write;
  rcr_enc->vartype_bits += bits_to_write;
  gt_bitoutstream_append(rcr_enc->bitstream, code, bits_to_write);
}

static void rcr_golomb_encode_write(GtRcrEncoder *rcr_enc,
                                    GtGolomb *gol,
                                    GtUword val)
{
  GtBittab *code = gt_golomb_encode(gol, val);
  GtUword size = gt_bittab_size(code);
  rcr_enc->all_bits += size;
  rcr_enc->varpos_bits += size;
  gt_bitoutstream_append_bittab(rcr_enc->bitstream, code);
  gt_bittab_delete(code);
}

static void rcr_elias_encode_write(GtRcrEncoder *rcr_enc,
                                   GtUword val)
{
  GtBittab *code = gt_elias_gamma_encode(val);
  GtUword size = gt_bittab_size(code);
  rcr_enc->all_bits += size;
  rcr_enc->dellen_bits += size;
  gt_bitoutstream_append_bittab(rcr_enc->bitstream, code);
  gt_bittab_delete(code);
}

static void rcr_encode_write_var_type(GtRcrEncoder *rcr_enc,
                                      GtUword cigar_op)
{
  rcr_huff_encode_write(rcr_enc, rcr_enc->cigar_ops_huff, cigar_op);
}

static void rcr_encode_write_var_pos(GtRcrEncoder *rcr_enc,
                                     GtUword rel_varpos)
{
  rcr_golomb_encode_write(rcr_enc, rcr_enc->varpos_golomb, rel_varpos);
}

#define RCR_UPDATE_VAR_POS(rel, pos, prev)                                    \
  do {                                                                        \
  rel = pos - prev; prev = pos;                                               \
  } while (false)

static int rcr_write_read_encoding(const bam1_t *alignment,
                                   GtRcrEncoder *rcr_enc)
{
  int had_err = 0;
  GtUchar ref,
          base;
  GtUword cigar_op,
                cigar_len;
  unsigned bits_to_write,
           one_bit = 1U;
  GtUword alpha_size,
                i,
                j,
                prev_varpos,
                qual,
                read_i,
                readlength,
                readpos,
                ref_i,
                rel_readpos = 0,
                rel_varpos,
                varpos;
  GtBitsequence code,
                one = (GtBitsequence) 1,
                zero = 0;
  uint32_t *cigar;
  uint8_t *qual_string,
          *seq_string;
  const bam1_core_t *core;
  GtAlphabet *encseq_alpha = gt_encseq_alphabet(rcr_enc->encseq);

  alpha_size = (GtUword) gt_alphabet_size(encseq_alpha);

  core = &alignment->core;
  qual_string = bam1_qual(alignment);
  seq_string = bam1_seq(alignment);
  cigar = bam1_cigar(alignment);

  /* only one option should be set */
  gt_assert(!(rcr_enc->store_all_qual && rcr_enc->store_var_qual));

  /* read is unmapped */
  if (core->flag & BAM_FUNMAP) {
    if (rcr_enc->store_unmmaped_reads)
      rcr_write_read_to_file(rcr_enc->unmapped_reads_ptr,
                             seq_string,
                             qual_string,
                             bam1_qname(alignment),
                             (GtUword) core->l_qseq);
    gt_bitoutstream_append(rcr_enc->bitstream, one, one_bit);
    return 0;
  }
  else
    gt_bitoutstream_append(rcr_enc->bitstream, zero, one_bit);

  /* encode read length */
  if (!rcr_enc->cons_readlength) {
    readlength = (GtUword) core->l_qseq;
    rcr_huff_encode_write(rcr_enc, rcr_enc->readlenghts_huff, readlength);
  }
  else
    readlength = rcr_enc->readlength;

  rcr_enc->encodedbases += readlength;

  gt_safe_assign(readpos, core->pos);
  ref_i = readpos + rcr_enc->cur_seq_startpos;
  read_i = 0;

  /* encode relative read position */
  if (!had_err) {
    gt_assert(readpos >= rcr_enc->prev_readpos);
    gt_safe_sub(rel_readpos, readpos, rcr_enc->prev_readpos);
    rcr_enc->prev_readpos = readpos;
    rcr_golomb_encode_write(rcr_enc, rcr_enc->readpos_golomb, rel_readpos);
  }

  /* write mapping qual */
  if (rcr_enc->store_mapping_qual) {
    qual = (GtUword) core->qual;
    rcr_huff_encode_write(rcr_enc, rcr_enc->qual_mapping_huff, qual);
  }
  /* encode qual string */
  if (rcr_enc->store_all_qual) {
    for (i = 0; i < readlength; i++) {
      qual = ((GtUword) qual_string[i]) + PHREDOFFSET;
      rcr_huff_encode_write(rcr_enc, rcr_enc->qual_huff, qual);
    }
  }

  /* write strand */
  if (core->flag & BAM_FREVERSE)
    gt_bitoutstream_append(rcr_enc->bitstream, one, one_bit);
  else
    gt_bitoutstream_append(rcr_enc->bitstream, zero, one_bit);
  rcr_enc->all_bits++;
  rcr_enc->strand_bits++;

  /* exact match? */
  if (gt_queue_size(rcr_enc->not_exact_matches) > 0 &&
      (void*) rcr_enc->cur_read ==
      gt_queue_head(rcr_enc->not_exact_matches)) {
    (void) gt_queue_get(rcr_enc->not_exact_matches);
    gt_bitoutstream_append(rcr_enc->bitstream, zero, one_bit);
    rcr_enc->all_bits++;
    rcr_enc->exact_match_flag_bits++;

    prev_varpos = 0;

    for (i = 0; i < (GtUword) core->n_cigar; i++) {
      gt_safe_assign(cigar_op, (cigar[i] & 0xf));
      gt_safe_assign(cigar_len, (cigar[i]>>4));

      if (cigar_op == (GtUword) BAM_CEQUAL ||
          cigar_op == (GtUword) BAM_CDIFF)
        cigar_op = BAM_CMATCH;

      switch (cigar_op) {
        case BAM_CMATCH:
          for (j = 0; j < cigar_len; j++) {
            ref = gt_encseq_get_encoded_char(rcr_enc->encseq,
                                     (ref_i + j),
                                     GT_READMODE_FORWARD);
            base =
              rcr_bambase2gtbase((uint8_t) bam1_seqi(seq_string, read_i + j),
                                 encseq_alpha);
            if (ref != base) {
              rcr_encode_write_var_type(rcr_enc, (GtUword) cigar_op);

              /* encode variation position */
              varpos = read_i + j;

              RCR_UPDATE_VAR_POS(rel_varpos, varpos, prev_varpos);
              rcr_encode_write_var_pos(rcr_enc, rel_varpos);

              /* write transition code */
              code = rcr_transencode(ref, base, encseq_alpha);
              if (code == (GtBitsequence) GT_UNDEF_UINT)
                return -1;
              rcr_enc->all_bits += 2;
              rcr_enc->subs_bits += 2;
              bits_to_write = 2U;
              gt_bitoutstream_append(rcr_enc->bitstream, code, bits_to_write);

              if (rcr_enc->store_var_qual) {
                qual = ((GtUword) qual_string[varpos]) + PHREDOFFSET;
                rcr_huff_encode_write(rcr_enc, rcr_enc->qual_huff, qual);
              }
            }
          }
          read_i += cigar_len;
          ref_i += cigar_len;
          break;

        case BAM_CDEL:
        case BAM_CREF_SKIP:
          rcr_encode_write_var_type(rcr_enc, (GtUword) cigar_op);

          /* encode variation position */
          varpos = read_i;

          RCR_UPDATE_VAR_POS(rel_varpos, varpos, prev_varpos);
          rcr_encode_write_var_pos(rcr_enc, rel_varpos);

          /* encode length of skip/del */
          rcr_elias_encode_write(rcr_enc, cigar_len);
          ref_i += cigar_len;
          break;

        case BAM_CINS:
        case BAM_CSOFT_CLIP:
          rcr_encode_write_var_type(rcr_enc, (GtUword) cigar_op);

          /* encode varation position */
          varpos = read_i;

          RCR_UPDATE_VAR_POS(rel_varpos, varpos, prev_varpos);
          rcr_encode_write_var_pos(rcr_enc, rel_varpos);

          /* encode inserted bases */
          for (j = 0; j < cigar_len; j++) {
            base = rcr_bambase2gtbase(
                (uint8_t) bam1_seqi(seq_string, read_i + j), encseq_alpha);

            /* use alphabet_size - 1 as wildcard symbol */
            /* XXX change this, it doesn't matter if wildcard is in a continuous
               range with the rest */
            if (base == (GtUchar) WILDCARD)
              base = (GtUchar) (alpha_size - 1);

            rcr_huff_encode_write(rcr_enc, rcr_enc->bases_huff,
                                  (GtUword) base);
          }

          /* append end symbol */
          rcr_huff_encode_write(rcr_enc, rcr_enc->bases_huff, alpha_size);

          if (rcr_enc->store_var_qual) {
            for (j = 0; j < cigar_len; j++) {
              qual = ((GtUword) qual_string[read_i + j]) + PHREDOFFSET;
              rcr_huff_encode_write(rcr_enc, rcr_enc->qual_huff, qual);
            }
          }
          read_i += cigar_len;
          break;
        default:
          /* XXX gt_error nutzen */
          gt_log_log("encountered funny cigar_op: "GT_WU"", cigar_op);
          return -1;
          ;
      }
    }
    /* end symbol of a record */
    rcr_encode_write_var_type(rcr_enc, (GtUword) ENDOFRECORD);
    if (readlength != read_i) {
      /* XXX gt_error nutzen */
      gt_log_log("readlength: "GT_WU", read_i: "GT_WU"", readlength, read_i);
      return -1;
    }
  }
  else {
    gt_bitoutstream_append(rcr_enc->bitstream, one, one_bit);
    rcr_enc->all_bits++;
    rcr_enc->exact_match_flag_bits++;
  }
  rcr_enc->cur_read++;
  return 0;
}

static void rcr_get_m(GtUword key, GtUint64 value,
                      void *data)
{
  MedianData *md = (MedianData*) data;
  md->x = md->x + value;
  if (md->x > (md->n) / 2 && (md->median == GT_UNDEF_UWORD))
    md->median = key;
}

static void rcr_get_n(GT_UNUSED GtUword key, GtUint64 value,
                      void *data)
{
  MedianData *md = (MedianData*) data;
  md->n = md->n + value;
}

static GtUword rcr_get_median(GtDiscDistri *distr)
{
  GtUword median;
  MedianData *md;

  md = gt_malloc(sizeof (MedianData));
  md->n = 0;
  md->x = 0;
  md->median = GT_UNDEF_UWORD;
  gt_disc_distri_foreach(distr, rcr_get_n, md);
  gt_disc_distri_foreach(distr, rcr_get_m, md);
  median = md->median;
  gt_free(md);
  return median;
}

static GtUint64 rcr_disc_distri_func(const void *data,
                                               GtUword symbol)
{
  GtDiscDistri *distr = (GtDiscDistri*)data;
  return gt_disc_distri_get(distr, symbol);
}

static GtUint64 rcr_array_func(const void *data, GtUword symbol)
{
  GtUword *distr = (GtUword*)data;
  return (GtUint64) distr[symbol];
}

static int rcr_initialize_encoders(GtRcrEncoder *rcr_enc,
                                   GtTimer *timer,
                                   GtError *err)
{
  bool has_var = true;
  GtUword median;

  if (timer != NULL)
    gt_timer_show_progress(timer, "initializing encoders", stdout);

  median = rcr_get_median(rcr_enc->readpos_distr);
  if (median == GT_UNDEF_UWORD) {
    gt_error_set(err, "no mapped reads present in %s",rcr_enc->samfilename);
    return -1;
  }
  /* parameter m must not be zero in golomb coding */
  if (median == 0)
    median = 1UL;

  rcr_enc->readpos_golomb = gt_golomb_new(median);

  median = rcr_get_median(rcr_enc->varpos_distr);

  /* no variations present -> have only exact matches */
  if (median == GT_UNDEF_UWORD)
    has_var = false;

  if (median == 0)
    median = 1UL;

  if (has_var) {
    rcr_enc->varpos_golomb = gt_golomb_new(median);

  }

  if (!(rcr_enc->cons_readlength)) {
    rcr_enc->readlenghts_huff = gt_huffman_new(rcr_enc->readlength_distr,
                                            rcr_disc_distri_func,
                                            rcr_enc->max_read_length + 1);
  }

  if (rcr_enc->store_all_qual || rcr_enc->store_var_qual) {
    rcr_enc->qual_huff = gt_huffman_new(rcr_enc->qual_distr,
                                     rcr_disc_distri_func, 256UL);
  }

  if (rcr_enc->store_mapping_qual) {

    rcr_enc->qual_mapping_huff = gt_huffman_new(rcr_enc->qual_mapping_distr,
                                             rcr_disc_distri_func, 256UL);
  }

  rcr_enc->cigar_ops_huff = gt_huffman_new(rcr_enc->present_cigar_ops,
                                           rcr_array_func, ENDOFRECORD + 1UL);
  if (!(rcr_enc->cigar_ops_huff)) {
    gt_error_set(err, "unable to initialize GtHuffmanCoding object");
    return -1;
  }

  rcr_enc->bases_huff =
    gt_huffman_new(rcr_enc->ins_bases,
                   rcr_array_func,
                   (GtUword)
                     gt_alphabet_size(gt_encseq_alphabet(rcr_enc->encseq)) + 1);
  if (!(rcr_enc->bases_huff)) {
    gt_error_set(err, "unable to initialize GtHuffmanCoding object");
    return -1;
  }

  return 0;
}

static int rcr_get_read_infos(const bam1_t *alignment, GtRcrEncoder *rcr_enc)
{
  uint32_t *cigar;
  const bam1_core_t *bam_core;
  GtUword i,
                j,
                rel_readpos,
                readpos,
                readlength,
                varpos,
                prev_varpos,
                rel_varpos,
                q,
                read_i,
                ref_i;
  int cigar_op,
      cigar_len;
  uint8_t *qual_string,
          *seq_string;
  bool exact_match = true;
  GtUchar base,
          ref;
  bam_core = &alignment->core;
  qual_string = bam1_qual(alignment);
  seq_string = bam1_seq(alignment);
  cigar = bam1_cigar(alignment);

  /* read is unmapped */
  if (bam_core->flag & BAM_FUNMAP)
    return 0;
  if (rcr_enc->store_mapping_qual)
    gt_disc_distri_add(rcr_enc->qual_mapping_distr,
                       (GtUword) bam_core->qual);

  readlength = (GtUword) bam_core->l_qseq;
  gt_disc_distri_add(rcr_enc->readlength_distr, readlength);
  if (rcr_enc->readlength == 0) {
    rcr_enc->readlength = readlength;
    rcr_enc->max_read_length = readlength;
  }
  else {
    if (rcr_enc->readlength != readlength) {
      if (readlength > rcr_enc->max_read_length)
        rcr_enc->max_read_length = readlength;
      rcr_enc->cons_readlength = false;
    }
  }

  if (rcr_enc->store_all_qual) {
    for (j = 0; j < readlength; j++) {
      q = ((GtUword) qual_string[j]) + PHREDOFFSET;
      gt_disc_distri_add(rcr_enc->qual_distr, q);
    }
  }

  readpos = (GtUword) bam_core->pos;
  ref_i = readpos;
  read_i = 0;

  rel_readpos = ref_i - rcr_enc->prev_readpos;
  gt_disc_distri_add(rcr_enc->readpos_distr, rel_readpos);
  rcr_enc->prev_readpos = readpos;

  varpos = 0;
  prev_varpos = 0;

  for (i = 0; i < (GtUword) bam_core->n_cigar; i++) {
    GtAlphabet *alpha = gt_encseq_alphabet(rcr_enc->encseq);
    unsigned int alpha_size = gt_alphabet_size(alpha);

    gt_safe_assign(cigar_op, (cigar[i] & 0xf));
    gt_safe_assign(cigar_len, (cigar[i] >> 4));

    /* cequal and cdiff are matches in the alignment */
    if (cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF)
      cigar_op = BAM_CMATCH;
    switch (cigar_op) {
      case BAM_CMATCH:
        for (j = 0; j < (GtUword) cigar_len; j++) {
          base = rcr_bambase2gtbase((uint8_t) bam1_seqi(seq_string, read_i + j),
                                alpha);
          ref = gt_encseq_get_encoded_char(rcr_enc->encseq,
                                         ref_i + j + rcr_enc->cur_seq_startpos,
                                         GT_READMODE_FORWARD);
          if (ref != base) {
            rcr_enc->present_cigar_ops[BAM_CMATCH]++;
            exact_match = false;
            varpos = read_i + j;

            if (rcr_enc->store_var_qual) {
              q = ((GtUword) qual_string[read_i + j]) + PHREDOFFSET;
              gt_disc_distri_add(rcr_enc->qual_distr, q);
            }
          }
        }
        read_i += cigar_len;
        ref_i += cigar_len;
        break;

      case BAM_CSOFT_CLIP:
        rcr_enc->present_cigar_ops[BAM_CSOFT_CLIP]++;
        exact_match = false;
        varpos = read_i;

        for (j = 0; j < (GtUword) cigar_len; j++) {
          base = rcr_bambase2gtbase((uint8_t) bam1_seqi(seq_string, read_i + j),
                                   alpha);

          if (base ==
                gt_alphabet_encode(alpha,
                                   (char) gt_alphabet_wildcard_show(
                                     alpha)))
            base = (GtUchar) (alpha_size - 1);
          rcr_enc->ins_bases[base]++;
        }

        if (rcr_enc->store_var_qual) {
          for (j = 0; j < (GtUword) cigar_len; j++) {
            q = ((GtUword) qual_string[read_i + j]) + PHREDOFFSET;
            gt_disc_distri_add(rcr_enc->qual_distr, q);
          }
        }
        read_i += cigar_len;
        break;

      case BAM_CDEL:
        rcr_enc->present_cigar_ops[BAM_CDEL]++;
        exact_match = false;
        varpos = read_i;
        ref_i += cigar_len;
        break;

      case BAM_CREF_SKIP:
        rcr_enc->present_cigar_ops[BAM_CREF_SKIP]++;
        exact_match = false;
        varpos = read_i;
        ref_i += cigar_len;
        break;

      case BAM_CINS:
        rcr_enc->present_cigar_ops[BAM_CINS]++;
        exact_match = false;
        varpos = read_i;

        for (j = 0; j < (GtUword) cigar_len; j++) {
          base = rcr_bambase2gtbase((uint8_t) bam1_seqi(seq_string, read_i + j),
                                   alpha);

          if (base == gt_alphabet_encode(alpha,
                                         (char) gt_alphabet_wildcard_show(
                                         alpha)))
            base = (GtUchar) (alpha_size - 1);
          rcr_enc->ins_bases[base]++;
        }

        rcr_enc->ins_bases[alpha_size]++;

        if (rcr_enc->store_var_qual) {
          for (j = 0; j < (GtUword) cigar_len; j++) {
            q = ((GtUword) qual_string[read_i + j]) + PHREDOFFSET;
            gt_disc_distri_add(rcr_enc->qual_distr, q);
          }
        }
        read_i += cigar_len;
        break;

      default:
        ;
    }
    rel_varpos = varpos - prev_varpos;
    gt_disc_distri_add(rcr_enc->varpos_distr, rel_varpos);
    prev_varpos = varpos;
  }

  /* store read number of inexact matches */
  if (!exact_match)
    gt_queue_add(rcr_enc->not_exact_matches, (void*) rcr_enc->cur_read);

  gt_safe_assign(rcr_enc->prev_readpos, bam_core->pos);
  rcr_enc->cur_read = rcr_enc->cur_read + 1;
  return 0;
}

static inline GtRcrEncoder *gt_rcr_encoder_init(const char *filename,
                                                const GtEncseq *ref)
{
  unsigned alpha_size= gt_alphabet_size(gt_encseq_alphabet(ref));
  GtUword i;

  GtRcrEncoder *rcr_enc = gt_malloc(sizeof (*rcr_enc));

  rcr_enc->encseq = ref;
  rcr_enc->samfilename = filename;
  rcr_enc->not_exact_matches = gt_queue_new();
  rcr_enc->qual_distr = gt_disc_distri_new();
  rcr_enc->qual_mapping_distr = gt_disc_distri_new();
  rcr_enc->readlength_distr = gt_disc_distri_new();
  rcr_enc->readpos_distr = gt_disc_distri_new();
  rcr_enc->varpos_distr = gt_disc_distri_new();
  rcr_enc->sam_align = bam_init1();
  rcr_enc->bases_huff = NULL;
  rcr_enc->cigar_ops_huff = NULL;
  rcr_enc->encdesc_enc = NULL;
  rcr_enc->qual_huff = NULL;
  rcr_enc->qual_mapping_huff = NULL;
  rcr_enc->readlenghts_huff = NULL;
  rcr_enc->readpos_golomb = NULL;
  rcr_enc->sam_iter = NULL;
  rcr_enc->cstr_iterator = NULL;
  rcr_enc->varpos_golomb = NULL;
  rcr_enc->cons_readlength = true;
  rcr_enc->is_verbose = false;
  rcr_enc->store_mapping_qual = false;
  rcr_enc->store_var_qual = false;
  rcr_enc->store_all_qual = false;
  rcr_enc->store_unmmaped_reads = false;
  rcr_enc->cur_read = 0;
  rcr_enc->cur_seq_startpos = 0;
  rcr_enc->max_read_length = 0;
  rcr_enc->numofreads = 0;
  rcr_enc->numofunmappedreads = 0;
  rcr_enc->prev_readpos = 0;
  rcr_enc->readlength = 0;

  rcr_enc->ins_bases = gt_calloc((size_t) (alpha_size + 1),
                              sizeof (GtUint64));
  for (i = 0; i < (GtUword) (alpha_size + 1); i ++)
    rcr_enc->ins_bases[i] = 0;

  for (i = 0; i <= (GtUword) ENDOFRECORD; i++)
    rcr_enc->present_cigar_ops[i] = 0;

  return rcr_enc;
}

static inline int gt_rcr_analyse_alignment_data(GtRcrEncoder *rcr_enc,
                                                GtTimer *timer,
                                                GtError *err)
{
  int had_err = 0;
  int32_t seq_id = 0;
  samfile_t *samfile = samopen(rcr_enc->samfilename, "rb", NULL);

  gt_assert(rcr_enc->sam_align != NULL);

  if (timer != NULL)
    gt_timer_show_progress(timer, "analyse sam/bam alignment data",
                           stdout);

  if (samfile == NULL) {
    gt_error_set(err, "Cannot open BAM file %s", rcr_enc->samfilename);
    had_err = -1;
  }

  while (!had_err &&
         samread(samfile, rcr_enc->sam_align) >= 0) {
    gt_assert(rcr_enc->sam_align != NULL);
    if (seq_id != rcr_enc->sam_align->core.tid) {
      rcr_enc->prev_readpos = 0;
      seq_id = rcr_enc->sam_align->core.tid;
      rcr_enc->cur_seq_startpos = gt_encseq_seqstartpos(rcr_enc->encseq,
                                                        (GtUword) seq_id);
    }
    if (rcr_enc->prev_readpos > (GtUword) rcr_enc->sam_align->core.pos) {
      gt_error_set(err, "file %s is not sorted", rcr_enc->samfilename);
      had_err = -1;
    }
    else {
      seq_id = rcr_enc->sam_align->core.tid;
      if (!(rcr_enc->sam_align->core.flag & BAM_FUNMAP))
        rcr_enc->numofreads++;
      else
        rcr_enc->numofunmappedreads++;

      had_err = rcr_get_read_infos(rcr_enc->sam_align, rcr_enc);
    }
  }

  /* end symbol */
  rcr_enc->present_cigar_ops[ENDOFRECORD] = (GtUint64)
                                              rcr_enc->numofreads;
  samclose(samfile);

  return had_err;
}

static int gt_rcr_init_desc_encoding(GtRcrEncoder *rcr_enc,
                                     GtTimer *timer,
                                     GtError *err)
{
  int had_err = 0;
  if (timer != NULL)
    gt_timer_show_progress(timer, "init desc encoding",
                           stdout);

  rcr_enc->encdesc_enc = gt_encdesc_encoder_new();
  gt_encdesc_encoder_set_sampling_none(rcr_enc->encdesc_enc);
  rcr_enc->sam_iter =
    gt_samfile_iterator_new_bam(rcr_enc->samfilename,
                                gt_encseq_alphabet(rcr_enc->encseq),
                                err);
  if (rcr_enc->sam_iter == NULL)
    had_err = -1;
  if (!had_err) {
    rcr_enc->cstr_iterator = gt_sam_query_name_iterator_new(rcr_enc->sam_iter,
                                                            err);
    if (rcr_enc->cstr_iterator == NULL)
      had_err = -1;
  }
  return had_err;
}

GtRcrEncoder *gt_rcr_encoder_new(const GtEncseq *ref,
                                 const char *filename,
                                 bool vquals,
                                 bool mquals,
                                 bool quals,
                                 bool ureads,
                                 bool descs,
                                 GtTimer *timer,
                                 GtError *err)
{
  int had_err = 0;
  GtRcrEncoder *rcr_enc;

  gt_error_check(err);
  gt_assert(ref && filename);

  if (timer != NULL)
    gt_timer_show_progress(timer, "Initializing RcrEncoder", stdout);

  rcr_enc = gt_rcr_encoder_init(filename, ref);

  if (quals) {
    gt_assert(!vquals);
    rcr_enc->store_all_qual = true;
  }
  if (vquals) {
    gt_assert(!quals);
    rcr_enc->store_var_qual = true;
  }
  if (mquals)
    rcr_enc->store_mapping_qual = true;
  if (ureads)
    rcr_enc->store_unmmaped_reads = true;

  if (!(gt_alphabet_is_dna(gt_encseq_alphabet(ref)))) {
    gt_error_set(err, "alphabet has to be DNA");
    had_err = -1;
  }

  if (!had_err)
    had_err = gt_rcr_analyse_alignment_data(rcr_enc, timer, err);

  if (!had_err)
    had_err = rcr_initialize_encoders(rcr_enc, timer, err);

  if (!had_err && descs)
    had_err = gt_rcr_init_desc_encoding(rcr_enc, timer, err);

  if (!had_err)
    return rcr_enc;
  else {
    gt_rcr_encoder_delete(rcr_enc);
    return NULL;
  }
}

typedef struct{
  FILE *fp;
} RcrData;

static int rcr_write_node(GtUword symbol,
                          GtUint64 freq,
                          GT_UNUSED GtBitsequence code,
                          GT_UNUSED unsigned length,
                          void *pt)
{
  RcrData *data = (RcrData*) pt;
  if (fwrite(&symbol, sizeof (symbol), (size_t) 1, data->fp) == (size_t) 1)
    if (fwrite(&freq, sizeof (freq), (size_t) 1, data->fp) == (size_t) 1)
      return 0;
  return -1;
}

static int rcr_write_distr_to_file(FILE *fp, GtHuffman *huff)
{
  RcrData data;
  data.fp = fp;
  return gt_huffman_iterate(huff, rcr_write_node, &data);
}

static int rcr_write_header_to_file(GtRcrEncoder *rcr_enc)
{
  GtUword numofleaves,
                m;
  FILE *fp = rcr_enc->output;

  gt_xfwrite_one(&rcr_enc->numofreads, fp);
  gt_xfwrite_one(&rcr_enc->cons_readlength, fp);

  if (!rcr_enc->cons_readlength) {
    numofleaves = gt_huffman_numofsymbols(rcr_enc->readlenghts_huff);
    gt_xfwrite_one(&numofleaves, fp);
    gt_xfwrite_one(&rcr_enc->max_read_length, fp);
    if (rcr_write_distr_to_file(fp, rcr_enc->readlenghts_huff) != 0)
      return -1;
  }
  else {
    gt_xfwrite_one(&rcr_enc->readlength, fp);
  }

  gt_xfwrite_one(&rcr_enc->store_all_qual, fp);
  gt_xfwrite_one(&rcr_enc->store_var_qual, fp);

  if (rcr_enc->store_all_qual || rcr_enc->store_var_qual) {
    numofleaves = gt_huffman_numofsymbols(rcr_enc->qual_huff);

    gt_xfwrite_one(&numofleaves, fp);
    if (rcr_write_distr_to_file(fp, rcr_enc->qual_huff) != 0)
      return -1;
  }
  gt_xfwrite_one(&rcr_enc->store_mapping_qual, fp);

  if (rcr_enc->store_mapping_qual) {
    numofleaves = gt_huffman_numofsymbols(rcr_enc->qual_mapping_huff);
    gt_xfwrite_one(&numofleaves, fp);
    if (rcr_write_distr_to_file(fp, rcr_enc->qual_mapping_huff) != 0)
      return -1;
  }

  m = gt_golomb_get_m(rcr_enc->readpos_golomb);
  gt_xfwrite_one(&m, fp);

  if (rcr_enc->varpos_golomb != NULL)
    m = gt_golomb_get_m(rcr_enc->varpos_golomb);
  else
    m = GT_UNDEF_UWORD;

  gt_xfwrite_one(&m, fp);

  gt_xfwrite(rcr_enc->present_cigar_ops, sizeof (*rcr_enc->present_cigar_ops),
             (size_t) (ENDOFRECORD + 1), fp);

  gt_xfwrite(rcr_enc->ins_bases, sizeof (GtUint64),
             (size_t) gt_alphabet_size(gt_encseq_alphabet(rcr_enc->encseq)) + 1,
             fp);

  return 0;
}

static int rcr_write_encoding_to_file(GtRcrEncoder *rcr_enc, GtError *err)
{
  samfile_t *samfile;
  int32_t tid = (int32_t) -1;
  unsigned one_bit = 1U;
  GtBitsequence new_ref = (GtBitsequence) 1,
                old_ref = 0;

  gt_error_check(err);
  gt_assert(rcr_enc);

  rcr_enc->all_bits = 0;
  rcr_enc->qual_bits = 0;
  rcr_enc->mapqual_bits = 0;
  rcr_enc->dellen_bits = 0;
  rcr_enc->ins_bases_bits = 0;
  rcr_enc->vartype_bits = 0;
  rcr_enc->readlen_bits = 0;
  rcr_enc->pos_bits = 0;
  rcr_enc->varpos_bits = 0;
  rcr_enc->strand_bits = 0;
  rcr_enc->subs_bits = 0;
  rcr_enc->skiplen_bits = 0;
  rcr_enc->exact_match_flag_bits =0 ;
  rcr_enc->sclip_bits = 0;
  rcr_enc->encodedbases = 0;

  samfile = samopen(rcr_enc->samfilename, "rb", NULL);
  if (samfile == NULL) {
    gt_error_set(err, "Cannot open BAM file %s", rcr_enc->samfilename);
    return -1;
  }
  rcr_enc->bitstream = gt_bitoutstream_new(rcr_enc->output);

  while (samread(samfile, rcr_enc->sam_align) >= 0) {
    gt_assert(rcr_enc->sam_align != NULL);
    if (tid != rcr_enc->sam_align->core.tid) {
      tid = rcr_enc->sam_align->core.tid;
      rcr_enc->prev_readpos = 0;
      rcr_enc->cur_seq_startpos =
        gt_encseq_seqstartpos(rcr_enc->encseq, (GtUword) tid);
      gt_bitoutstream_append(rcr_enc->bitstream, new_ref, one_bit);
      gt_log_log("reset pos for new ref "GT_WU"", rcr_enc->cur_seq_startpos);
    }
    else
      gt_bitoutstream_append(rcr_enc->bitstream, old_ref, one_bit);

    if (rcr_write_read_encoding(rcr_enc->sam_align, rcr_enc) != 0) {
      return -1;
    }
  }
  gt_bitoutstream_flush(rcr_enc->bitstream);
  gt_bitoutstream_delete(rcr_enc->bitstream);
  samclose(samfile);

#ifndef S_SPLINT_S
  if (rcr_enc->is_verbose) {
    printf("encoded "GT_WU" BAM records, "GT_WU" reads(s) unmapped\n",
           rcr_enc->numofreads, rcr_enc->numofunmappedreads);
    printf("encoded "GT_LLU" bases\n",rcr_enc->encodedbases);
    printf("total number of bits used for encoding: "GT_LLU"\n",
           rcr_enc->all_bits);
    printf("%% bits used for quality values: %.3f\n",
           rcr_enc->qual_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for mapping quality values: %.3f\n",
           rcr_enc->mapqual_bits * 100.0 / rcr_enc->all_bits );
    printf("%% bits used for deletion length: %.3f\n",
           rcr_enc->dellen_bits * 100.0 / rcr_enc->all_bits );
    printf("%% bits used for inserted bases: %.3f\n",
           rcr_enc->ins_bases_bits * 100.0 / rcr_enc->all_bits );
    printf("%% bits used for variation type: %.3f\n",
           rcr_enc->vartype_bits * 100.0 / rcr_enc->all_bits );
    printf("%% bits used for read length: %.3f\n",
           rcr_enc->readlen_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for read position on reference: %.3f\n",
           rcr_enc->pos_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for variation position on read: %.3f\n",
           rcr_enc->varpos_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for strand information: %.3f\n",
           rcr_enc->strand_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for substituted bases: %.3f\n",
           rcr_enc->subs_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for skip length: %.3f\n",
           rcr_enc->skiplen_bits * 100.0 / rcr_enc->all_bits);
    printf("%% bits used for soft clipped bases: %.3f\n",
           rcr_enc->sclip_bits * 100.0 / rcr_enc->all_bits );
    printf("%% bits used for exact match flag: %.3f\n",
           rcr_enc->exact_match_flag_bits * 100.0 / rcr_enc->all_bits );
  }
#endif

  return 0;
}

static int rcr_write_data(const char *name, GtRcrEncoder *rcr_enc, GtError *err)
{
  bool is_not_at_pageborder;
  int had_err = 0;
  GtWord fpos;
  GtUword pagesize = gt_pagesize();
  GtStr *unmapped_reads_filename;
  gt_error_check(err);

  rcr_enc->output = gt_fa_fopen_with_suffix(name, RCRFILESUFFIX, "wb", err);
  if (rcr_enc->output == NULL)
    return -1;
  if (rcr_write_header_to_file(rcr_enc) != 0)
    return -1;
  fpos = ftell(rcr_enc->output);

  /* jump to beginning of next page */
  is_not_at_pageborder = (fpos % pagesize) != 0;
  if (is_not_at_pageborder)
    had_err = fseek(rcr_enc->output,
                    (fpos / pagesize + 1) * pagesize,
                    SEEK_SET);
  if (had_err)
    gt_error_set(err, "error calling fseek: %s", strerror(errno));
  else {
    if (rcr_enc->store_unmmaped_reads) {
      unmapped_reads_filename = gt_str_new_cstr(name);
      gt_str_append_cstr(unmapped_reads_filename, "_unmapped");
      rcr_enc->unmapped_reads_ptr =
        gt_fa_fopen_with_suffix(gt_str_get(unmapped_reads_filename),
                                ".fastq",
                                "w",
                                err);
      if (rcr_enc->unmapped_reads_ptr == NULL)
        had_err = -1;
      gt_str_delete(unmapped_reads_filename);
    }
    else
      rcr_enc->unmapped_reads_ptr = NULL;

    if (!had_err)
      had_err = rcr_write_encoding_to_file(rcr_enc, err);
    gt_fa_xfclose(rcr_enc->output);
    gt_fa_xfclose(rcr_enc->unmapped_reads_ptr);
  }
  return had_err;
}

int gt_rcr_encoder_encode(GtRcrEncoder *rcr_enc, const char *name,
                          GtTimer *timer, GtError *err)
{
  int had_err = 0;
  gt_assert(rcr_enc && name);
  gt_error_check(err);
  if (timer != NULL)
    gt_timer_show_progress(timer, "write encoding", stdout);
  rcr_enc->prev_readpos = 0;
  rcr_enc->cur_read = 0;
  had_err = rcr_write_data(name, rcr_enc, err);
  if (!had_err &&
      rcr_enc->encdesc_enc != NULL) {
    had_err = gt_encdesc_encoder_encode(rcr_enc->encdesc_enc,
                                        rcr_enc->cstr_iterator,
                                        name,
                                        err);
  }

  return had_err;
}

void gt_rcr_encoder_enable_verbosity(GtRcrEncoder *rcr_enc)
{
  gt_assert(rcr_enc);
  rcr_enc->is_verbose = true;
}

void gt_rcr_encoder_disable_verbosity(GtRcrEncoder *rcr_enc)
{
  gt_assert(rcr_enc);
  rcr_enc->is_verbose = false;
}

static void rcr_read_header(GtRcrDecoder *rcr_dec)
{
  unsigned alpha_size;
  GtUword numofleaves,
                m,
                i,
                symbol,
                max_read_length;
  GtUint64 freq;
  GT_UNUSED size_t read,
            one = (size_t) 1;

  GtDiscDistri *readlength_distr,
               *qual_distr,
               *qual_mapping_distr = NULL;

  read = gt_xfread_one(&rcr_dec->numofreads, rcr_dec->fp);
  gt_assert(read == one);
  read = gt_xfread_one(&rcr_dec->cons_readlength, rcr_dec->fp);
  gt_assert(read == one);

  if (rcr_dec->cons_readlength) {
    read = gt_xfread_one(&rcr_dec->readlength, rcr_dec->fp);
    gt_assert(read == one);
  }
  else {
    readlength_distr = gt_disc_distri_new();
    read = gt_xfread_one(&numofleaves, rcr_dec->fp);
    gt_assert(read == one);

    read = gt_xfread_one(&max_read_length, rcr_dec->fp);
    gt_assert(read == one);

    for (i = 0; i < numofleaves; i++) {
      read = gt_xfread_one(&symbol, rcr_dec->fp);
      gt_assert(read == one);
      read = gt_xfread_one(&freq, rcr_dec->fp);
      gt_assert(read == one);
      gt_disc_distri_add_multi(readlength_distr, symbol, freq);
    }
    rcr_dec->readlenghts_huff = gt_huffman_new(readlength_distr,
                                            rcr_disc_distri_func,
                                            max_read_length + 1);
    gt_disc_distri_delete(readlength_distr);
  }

  read = gt_xfread_one(&rcr_dec->store_all_qual, rcr_dec->fp);
  gt_assert(read == one);

  read = gt_xfread_one(&rcr_dec->store_var_qual, rcr_dec->fp);
  gt_assert(read == one);

  if (rcr_dec->store_all_qual || rcr_dec->store_var_qual) {
    qual_distr = gt_disc_distri_new();

    read = gt_xfread_one(&numofleaves, rcr_dec->fp);
    gt_assert(read == one);

    for (i = 0; i < numofleaves; i++) {
      read = gt_xfread_one(&symbol, rcr_dec->fp);
      gt_assert(read == one);
      read = gt_xfread_one(&freq, rcr_dec->fp);
      gt_assert(read == one);
      gt_disc_distri_add_multi(qual_distr, symbol, freq);
    }
    rcr_dec->qual_huff = gt_huffman_new(qual_distr, rcr_disc_distri_func,
                                        256UL);
    gt_disc_distri_delete(qual_distr);
  }

  read = gt_xfread_one(&rcr_dec->store_mapping_qual, rcr_dec->fp);
  gt_assert(read == one);

  if (rcr_dec->store_mapping_qual) {
    qual_mapping_distr= gt_disc_distri_new();
    read = gt_xfread_one(&numofleaves, rcr_dec->fp);
    gt_assert(read == one);
    for (i = 0; i < numofleaves; i++) {
      read = gt_xfread_one(&symbol, rcr_dec->fp);
      gt_assert(read == one);
      read = gt_xfread_one(&freq, rcr_dec->fp);
      gt_assert(read == one);
      gt_disc_distri_add_multi(qual_mapping_distr, symbol, freq);
    }
    rcr_dec->qual_mapping_huff = gt_huffman_new(qual_mapping_distr,
                                             rcr_disc_distri_func, 256UL);
    gt_disc_distri_delete(qual_mapping_distr);
  }
  read = gt_xfread_one(&m, rcr_dec->fp);
  gt_assert(read == one);

  rcr_dec->readpos_golomb = gt_golomb_new(m);

  read = gt_xfread_one(&m, rcr_dec->fp);
  gt_assert(read == one);
  if (m != GT_UNDEF_UWORD) {
    rcr_dec->varpos_golomb = gt_golomb_new(m);
  }

  read = gt_xfread(rcr_dec->present_cigar_ops,
                   sizeof (*rcr_dec->present_cigar_ops),
                   (size_t) (ENDOFRECORD + 1),
                   rcr_dec->fp);
  gt_assert(read == (size_t) (ENDOFRECORD + 1));

  rcr_dec->cigar_ops_huff =
    gt_huffman_new(rcr_dec->present_cigar_ops,
                   rcr_array_func,
                   (GtUword) (ENDOFRECORD + 1));

  alpha_size = gt_alphabet_size(gt_encseq_alphabet(rcr_dec->encseq));
  read = gt_xfread(rcr_dec->ins_bases,
            sizeof (*rcr_dec->ins_bases),
            (size_t) (alpha_size + 1),
            rcr_dec->fp);
  gt_assert(read == (size_t) (alpha_size + 1));

  rcr_dec->bases_huff =
    gt_huffman_new(rcr_dec->ins_bases,
                   rcr_array_func,
                   (GtUword) (alpha_size + 1));
  gt_assert(rcr_dec->bases_huff != NULL);

}

#define RCR_NEXT_BIT(bit)                                                      \
  ({((gt_bitinstream_get_next_bit(bitstream, &bit) != 1) ?                     \
    gt_error_set(err, "could not get next bit on line %d in file %s",          \
                 __LINE__, __FILE__),                                          \
    had_err = -1,                                                              \
    false :                                                                    \
    true);                                                                     \
   })

typedef struct RcrDecodeInfo {
  GtAlphabet                 *alphabet;
  GtEliasGammaBitwiseDecoder *del_len_ebwd;
  GtGolombBitwiseDecoder     *varpos_gbwd;
  GtHuffmanBitwiseDecoder    *base_hbwd,
                             *qual_hbwd,
                             *cigar_hbwd;
  GtStr                      *base_string,
                             *qual_string,
                             *cigar_string;
  GtUword               alpha_size,
                              offset,
                              inserted_bases;
} RcrDecodeInfo;

static int rcr_huff_read(GtHuffmanBitwiseDecoder *hbwd,
                         GtBitInStream *bitstream,
                         GtUword *val,
                         GtError *err)
{
  bool bit;
  int had_err = 0,
      stat = -1;
  while (!had_err && stat != 0) {
    if (RCR_NEXT_BIT(bit)) {
      stat = gt_huffman_bitwise_decoder_next(hbwd, bit, val, err);
      if (stat == -1)
        had_err = stat;
    }
  }
  return had_err;
}

static int rcr_huff_read_string(GtHuffmanBitwiseDecoder *hbwd,
                                GtBitInStream *bitstream,
                                GtUword length,
                                GtStr *str,
                                GtError *err)
{
  int had_err = 0;
  GtUword symbol,
                read = 0;
  while (!had_err && read < length) {
    had_err = rcr_huff_read(hbwd, bitstream, &symbol, err);
    if (!had_err) {
      gt_str_append_char(str, (char) symbol);
      read++;
    }
  }
  return had_err;
}

static int rcr_golomb_read(GtGolombBitwiseDecoder *gbwd,
                           GtBitInStream *bitstream,
                           GtUword *val,
                           GtError *err)
{
  bool bit;
  int had_err = 0,
      stat = -1;
  while (!had_err && stat != 0) {
    if (RCR_NEXT_BIT(bit)) {
      stat = gt_golomb_bitwise_decoder_next(gbwd, bit, val);
      if (stat == -1)
        had_err = stat;
    }
  }
  return had_err;
}

/* XXX struct with all strings, and struct with info */
static void rcr_decode_exact_range(GtRcrDecoder *rcr_dec,
                                   RcrDecodeInfo *info,
                                   GtUword start,
                                   GtUword end)
{
  GtUchar ref;
  GtUword i;
  for (i = start; i < end ; i++) {
    ref = (GtUchar) gt_encseq_get_decoded_char(rcr_dec->encseq, i,
                                               GT_READMODE_FORWARD);
    gt_str_append_char(info->base_string, toupper((int) ref));
    gt_str_append_char(info->cigar_string, '=');
    if (!rcr_dec->store_all_qual)
      gt_str_append_char(info->qual_string, DEFAULTQUAL);
  }
}

static void rcr_decode_exact(GtRcrDecoder *rcr_dec,
                             RcrDecodeInfo *info,
                             GtUword readlength,
                             GtUword seqstart,
                             GtUword readpos)
{
  rcr_decode_exact_range(rcr_dec, info, seqstart + readpos,
                         seqstart + readpos + readlength);
}

static int rcr_decode_mismatch(GtRcrDecoder *rcr_dec, GtUword pos,
                               RcrDecodeInfo *info, GtBitInStream *bitstream,
                               GtError *err)
{
  bool bit;
  int had_err = 0;
  GtBitsequence coded = 0,
                first_bit = (GtBitsequence) 1 << 1,
                second_bit = (GtBitsequence) 1;
  GtUchar base,
          ref = gt_encseq_get_encoded_char(rcr_dec->encseq, pos,
                                           GT_READMODE_FORWARD);

  if (RCR_NEXT_BIT(bit)) {
    if (bit)
      coded = first_bit;
  }
  if (!had_err && RCR_NEXT_BIT(bit)) {
    if (bit)
      coded |= second_bit;

    base = rcr_transdecode(ref, coded, info->alphabet);

    gt_str_append_char(info->base_string,
                       toupper(gt_alphabet_decode(info->alphabet, base)));
  }
  return had_err;
}

static int rcr_decode_mismatch_qual(GtRcrDecoder *rcr_dec,
                                    GtBitInStream *bitstream,
                                    RcrDecodeInfo *info,
                                    GtError *err)
{
  int had_err = 0;
  GtUword symbol = 0;

  if (rcr_dec->store_var_qual) {
    had_err = rcr_huff_read(info->qual_hbwd, bitstream, &symbol, err);
    gt_str_append_char(info->qual_string, (char) symbol);
  }
  else {
    if (!rcr_dec->store_all_qual)
      gt_str_append_char(info->qual_string, DEFAULTQUAL);
  }
  return had_err;
}

static int rcr_decode_insert_var(GtRcrDecoder *rcr_dec,
                                 GtBitInStream *bitstream,
                                 RcrDecodeInfo *info,
                                 char type,
                                 GtError *err)
{
  int had_err = 0;
  GtUword symbol = 0, i;
  GtUchar base;

  info->inserted_bases = 0;

  had_err = rcr_huff_read(info->base_hbwd, bitstream, &symbol, err);
  while (!had_err && symbol != info->alpha_size) {
    if (symbol == (info->alpha_size - 1)) {
      base = (GtUchar) WILDCARD;
    }
    else
      base = (GtUchar) symbol;

    gt_str_append_char(info->base_string,
                       toupper(gt_alphabet_decode(info->alphabet, base)));
    info->inserted_bases++;
    had_err = rcr_huff_read(info->base_hbwd, bitstream, &symbol, err);
  }

  if (!had_err) {
    for (i = 0; i < info->inserted_bases; i++)
      gt_str_append_char(info->cigar_string, type);

    if (rcr_dec->store_var_qual) {
      had_err = rcr_huff_read_string(info->qual_hbwd, bitstream,
                                     info->inserted_bases, info->qual_string,
                                     err);
    }
    else {
      if (!rcr_dec->store_all_qual) {
        for (i = 0; i < info->inserted_bases; i++)
          gt_str_append_char(info->qual_string, DEFAULTQUAL);
      }
    }
  }
  return had_err;
}

static inline int rcr_elias_read(GtEliasGammaBitwiseDecoder *ebwd,
                                 GtBitInStream *bitstream,
                                 GtUword *symbol,
                                 GtError *err)
{
  bool bit;
  int had_err = 0,
      stat = -1;
  while (!had_err && stat != 0) {
    if (RCR_NEXT_BIT(bit))
      stat = gt_elias_gamma_bitwise_decoder_next(ebwd, bit, symbol);
  }
  return had_err;
}

static int rcr_decode_delete_var(RcrDecodeInfo *info,
                                 GtBitInStream *bitstream,
                                 char type,
                                 GtError *err)
{
  /* XXX remove loop and just set in as xD */
  int had_err = 0;
  GtUword del_length , i;
  had_err = rcr_elias_read(info->del_len_ebwd, bitstream, &del_length, err);
  if (!had_err) {
    for (i = 0; i < del_length; i++)
      gt_str_append_char(info->cigar_string, type);
    info->offset = del_length;
  }
  return had_err;
}

static RcrDecodeInfo *rcr_init_decode_info(GtRcrDecoder *rcr_dec, GtError *err)
{
  RcrDecodeInfo *info = gt_malloc(sizeof (*info));

  gt_error_check(err);

  info->base_hbwd = NULL;
  info->cigar_hbwd = NULL;
  info->qual_hbwd = NULL;
  info->del_len_ebwd = NULL;

  info->base_hbwd = gt_huffman_bitwise_decoder_new(rcr_dec->bases_huff, err);
  if (info->base_hbwd != NULL)
    info->cigar_hbwd =
      gt_huffman_bitwise_decoder_new(rcr_dec->cigar_ops_huff, err);
  if (info->cigar_hbwd != NULL) {
    info->del_len_ebwd = gt_elias_gamma_bitwise_decoder_new();

    info->alphabet = gt_encseq_alphabet(rcr_dec->encseq);
    info->alpha_size = (GtUword) gt_alphabet_size(info->alphabet);
    info->qual_hbwd = NULL;
    info->varpos_gbwd = NULL;
    info->base_string = gt_str_new();
    info->qual_string = gt_str_new();
    info->cigar_string = gt_str_new();
    info->inserted_bases = 0;

    if (rcr_dec->varpos_golomb != NULL)
      info->varpos_gbwd = gt_golomb_bitwise_decoder_new(rcr_dec->varpos_golomb);

    if (rcr_dec->store_var_qual || rcr_dec->store_all_qual) {
      info->qual_hbwd = gt_huffman_bitwise_decoder_new(rcr_dec->qual_huff, err);
      if (info->qual_hbwd == NULL)
        info = NULL;
    }
  }
  else {
    info = NULL;
  }
  return info;
}

static void rcr_delete_decode_info(RcrDecodeInfo *info)
{
  if (info != NULL) {
    gt_elias_gamma_bitwise_decoder_delete(info->del_len_ebwd);
    gt_golomb_bitwise_decoder_delete(info->varpos_gbwd);
    gt_huffman_bitwise_decoder_delete(info->base_hbwd);
    gt_huffman_bitwise_decoder_delete(info->cigar_hbwd);
    gt_huffman_bitwise_decoder_delete(info->qual_hbwd);
    gt_str_delete(info->base_string);
    gt_str_delete(info->cigar_string);
    gt_str_delete(info->qual_string);
    gt_free(info);
  }
}

static inline int rcr_decode_inexact(GtRcrDecoder *rcr_dec,
                                     GtBitInStream *bitstream,
                                     RcrDecodeInfo *info,
                                     GtUword seq_i,
                                     GtUword readlength,
                                     GtError *err)
{
  int had_err = 0,
      cigar_op = GT_UNDEF_INT;
  GtUword end,
                prev_varpos = 0,
                read_i = 0,
                rel_varpos = 0,
                symbol,
                varpos;

  info->offset = 0;

  /* read variation type */
  had_err = rcr_huff_read(info->cigar_hbwd, bitstream, &symbol, err);
  if (!had_err)
    gt_safe_assign(cigar_op, symbol);

  while (!had_err && cigar_op != ENDOFRECORD) {

    /* read variation position */
    if (!had_err)
      had_err = rcr_golomb_read(info->varpos_gbwd, bitstream, &rel_varpos, err);

    varpos = rel_varpos + prev_varpos;

    /* are there matches between the current and the previous variation? */
    if (!had_err && (read_i < varpos)) {
      end = seq_i + varpos - read_i;
      rcr_decode_exact_range(rcr_dec, info, seq_i, end);
      seq_i = end;
      read_i += varpos - read_i;
    }

    /* variation handling */
    if (!had_err) {
      switch (cigar_op) {
        case  BAM_CMATCH:
          gt_str_append_char(info->cigar_string, 'X');
          had_err = rcr_decode_mismatch(rcr_dec, seq_i, info, bitstream, err);
          if (!had_err)
            had_err = rcr_decode_mismatch_qual(rcr_dec, bitstream, info, err);
          read_i++; seq_i++;
          break;
        case BAM_CSOFT_CLIP:
          had_err = rcr_decode_insert_var(rcr_dec, bitstream, info, 'S', err);
          read_i += info->inserted_bases;
          break;
        case BAM_CINS:
          had_err = rcr_decode_insert_var(rcr_dec, bitstream, info, 'I', err);
          read_i += info->inserted_bases;
          break;
        case BAM_CDEL:
          had_err = rcr_decode_delete_var(info, bitstream, 'D', err);
          seq_i += info->offset;
          break;
        case BAM_CREF_SKIP:
          had_err = rcr_decode_delete_var(info, bitstream, 'N', err);
          seq_i += info->offset;
          break;
        default:
          gt_error_set(err, "encountered funny cigar op: %d", cigar_op);
          had_err = -1;
      }
    }
    prev_varpos = varpos;
    /* read next variation type */
    if (!had_err)
      had_err = rcr_huff_read(info->cigar_hbwd, bitstream, &symbol, err);
    if (!had_err)
      gt_safe_assign(cigar_op, symbol);
  }
  /* end symbol was read are there trailing matches?*/
  if (!had_err && (read_i < readlength)) {
      end = seq_i + readlength - read_i;
      rcr_decode_exact_range(rcr_dec, info, seq_i, end);
  }
  return had_err;
}

static int rcr_write_decoding_to_file(GtRcrDecoder *rcr_dec, GtError *err)
{
  bool bit,
       strand;
  int had_err = 0;
  uint32_t mapping_qual = 0;
  GtUword
                cur_read = 0,
                i,
                l,
                prev_readpos = 0,
                readlength,
                readpos = 0,
                refnum = 0,
                rel_readpos,
                seqstart = 0,
                symbol;
  GtStr *qname;
  GtHuffmanBitwiseDecoder *readlen_hbwd = NULL,
                          *mapping_qual_hbwd = NULL;
  GtGolombBitwiseDecoder *readpos_gbwd = NULL;
  GtBitInStream *bitstream;
  RcrDecodeInfo *info = rcr_init_decode_info(rcr_dec, err);

  if (info == NULL)
    /* had_err = -1 is nicer, but splint does not like it */
    return -1;

  if (!had_err && !rcr_dec->cons_readlength) {
    readlen_hbwd =
      gt_huffman_bitwise_decoder_new(rcr_dec->readlenghts_huff, err);
  }

  if (!had_err && rcr_dec->store_mapping_qual) {
    mapping_qual_hbwd =
      gt_huffman_bitwise_decoder_new(rcr_dec->qual_mapping_huff, err);
  }

  if (!had_err) {
    readpos_gbwd =
      gt_golomb_bitwise_decoder_new(rcr_dec->readpos_golomb);

    qname = gt_str_new();

    bitstream = gt_bitinstream_new(gt_str_get(rcr_dec->inputname),
                                   (size_t) rcr_dec->startofencoding,
                                   1UL);

    for (i = 0; i < gt_encseq_num_of_sequences(rcr_dec->encseq); i++) {
      const char *seqname = gt_encseq_description(rcr_dec->encseq, &l, i);
      GtUword len = gt_encseq_seqlength(rcr_dec->encseq, i);
      fprintf(rcr_dec->fp, "@SQ\tSN:%.*s\tLN:"GT_WU"\n", (int) l, seqname, len);
    }

    gt_log_log("start to decode "GT_WU" reads", rcr_dec->numofreads);
  }

  while (!had_err &&
         cur_read < rcr_dec->numofreads) {

    /* check if there is a new seq in encseq */
    if (RCR_NEXT_BIT(bit)) {
      if (bit) {
        seqstart = gt_encseq_seqstartpos(rcr_dec->encseq, refnum);
        gt_log_log("get start for new ref "GT_WU": "GT_WU"", refnum, seqstart);
        refnum++;
        prev_readpos = 0;
      }
    }
    /* check if read was unmapped */
    if (!had_err && RCR_NEXT_BIT(bit)) {
      if (bit) {
        /* XXX this is ugly */
        continue;
      }
    }

    gt_str_reset(qname);
    /* read read name */
    if (!had_err && rcr_dec->encdesc != NULL) {
      if (gt_encdesc_decode(rcr_dec->encdesc, cur_read, qname, err) != 1) {
        had_err = -1;
      }
    }
    else
      gt_str_append_ulong(qname, cur_read);

    /* read read length */
    if (!had_err) {
      if (rcr_dec->cons_readlength)
        readlength = rcr_dec->readlength;
      else
        had_err = rcr_huff_read(readlen_hbwd, bitstream, &readlength, err);
    }

    /* read read position */
    if (!had_err) {
      had_err = rcr_golomb_read(readpos_gbwd, bitstream, &rel_readpos, err);
      if (!had_err) {
        readpos = rel_readpos + prev_readpos;
        prev_readpos = readpos;
      }
    }

    /* read mapping qual */
    if (!had_err && rcr_dec->store_mapping_qual) {
      had_err = rcr_huff_read(mapping_qual_hbwd, bitstream, &symbol, err);
      if (!had_err) {
        gt_safe_assign(mapping_qual, symbol);
      }
    }

    /* read qual string */
    if (!had_err && rcr_dec->store_all_qual)
      had_err = rcr_huff_read_string(info->qual_hbwd, bitstream, readlength,
                                     info->qual_string, err);

    if (!had_err) {
      /* read strand */
      if (RCR_NEXT_BIT(bit)) {
        strand = bit;
      }
    }
    if (!had_err) {
      /* exact match? */
      GtUword seq_i = seqstart + readpos;
      if (RCR_NEXT_BIT(bit)) {
        if (bit)
          rcr_decode_exact(rcr_dec, info, readlength, seqstart, readpos);
        else
          had_err = rcr_decode_inexact(rcr_dec, bitstream, info, seq_i,
                                       readlength, err);

        /* write read to file */
        if (readlength != gt_str_length(info->base_string)) {
          gt_log_log("readlen: "GT_WU", stringlen: "GT_WU", read: "GT_WU"",
                     readlength, gt_str_length(info->base_string), cur_read);
        }
        gt_assert(readlength == gt_str_length(info->base_string));
        gt_assert(readlength == gt_str_length(info->qual_string));
        fprintf(rcr_dec->fp, "%s", gt_str_get(qname));
        fprintf(rcr_dec->fp, "\t%c", strand?'-':'+');
        fprintf(rcr_dec->fp, "\t"GT_WU"", readpos + 1);
        if (rcr_dec->store_mapping_qual)
          fprintf(rcr_dec->fp, "\t%u", (unsigned) mapping_qual);
        else
          fprintf(rcr_dec->fp, "\t%u", DEFAULTMQUAL);

        rcr_convert_cigar_string(info->cigar_string);
        fprintf(rcr_dec->fp, "\t%s", gt_str_get(info->cigar_string));
        fprintf(rcr_dec->fp,"\t%s", gt_str_get(info->base_string));
        fprintf(rcr_dec->fp, "\t%s\n", gt_str_get(info->qual_string));
        gt_str_reset(info->cigar_string);
        gt_str_reset(info->qual_string);
        gt_str_reset(info->base_string);
        cur_read++;
      }
    }
  }
  gt_str_delete(qname);
  gt_huffman_bitwise_decoder_delete(readlen_hbwd);
  gt_huffman_bitwise_decoder_delete(mapping_qual_hbwd);
  gt_golomb_bitwise_decoder_delete(readpos_gbwd);
  gt_bitinstream_delete(bitstream);
  rcr_delete_decode_info(info);
  gt_log_log("decoded "GT_WU" reads", cur_read);
  return had_err;
}

static GtRcrDecoder *gt_rcr_decoder_init(const char *name,
                                         const GtEncseq *ref,
                                         GtError *err)
{
  size_t alpha_size =
    (size_t) gt_alphabet_size(gt_encseq_alphabet(ref));
  GtRcrDecoder *rcr_dec = gt_malloc(sizeof (GtRcrDecoder));

  rcr_dec->basename = name;
  rcr_dec->inputname = gt_str_new_cstr(name);
  gt_str_append_cstr(rcr_dec->inputname,RCRFILESUFFIX);
  rcr_dec->encseq = ref;

  rcr_dec->encdesc = NULL;
  rcr_dec->qual_huff = NULL;
  rcr_dec->qual_mapping_huff = NULL;
  rcr_dec->readlenghts_huff = NULL;
  rcr_dec->readpos_golomb = NULL;
  rcr_dec->varpos_golomb = NULL;

  rcr_dec->ins_bases = gt_calloc(alpha_size + 1, sizeof (GtUint64));

  rcr_dec->fp = gt_fa_fopen(gt_str_get(rcr_dec->inputname),
                            "rb", err);
  if (rcr_dec->fp == NULL) {
    gt_rcr_decoder_delete(rcr_dec);
    return NULL;
  }
  return rcr_dec;
}

GtRcrDecoder *gt_rcr_decoder_new(const char *name, const GtEncseq *ref,
                                 GtTimer *timer, GtError *err)
{
  bool is_not_at_pageborder;
  GtUword pagesize = gt_pagesize();
  GtWord filepos;
  GtRcrDecoder *rcr_dec;

  gt_assert(name);

  if (timer != NULL)
    gt_timer_show_progress(timer, "initializing RcrDecoder", stdout);
  gt_error_check(err);
  if (!(gt_alphabet_is_dna(gt_encseq_alphabet(ref)))) {
    gt_error_set(err, "alphabet has to be DNA");
    return NULL;
  }
  rcr_dec = gt_rcr_decoder_init(name, ref, err);

  rcr_read_header(rcr_dec);

  gt_assert(rcr_dec != NULL && rcr_dec->fp != NULL);
  filepos = ftell(rcr_dec->fp);
  is_not_at_pageborder = (filepos % pagesize) != 0;
  if (is_not_at_pageborder)
    filepos = (filepos / pagesize + 1) * pagesize;

  gt_safe_assign(rcr_dec->startofencoding, filepos);
  gt_fa_fclose(rcr_dec->fp);
  return rcr_dec;
}

int gt_rcr_decoder_enable_description_support(GtRcrDecoder *rcr_dec,
                                              GtError *err)
{
  gt_assert(rcr_dec != NULL);
  rcr_dec->encdesc = gt_encdesc_load(rcr_dec->basename, err);
  if (rcr_dec->encdesc == NULL)
    return -1;
  return 0;
}

void gt_rcr_decoder_disable_description_support(GtRcrDecoder *rcr_dec)
{
  gt_assert(rcr_dec != NULL);
  gt_encdesc_delete(rcr_dec->encdesc);
  rcr_dec->encdesc = NULL;
}

int gt_rcr_decoder_decode(GtRcrDecoder *rcr_dec,
                          const char *name,
                          GtTimer *timer,
                          GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(rcr_dec && name);

  if (timer != NULL)
    gt_timer_show_progress(timer, "decode rcr", stdout);
  if (!had_err) {
    rcr_dec->fp = gt_fa_fopen_with_suffix(name, ".rcr.decoded", "w", err);
    if (rcr_dec->fp == NULL)
      had_err = -1;
  }
  if (!had_err)
    had_err = rcr_write_decoding_to_file(rcr_dec, err);

  gt_fa_xfclose(rcr_dec->fp);
  return had_err;
}

void gt_rcr_encoder_delete(GtRcrEncoder *rcr_enc)
{
  if (rcr_enc != NULL) {
    gt_huffman_delete(rcr_enc->readlenghts_huff);
    gt_huffman_delete(rcr_enc->qual_mapping_huff);
    gt_huffman_delete(rcr_enc->qual_huff);
    gt_huffman_delete(rcr_enc->cigar_ops_huff);
    gt_huffman_delete(rcr_enc->bases_huff);

    gt_golomb_delete(rcr_enc->readpos_golomb);
    gt_golomb_delete(rcr_enc->varpos_golomb);

    gt_queue_delete(rcr_enc->not_exact_matches);

    bam_destroy1(rcr_enc->sam_align);

    gt_disc_distri_delete(rcr_enc->qual_mapping_distr);
    gt_disc_distri_delete(rcr_enc->readlength_distr);
    gt_disc_distri_delete(rcr_enc->readpos_distr);
    gt_disc_distri_delete(rcr_enc->varpos_distr);
    gt_disc_distri_delete(rcr_enc->qual_distr);

    gt_samfile_iterator_delete(rcr_enc->sam_iter);
    gt_cstr_iterator_delete(rcr_enc->cstr_iterator);
    gt_encdesc_encoder_delete(rcr_enc->encdesc_enc);

    gt_free(rcr_enc->ins_bases);
    gt_free(rcr_enc);
  }
}

void gt_rcr_decoder_delete(GtRcrDecoder *rcr_dec)
{
  if (rcr_dec != NULL) {
    gt_huffman_delete(rcr_dec->readlenghts_huff);
    gt_huffman_delete(rcr_dec->qual_mapping_huff);
    gt_huffman_delete(rcr_dec->cigar_ops_huff);
    gt_huffman_delete(rcr_dec->qual_huff);
    gt_huffman_delete(rcr_dec->bases_huff);

    gt_golomb_delete(rcr_dec->readpos_golomb);
    gt_golomb_delete(rcr_dec->varpos_golomb);

    gt_str_delete(rcr_dec->inputname);

    gt_encdesc_delete(rcr_dec->encdesc);

    gt_free(rcr_dec->ins_bases);

    gt_free(rcr_dec);
  }
}
