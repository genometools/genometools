/*
  Copyright (c) 2011-2012 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#ifndef S_SPLINT_S
#include <ctype.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "core/array2dim_api.h"
#include "core/array_api.h"
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/compat.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/safearith.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/bitoutstream.h"
#include "extended/encdesc.h"
#include "extended/fasta_header_iterator.h"
#include "extended/hcr.h"
#include "extended/rbtree.h"
#include "extended/cstr_iterator.h"

#define HCR_LOWESTQUALVALUE 0
#define HCR_HIGHESTQUALVALUE 127U
#define HCR_LINEWIDTH 80UL
#define HCR_DESCSEPSEQ '@'
#define HCR_DESCSEPQUAL '+'
#define HCR_PAGES_PER_CHUNK 10UL

typedef struct GtBaseQualDistr {
  GtUint64 **distr;
  GtAlphabet          *alpha;
  unsigned int         ncols,
                       nrows,
                       max_qual,
                       min_qual,
                       wildcard_indx,
                       qrange_start,
                       qual_offset,
                       qrange_end;
} GtBaseQualDistr;

typedef struct FastqFileInfo {
  GtUword readlength,
                readnum;
} FastqFileInfo;

typedef struct GtHcrSeqEncoder {
  GtQualRange        qrange;
  FastqFileInfo     *fileinfos;
  GtAlphabet        *alpha;
  GtHuffman         *huffman;
  GtSampling        *sampling;
  GtUint64 total_num_of_symbols;
  GtUword      num_of_files;
  GtWord               startofsamplingtab,
                     start_of_encoding;
  unsigned int       qual_offset;
} GtHcrSeqEncoder;

struct GtHcrEncoder {
  GtEncdescEncoder *encdesc_encoder;
  GtHcrSeqEncoder  *seq_encoder;
  GtStrArray       *files;
  GtUword     num_of_reads,
                    num_of_files,
                    sampling_rate,
                    pagesize;
  bool              page_sampling,
                    regular_sampling;
};

typedef struct hcr_huff_mem_info {
  char         *path;
  void         *data;
  size_t        start,
                end,
                pos;
  GtUword pages_per_chunk,
                bitseq_per_chunk,
                pagesize,
                blocksize;
  bool          first;
} HcrHuffDataIterator;

typedef struct GtHcrSeqDecoder {
  FastqFileInfo       *fileinfos;
  GtAlphabet          *alpha;
  GtArray             *symbols;
  GtBitsequence       *bitseq_buffer;
  GtHuffman           *huffman;
  GtHuffmanDecoder    *huff_dec;
  GtRBTree            *file_info_rbt;
  GtSampling          *sampling;
  GtStr               *filename;
  HcrHuffDataIterator *data_iter;
  GtUword        readlength,
                       cur_read,
                       num_of_reads,
                       num_of_files;
  GtWord                 start_of_encoding;
  unsigned int         alphabet_size,
                       qual_offset;
} GtHcrSeqDecoder;

struct GtHcrDecoder {
  GtEncdesc       *encdesc;
  GtHcrSeqDecoder *seq_dec;
};

typedef struct WriteNodeInfo {
  FILE       *output;
  GtAlphabet *alpha;
  unsigned    qual_offset;
} WriteNodeInfo;

static GtUint64 hcr_base_qual_distr_func(const void *distr,
                                                   GtUword symbol)
{
  GtUint64 value;
  GtBaseQualDistr *bqd = (GtBaseQualDistr*)distr;
  value = bqd->distr[symbol / bqd->ncols][symbol % bqd->ncols];
  return value;
}

static void hcr_base_qual_distr_trim(GtBaseQualDistr *bqd)
{
  if (bqd->min_qual != 0) {
    GtUint64 **distr_trimmed;
    unsigned nrows_new,
             i,
             j;

    nrows_new = bqd->max_qual - bqd->min_qual + 1;
    gt_array2dim_calloc(distr_trimmed, (GtUword) nrows_new, bqd->ncols);

    for (i = 0; i < nrows_new; i++)
      for (j = 0; j < bqd->ncols; j++)
        distr_trimmed[i][j] = bqd->distr[i + bqd->min_qual][j];
    gt_array2dim_delete(bqd->distr);
    bqd->distr = distr_trimmed;
    bqd->nrows = nrows_new;
    bqd->qual_offset = bqd->min_qual;
  }
}

static GtBaseQualDistr* hcr_base_qual_distr_new_from_file(FILE *fp,
                                                          GtAlphabet *alpha)
{
  GtBaseQualDistr *bqd;
  char read_char_code;
  GtUchar cur_char_code;
  unsigned char cur_qual;
  unsigned alpha_size,
           min_qual = HCR_HIGHESTQUALVALUE,
           max_qual = HCR_LOWESTQUALVALUE;
  GtUword numofleaves,
                i;
  GtUint64 cur_freq;
  GT_UNUSED size_t read,
            one = (size_t) 1;

  alpha_size = gt_alphabet_size(alpha);
  bqd = gt_malloc(sizeof (GtBaseQualDistr));
  gt_array2dim_calloc(bqd->distr, HCR_HIGHESTQUALVALUE + 1UL, alpha_size)
  bqd->ncols = alpha_size;
  bqd->nrows = HCR_HIGHESTQUALVALUE + 1U;
  bqd->qual_offset = HCR_LOWESTQUALVALUE;
  bqd->wildcard_indx = alpha_size - 1;

  read = gt_xfread_one(&numofleaves, fp);
  gt_assert(read == one);
  for (i = 0; i < numofleaves; i++) {
    read = gt_xfread_one(&read_char_code, fp);
    gt_assert(read == one);
    read = gt_xfread_one(&cur_qual, fp);
    gt_assert(read == one);
    read = gt_xfread_one(&cur_freq, fp);
    gt_assert(read == one);
    cur_char_code = gt_alphabet_encode(alpha, read_char_code);
    if (cur_char_code == (GtUchar) WILDCARD)
      gt_safe_assign(cur_char_code, bqd->wildcard_indx);
    bqd->distr[cur_qual][cur_char_code] = cur_freq;
    if ((unsigned) cur_qual > max_qual)
      max_qual = cur_qual;
    if ((unsigned) cur_qual < min_qual)
      min_qual = cur_qual;
  }

  bqd->min_qual = min_qual;
  bqd->max_qual = max_qual;
  hcr_base_qual_distr_trim(bqd);
  return bqd;
}

static GtBaseQualDistr* hcr_base_qual_distr_new(GtAlphabet *alpha,
                                                GtQualRange qrange)
{
  GtBaseQualDistr *bqd;
  bqd = gt_calloc((size_t) 1, sizeof (GtBaseQualDistr));
  gt_array2dim_calloc(bqd->distr, HCR_HIGHESTQUALVALUE + 1UL,
                      gt_alphabet_size(alpha));

  bqd->ncols = gt_alphabet_size(alpha);
  bqd->nrows = HCR_HIGHESTQUALVALUE + 1U;
  bqd->qual_offset = HCR_LOWESTQUALVALUE;
  bqd->wildcard_indx = gt_alphabet_size(alpha) - 1;
  bqd->min_qual = HCR_HIGHESTQUALVALUE;
  bqd->max_qual = HCR_LOWESTQUALVALUE;
  gt_safe_assign(bqd->qrange_start, qrange.start);
  gt_safe_assign(bqd->qrange_end, qrange.end);
  bqd->alpha = alpha;
  return bqd;
}

static int hcr_base_qual_distr_add(GtBaseQualDistr *bqd, const GtUchar *qual,
                                   const GtUchar *seq, GtUword len)
{
  GtUword i;
  unsigned cur_char_code,
           cur_qual;

  for (i = 0; i < len; i++) {
    cur_char_code = (unsigned) gt_alphabet_encode(bqd->alpha,
                           (char) gt_alphabet_pretty_symbol(bqd->alpha,
                                                            (unsigned) seq[i]));
    cur_qual = (unsigned) qual[i];

    if (bqd->qrange_start != GT_UNDEF_UINT) {
      if (cur_qual <= bqd->qrange_start)
        cur_qual = bqd->qrange_start;
    }

    if (bqd->qrange_end != GT_UNDEF_UINT) {
      if (cur_qual >= bqd->qrange_end)
        cur_qual = bqd->qrange_end;
    }
    if (cur_char_code == WILDCARD)
      bqd->distr[cur_qual][bqd->wildcard_indx]++;
    else
      bqd->distr[cur_qual][cur_char_code]++;
    if (cur_qual > bqd->max_qual)
      bqd->max_qual = cur_qual;
    if (cur_qual < bqd->min_qual)
      bqd->min_qual = cur_qual;
  }
  return 0;
}

static void hcr_base_qual_distr_delete(GtBaseQualDistr *bqd)
{
  if (!bqd)
    return;
  gt_array2dim_delete(bqd->distr);
  gt_free(bqd);
}

static int hcr_cmp_FastqFileInfo(const void *node1, const void *node2,
                                 GT_UNUSED void *unused)
{
  FastqFileInfo *n1 = (FastqFileInfo*) node1,
                *n2 = (FastqFileInfo*) node2;

  gt_assert(n1 && n2);

  if (n1->readnum < n2->readnum)
    return -1;
  if (n1->readnum > n2->readnum)
    return 1;
  return 0;
}

static GtUword hcr_write_seq(GtHcrSeqEncoder *seq_encoder,
                                   const GtUchar *seq,
                                   const GtUchar *qual,
                                   GtUword len,
                                   GtBitOutStream *bitstream,
                                   bool dry)
{
  unsigned bits_to_write,
           cur_char_code,
           cur_qual,
           symbol;
  GtUword i,
                written_bits = 0;
  GtBitsequence code;

  for (i = 0; i < len; i++) {
    cur_char_code = (unsigned) seq[i];

    if (cur_char_code == WILDCARD)
      cur_char_code = gt_alphabet_size(seq_encoder->alpha) - 1;

    cur_qual = (unsigned) qual[i];

    if (seq_encoder->qrange.start != GT_UNDEF_UINT) {
      if (cur_qual <= seq_encoder->qrange.start)
        cur_qual = seq_encoder->qrange.start;
    }

    if (seq_encoder->qrange.end != GT_UNDEF_UINT) {
      if (cur_qual >= seq_encoder->qrange.end)
        cur_qual = seq_encoder->qrange.end;
    }

    cur_qual = cur_qual - seq_encoder->qual_offset;

    symbol = gt_alphabet_size(seq_encoder->alpha) * cur_qual + cur_char_code;
    gt_huffman_encode(seq_encoder->huffman, (GtUword) symbol,
                      &code, &bits_to_write);
    written_bits += bits_to_write;
    if (!dry) {
      gt_bitoutstream_append(bitstream, code, bits_to_write);
    }
  }
  return written_bits;
}

static int hcr_write_seqs(FILE *fp, GtHcrEncoder *hcr_enc, GtError *err)
{
  int had_err = 0, seqit_err;
  GtUword bits_to_write = 0,
                len,
                read_counter = 0,
                page_counter = 0,
                bits_left_in_page,
                cur_read = 0;
  GtWord filepos;
  GtSeqIterator *seqit;
  const GtUchar *seq,
                *qual;
  char *desc;
  GtBitOutStream *bitstream;

  gt_error_check(err);
  gt_assert(hcr_enc->seq_encoder->sampling);

  gt_safe_assign(bits_left_in_page, (hcr_enc->pagesize * 8));

  gt_xfseek(fp, hcr_enc->seq_encoder->start_of_encoding, SEEK_SET);
  bitstream = gt_bitoutstream_new(fp);

  seqit = gt_seq_iterator_fastq_new(hcr_enc->files, err);
  if (!seqit) {
    gt_assert(gt_error_is_set(err));
    had_err = -1;
  }

  if (!had_err) {
    gt_seq_iterator_set_quality_buffer(seqit, &qual);
    gt_seq_iterator_set_symbolmap(seqit,
                            gt_alphabet_symbolmap(hcr_enc->seq_encoder->alpha));
    hcr_enc->seq_encoder->total_num_of_symbols = 0;
    while (!had_err &&
           (seqit_err = gt_seq_iterator_next(seqit,
                                            &seq,
                                            &len,
                                            &desc, err)) == 1) {

      /* count the bits */
      bits_to_write = hcr_write_seq(hcr_enc->seq_encoder, seq, qual, len,
                                    bitstream, true);

      /* check if a new sample has to be added */
      if (gt_sampling_is_next_element_sample(hcr_enc->seq_encoder->sampling,
                                             page_counter,
                                             read_counter,
                                             bits_to_write,
                                             bits_left_in_page)) {
        gt_bitoutstream_flush_advance(bitstream);

        filepos = gt_bitoutstream_pos(bitstream);
        if (filepos < 0) {
          had_err = -1;
          gt_error_set(err, "error by ftell: %s", strerror(errno));
        }
        else {
        gt_sampling_add_sample(hcr_enc->seq_encoder->sampling,
                               (size_t) filepos,
                               cur_read);

        read_counter = 0;
        page_counter = 0;
        gt_safe_assign(bits_left_in_page, (hcr_enc->pagesize * 8));
        }
      }

      if (!had_err) {
      /* do the writing */
      bits_to_write = hcr_write_seq(hcr_enc->seq_encoder,
                                    seq, qual, len, bitstream, false);

      /* update counter for sampling */
      while (bits_left_in_page < bits_to_write) {
        page_counter++;
        bits_to_write -= bits_left_in_page;
        gt_safe_assign(bits_left_in_page, (hcr_enc->pagesize * 8));
      }
      bits_left_in_page -= bits_to_write;
      /* always set first page as written */
      if (page_counter == 0)
        page_counter++;
      read_counter++;
      hcr_enc->seq_encoder->total_num_of_symbols += len;
      cur_read++;
      }
    }
    gt_assert(hcr_enc->num_of_reads == cur_read);
    if (!had_err && seqit_err) {
      had_err = seqit_err;
      gt_assert(gt_error_is_set(err));
    }
  }

  if (!had_err) {
    gt_bitoutstream_flush(bitstream);
    filepos = gt_bitoutstream_pos(bitstream);
    if (filepos < 0) {
      had_err = -1;
      gt_error_set(err, "error by ftell: %s", strerror(errno));
    }
    else {
    hcr_enc->seq_encoder->startofsamplingtab = filepos;
    gt_log_log("start of samplingtab: "GT_WU"",
               hcr_enc->seq_encoder->startofsamplingtab);
    if (hcr_enc->seq_encoder->sampling != NULL)
      gt_sampling_write(hcr_enc->seq_encoder->sampling, fp);
    }
  }
  gt_bitoutstream_delete(bitstream);
  gt_seq_iterator_delete(seqit);
  return had_err;
}

static int hcr_huffman_write_base_qual_freq(GtUword symbol,
                                            GtUint64 freq,
                                            GT_UNUSED GtBitsequence code,
                                            GT_UNUSED unsigned code_length,
                                            void *pt)
{
  GtUchar base,
          qual;
  WriteNodeInfo *info = (WriteNodeInfo*)pt;

  gt_safe_assign(base, (symbol % gt_alphabet_size(info->alpha)));
  if (base == (GtUchar) gt_alphabet_size(info->alpha) - 1)
    base = (GtUchar) WILDCARD;
  gt_safe_assign(base, (toupper(gt_alphabet_decode(info->alpha, base))));

  gt_xfwrite_one(&base, info->output);

  gt_safe_assign(qual,
                 (symbol / gt_alphabet_size(info->alpha) + info->qual_offset));

  gt_xfwrite_one(&qual, info->output);
  gt_xfwrite_one(&freq, info->output);
  return 0;
}

static int hcr_write_seqdistrtab(FILE *fp, GtHcrEncoder *hcr_enc)
{
  WriteNodeInfo *info;
  int had_err = 0;
  GtUword numofleaves;

  info = gt_calloc((size_t) 1, sizeof (WriteNodeInfo));
  info->alpha = hcr_enc->seq_encoder->alpha;
  info->qual_offset = hcr_enc->seq_encoder->qual_offset;
  info->output = fp;

  numofleaves = gt_huffman_numofsymbols(hcr_enc->seq_encoder->huffman);
  gt_xfwrite_one(&numofleaves, fp);

  had_err = gt_huffman_iterate(hcr_enc->seq_encoder->huffman,
                               hcr_huffman_write_base_qual_freq,
                               info);
  gt_free(info);
  return had_err;
}

static void hcr_write_file_info(FILE *fp, GtHcrEncoder *hcr_enc)
{
  GtUword i;

  gt_xfwrite_one(&hcr_enc->num_of_files, fp);

  for (i = 0; i < hcr_enc->num_of_files; i++) {
    gt_xfwrite_one(&hcr_enc->seq_encoder->fileinfos[i].readnum, fp);
    gt_xfwrite_one(&hcr_enc->seq_encoder->fileinfos[i].readlength, fp);
  }
}

static int hcr_write_seq_qual_data(const char *name, GtHcrEncoder *hcr_enc,
                                   GtTimer *timer, GtError *err)
{
  int had_err = 0;
  FILE *fp;
  GtUword dummy = 0;
  GtWord pos;

  gt_error_check(err);

  fp = gt_fa_fopen_with_suffix(name, HCRFILESUFFIX, "wb", err);
  if (fp == NULL)
    had_err = -1;

  if (!had_err) {
    if (timer != NULL)
      gt_timer_show_progress(timer, "write sequences and qualities encoding",
                             stdout);
    hcr_write_file_info(fp, hcr_enc);

    had_err = hcr_write_seqdistrtab(fp, hcr_enc);

    if (!had_err) {
      bool is_not_at_pageborder;

      pos = ftell(fp);
      gt_xfwrite_one(&dummy, fp);

      is_not_at_pageborder = (ftell(fp) % hcr_enc->pagesize) != 0;

      if (is_not_at_pageborder)
        hcr_enc->seq_encoder->start_of_encoding =
          (ftell(fp) / hcr_enc->pagesize + 1) * hcr_enc->pagesize;
      else
        hcr_enc->seq_encoder->start_of_encoding = ftell(fp);

      if (hcr_enc->page_sampling)
        hcr_enc->seq_encoder->sampling =
          gt_sampling_new_page(hcr_enc->sampling_rate,
                               (off_t) hcr_enc->seq_encoder->start_of_encoding);
      else if (hcr_enc->regular_sampling)
        hcr_enc->seq_encoder->sampling =
          gt_sampling_new_regular(hcr_enc->sampling_rate,
                               (off_t) hcr_enc->seq_encoder->start_of_encoding);

      had_err = hcr_write_seqs(fp, hcr_enc, err);
    }
    if (!had_err) {
      gt_assert(fp);
      gt_xfseek(fp, pos, SEEK_SET);
      gt_xfwrite_one(&hcr_enc->seq_encoder->startofsamplingtab, fp);
    }
    gt_fa_xfclose(fp);
  }
  return 0;
}

static void hcr_read_file_info(GtHcrSeqDecoder *seq_dec, FILE *fp)
{
  GtUword i;
  GT_UNUSED size_t read,
            one = (size_t) 1;

  read = gt_xfread_one(&seq_dec->num_of_files, fp);
  gt_assert(read == one);

  seq_dec->fileinfos = gt_calloc((size_t) seq_dec->num_of_files,
                                          sizeof (*seq_dec->fileinfos));
  for (i = 0; i < seq_dec->num_of_files; i++) {
    read = gt_xfread_one(&seq_dec->fileinfos[i].readnum, fp);
    gt_assert(read == one);
    read = gt_xfread_one(&seq_dec->fileinfos[i].readlength, fp);
    gt_assert(read == one);
  }
  seq_dec->num_of_reads = seq_dec->fileinfos[seq_dec->num_of_files - 1].readnum;
}

static HcrHuffDataIterator *decoder_init_data_iterator(
                                                GtWord start_of_encoding,
                                                GtWord end_of_encoding,
                                                const GtStr *filename)
{
  HcrHuffDataIterator *data_iter = gt_malloc(sizeof (*data_iter));
  data_iter->path = gt_str_get(filename);
  data_iter->pos = data_iter->start = (size_t) start_of_encoding;
  data_iter->end = (size_t) end_of_encoding;
  data_iter->pages_per_chunk = HCR_PAGES_PER_CHUNK;
  data_iter->pagesize = gt_pagesize();
  gt_assert(data_iter->start % data_iter->pagesize == 0);
  data_iter->blocksize = data_iter->pagesize * data_iter->pages_per_chunk;
  gt_safe_assign(data_iter->bitseq_per_chunk,
                 (data_iter->blocksize / sizeof (GtBitsequence)));
  data_iter->data = NULL;
  return data_iter;
}

static int get_next_file_chunk_for_huffman(GtBitsequence **bits,
                                           GtUword *length,
                                           GtUword *offset,
                                           GtUword *pad_length,
                                           void *meminfo)
{
  const int empty = 0,
            success = 1;
  HcrHuffDataIterator *data_iter;
  gt_assert(meminfo);
  gt_assert(bits && length && offset && pad_length);
  data_iter = (HcrHuffDataIterator*) meminfo;

  gt_log_log("pos in iter: "GT_WU"", (GtUword) data_iter->pos);
  if (data_iter->pos < data_iter->end) {
    gt_fa_xmunmap(data_iter->data);
    data_iter->data = NULL;
    data_iter->data = gt_fa_xmmap_read_range(data_iter->path,
                                        (size_t) data_iter->blocksize,
                                        data_iter->pos);
    data_iter->pos += data_iter->blocksize;
    if (data_iter->pos > data_iter->end) {
      gt_safe_assign(*length,
                    (data_iter->blocksize - (data_iter->pos - data_iter->end)));
      *length /= sizeof (GtBitsequence);
    }
    else
      *length = data_iter->bitseq_per_chunk;

    *offset = 0;
    *pad_length = 0;
    *bits = data_iter->data;
    return success;
  }
  gt_fa_xmunmap(data_iter->data);
  data_iter->data = NULL;
  *bits = NULL;
  *length = 0;
  *offset = 0;
  *pad_length = 0;
  return empty;
}

static void reset_data_iterator_to_pos(HcrHuffDataIterator *data_iter,
                                       size_t pos)
{
  gt_assert(data_iter);
  gt_assert(pos < data_iter->end);
  gt_assert(data_iter->start <= pos);
  gt_fa_xmunmap(data_iter->data);
  gt_log_log("reset to pos: "GT_WU"", (GtUword) pos);
  data_iter->data = NULL;
  data_iter->pos = pos;
}

static void reset_data_iterator_to_start(HcrHuffDataIterator *data_iter)
{
  gt_assert(data_iter);
  reset_data_iterator_to_pos(data_iter, data_iter->start);
}

static void data_iterator_delete(HcrHuffDataIterator *data_iter)
{
  if (data_iter != NULL)
    gt_fa_xmunmap(data_iter->data);
  gt_free(data_iter);
}

static GtRBTree *seq_decoder_init_file_info(FastqFileInfo *fileinfos,
                                            GtUword num_off_files)
{
  GtRBTree *tree;
  bool nodecreated = false;
  GtUword i;

  tree = gt_rbtree_new(hcr_cmp_FastqFileInfo, NULL, NULL);
  for (i = 0; i < num_off_files; i++) {
    (void) gt_rbtree_search(tree, &fileinfos[i], &nodecreated);
  }
  gt_assert(nodecreated);
  return tree;
}

static inline GtWord decoder_calc_start_of_encoded_data(FILE *fp)
{
  bool is_not_at_pageborder;
  GtUword pagesize = gt_pagesize();

  is_not_at_pageborder = (ftell(fp) % pagesize) != 0;

  if (is_not_at_pageborder)
    return (ftell(fp) / pagesize + 1) * pagesize;
  else
    return ftell(fp);
}

static void seq_decoder_init_huffman(GtHcrSeqDecoder *seq_dec,
                                     GtWord end_of_encoding,
                                     GtBaseQualDistr *bqd,
                                     GtError *err)
{
  seq_dec->data_iter = decoder_init_data_iterator(seq_dec->start_of_encoding,
                                                  end_of_encoding,
                                                  seq_dec->filename);

  gt_assert(seq_dec->data_iter);
  seq_dec->huffman = gt_huffman_new(bqd, hcr_base_qual_distr_func,
                                    (GtUword) bqd->ncols * bqd->nrows);

  seq_dec->huff_dec = gt_huffman_decoder_new_from_memory(seq_dec->huffman,
                                                get_next_file_chunk_for_huffman,
                                                seq_dec->data_iter, err);
}

static void hcr_seq_decoder_delete(GtHcrSeqDecoder *seq_dec)
{
  if (seq_dec != NULL) {
    gt_free(seq_dec->fileinfos);
    gt_huffman_decoder_delete(seq_dec->huff_dec);
    gt_huffman_delete(seq_dec->huffman);
    gt_sampling_delete(seq_dec->sampling);
    gt_rbtree_delete(seq_dec->file_info_rbt);
    gt_str_delete(seq_dec->filename);
    data_iterator_delete(seq_dec->data_iter);
    gt_array_delete(seq_dec->symbols);
    gt_free(seq_dec);
  }
}

static GtHcrSeqDecoder *hcr_seq_decoder_new(GtAlphabet *alpha, const char *name,
                                            GtError *err)
{
  GtHcrSeqDecoder *seq_dec = gt_malloc(sizeof (GtHcrSeqDecoder));
  GtBaseQualDistr *bqd = NULL;
  GtWord end_enc_start_sampling = 0;
  FILE *fp = NULL;
  GT_UNUSED size_t read,
            one = (size_t) 1;

  seq_dec->alpha = alpha;
  seq_dec->alphabet_size = gt_alphabet_size(alpha);
  seq_dec->cur_read = 0;
  seq_dec->data_iter = NULL;
  seq_dec->file_info_rbt = NULL;
  seq_dec->fileinfos = NULL;
  seq_dec->filename = gt_str_new_cstr(name);
  seq_dec->huff_dec = NULL;
  seq_dec->huffman = NULL;
  seq_dec->sampling = NULL;
  seq_dec->symbols = NULL;
  gt_str_append_cstr(seq_dec->filename, HCRFILESUFFIX);

  fp = gt_fa_fopen_with_suffix(name, HCRFILESUFFIX, "rb", err);
  if (gt_error_is_set(err)) {
    hcr_seq_decoder_delete(seq_dec);
    seq_dec = NULL;
  }
  else {
    hcr_read_file_info(seq_dec, fp);

    bqd = hcr_base_qual_distr_new_from_file(fp, seq_dec->alpha);
    seq_dec->qual_offset = bqd->qual_offset;

    read = gt_xfread_one(&end_enc_start_sampling, fp);
    gt_assert(read == one);

    seq_dec->start_of_encoding = decoder_calc_start_of_encoded_data(fp);

    seq_decoder_init_huffman(seq_dec, end_enc_start_sampling, bqd, err);
    if (gt_error_is_set(err)) {
      hcr_seq_decoder_delete(seq_dec);
      seq_dec = NULL;
    }
  }

  if (seq_dec != NULL) {
    gt_xfseek(fp, end_enc_start_sampling, SEEK_SET);
    seq_dec->sampling = gt_sampling_read(fp);

    seq_dec->file_info_rbt = seq_decoder_init_file_info(seq_dec->fileinfos,
                                                        seq_dec->num_of_files);
  }

  hcr_base_qual_distr_delete(bqd);
  gt_fa_fclose(fp);
  return seq_dec;
}

GtHcrDecoder *gt_hcr_decoder_new(const char *name, GtAlphabet *alpha,
                                 bool descs, GtTimer *timer, GtError *err)
{
  GtHcrDecoder *hcr_dec;
  int had_err = 0;

  gt_error_check(err);
  if (timer != NULL)
    gt_timer_show_progress(timer, "initialize hcr decoder", stdout);

  hcr_dec = gt_malloc(sizeof (GtHcrDecoder));

  if (descs) {
    hcr_dec->encdesc = gt_encdesc_load(name, err);
    if (gt_error_is_set(err)) {
      had_err = -1;
    }
  }
  else
    hcr_dec->encdesc = NULL;

  if (!had_err) {
    hcr_dec->seq_dec = hcr_seq_decoder_new(alpha, name, err);
    if (!gt_error_is_set(err))
      return hcr_dec;
  }
  gt_hcr_decoder_delete(hcr_dec);
  return NULL;
}

bool gt_hcr_decoder_has_descs_support(GtHcrDecoder *hcr_dec)
{
  gt_assert(hcr_dec);
  if (!hcr_dec->encdesc)
    return false;
  else
    return true;
}

static inline char get_qual_from_symbol(GtHcrSeqDecoder *seq_dec,
                                        GtUword symbol)
{
  return (char) (symbol / seq_dec->alphabet_size + seq_dec->qual_offset);
}

static unsigned char get_base_from_symbol(GtHcrSeqDecoder *seq_dec,
                                 GtUword symbol)
{
  unsigned char base = (unsigned char) symbol % seq_dec->alphabet_size;
  if (base == (unsigned char) seq_dec->alphabet_size)
    return (unsigned char) WILDCARD;
  return base;
}

static int hcr_next_seq_qual(GtHcrSeqDecoder *seq_dec, char *seq, char *qual,
                             GtError *err)
{
  enum state {
    HCR_ERROR = -1,
    END,
    SUCCESS
  };
  unsigned char base;
  GtUword i,
                nearestsample,
                *symbol;
  size_t startofnearestsample = 0;
  enum state status = END;
  FastqFileInfo cur_read;
  FastqFileInfo *fileinfo = NULL;

  if (seq_dec->cur_read <= seq_dec->num_of_reads) {
    status = SUCCESS;
    if (seq_dec->symbols == NULL)
      seq_dec->symbols = gt_array_new(sizeof (GtUword));
    else
      gt_array_reset(seq_dec->symbols);

    cur_read.readnum = seq_dec->cur_read;
    gt_log_log("cur_read: "GT_WU"",seq_dec->cur_read);
    fileinfo = (FastqFileInfo *)gt_rbtree_next_key(seq_dec->file_info_rbt,
                                                   &cur_read,
                                                   hcr_cmp_FastqFileInfo,
                                                   NULL);
    gt_assert(fileinfo);

    /* reset huffman_decoder if next read is sampled */
    if (gt_sampling_get_next_elementnum(seq_dec->sampling) ==
          seq_dec->cur_read) {
      gt_log_log("reset because sampled read is next");
      (void) gt_sampling_get_next_sample(seq_dec->sampling,
                                         &nearestsample,
                                         &startofnearestsample);
      reset_data_iterator_to_pos(seq_dec->data_iter, startofnearestsample);
      (void) gt_huffman_decoder_get_new_mem_chunk(seq_dec->huff_dec, err);
      if (gt_error_is_set(err))
        status = HCR_ERROR;
    }
    if (status != HCR_ERROR) {
      int ret;
      ret =  gt_huffman_decoder_next(seq_dec->huff_dec, seq_dec->symbols,
                                     fileinfo->readlength, err);
      if (ret != 1)
        status = HCR_ERROR;
      if (ret == 0)
        gt_error_set(err, "reached end of file");
    }
    if (qual || seq) {
      gt_log_log("set strings");
      for (i = 0; i < gt_array_size(seq_dec->symbols); i++) {
        symbol = (GtUword*) gt_array_get(seq_dec->symbols, i);
        if (qual != NULL)
          qual[i] = get_qual_from_symbol(seq_dec, *symbol);
        if (seq != NULL) {
          base = get_base_from_symbol(seq_dec, *symbol);
          seq[i] = (char)toupper(gt_alphabet_decode(seq_dec->alpha,
                                                    (GtUchar) base));
        }
      }
      if (qual != NULL)
        qual[gt_array_size(seq_dec->symbols)] = '\0';
      if (seq != NULL)
        seq[gt_array_size(seq_dec->symbols)] = '\0';
    }
    seq_dec->cur_read++;
  }
  return (int) status;
}

int gt_hcr_decoder_decode(GtHcrDecoder *hcr_dec, GtUword readnum,
                          char *seq, char *qual, GtStr *desc, GtError *err)
{
  GtUword nearestsample = 0,
                reads_to_read = 0,
                idx,
                current_read = hcr_dec->seq_dec->cur_read ;
  size_t startofnearestsample = 0;
  GtSampling *sampling;
  HcrHuffDataIterator *data_iter;
  GtHuffmanDecoder *huff_dec;

  gt_error_check(err);
  gt_assert(hcr_dec);
  gt_assert(readnum < hcr_dec->seq_dec->num_of_reads);
  gt_assert(seq != NULL && qual != NULL);

  if (current_read == readnum) {
    if (hcr_next_seq_qual(hcr_dec->seq_dec, seq, qual, err) == -1) {
      gt_assert(gt_error_is_set(err));
      return -1;
    }
  }
  else {
    sampling = hcr_dec->seq_dec->sampling;
    data_iter = hcr_dec->seq_dec->data_iter;
    huff_dec = hcr_dec->seq_dec->huff_dec;

    if (sampling != NULL) {
      gt_sampling_get_page(sampling,
                           readnum,
                           &nearestsample,
                           &startofnearestsample);
      /* nearestsample <= cur_read < readnum: current sample is the right one */
      if (nearestsample <= current_read && current_read <= readnum)
        reads_to_read = readnum - current_read;
      else { /* reset decoder to new sample */
        reset_data_iterator_to_pos(data_iter, startofnearestsample);
        (void) gt_huffman_decoder_get_new_mem_chunk(huff_dec, err);
        if (gt_error_is_set(err))
          return -1;
        reads_to_read = readnum - nearestsample;
        hcr_dec->seq_dec->cur_read = nearestsample;
      }
      gt_log_log("reads to read: "GT_WU", nearest sample: "GT_WU"",
                 reads_to_read,nearestsample);
      gt_log_log("start of nearest: "GT_WU"", (GtUword) startofnearestsample);
    }
    else {
      if (current_read <= readnum)
        reads_to_read = readnum - current_read;
      else {
        reset_data_iterator_to_start(data_iter);
        (void) gt_huffman_decoder_get_new_mem_chunk(huff_dec, err);
        if (gt_error_is_set(err))
          return -1;
        reads_to_read = readnum;
        hcr_dec->seq_dec->cur_read = 0;
      }
    }

    for (idx = 0; idx < reads_to_read; idx++) {
      if (hcr_next_seq_qual(hcr_dec->seq_dec, seq,qual, err) == -1) {
        gt_assert(gt_error_is_set(err));
        return -1;
      }
      gt_log_log("seq:\n%s\nqual:\n%s", seq, qual);
    }

    if (hcr_next_seq_qual(hcr_dec->seq_dec, seq, qual, err) == -1) {
      gt_assert(gt_error_is_set(err));
      return -1;
    }
  }

  if (hcr_dec->encdesc != NULL) {
    if (gt_encdesc_decode(hcr_dec->encdesc, readnum, desc, err) == -1) {
      gt_error_set(err, "cannot retrieve description with number " GT_WU "."
                   "(%d)", readnum, __LINE__);
      return -1;
    }
  }
  return 0;
}

int gt_hcr_decoder_decode_range(GtHcrDecoder *hcr_dec, const char *name,
                                GtUword start, GtUword end,
                                GtTimer *timer, GtError *err)
{
  char qual[BUFSIZ] = {0},
       seq[BUFSIZ] = {0};
  GtStr *desc = gt_str_new();
  int had_err = 0;
  GtUword cur_width,
                cur_read;
  size_t i;
  FILE *output;
  GT_UNUSED GtHcrSeqDecoder *seq_dec;

  gt_error_check(err);
  gt_assert(hcr_dec && name);
  seq_dec = hcr_dec->seq_dec;
  gt_assert(start <= end);
  gt_assert(start < seq_dec->num_of_reads && end < seq_dec->num_of_reads);
  if (timer != NULL)
    gt_timer_show_progress(timer, "decode hcr", stdout);
  output = gt_fa_fopen_with_suffix(name, HCRFILEDECODEDSUFFIX, "w", err);
  if (output == NULL)
    had_err = -1;

  for (cur_read = start; had_err == 0 && cur_read <= end; cur_read++) {
    if (gt_hcr_decoder_decode(hcr_dec, cur_read, seq, qual, desc, err) != 0)
      had_err = -1;
    else {
      gt_xfputc(HCR_DESCSEPSEQ, output);
      if (hcr_dec->encdesc != NULL)
        gt_xfputs(gt_str_get(desc), output);
      else
        fprintf(output, ""GT_WU"", cur_read);
      gt_xfputc('\n', output);
      for (i = 0, cur_width = 0; i < strlen(seq); i++, cur_width++) {
        if (cur_width == HCR_LINEWIDTH) {
          cur_width = 0;
          gt_xfputc('\n', output);
        }
        gt_xfputc(seq[i], output);
      }
      gt_xfputc('\n', output);
      gt_xfputc(HCR_DESCSEPQUAL, output);
      gt_xfputc('\n', output);
      for (i = 0, cur_width = 0; i < strlen(qual); i++, cur_width++) {
        if (cur_width == HCR_LINEWIDTH) {
          cur_width = 0;
          gt_xfputc('\n', output);
        }
        gt_xfputc(qual[i], output);
      }
      gt_xfputc('\n', output);
    }
  }
  gt_fa_xfclose(output);
  gt_str_delete(desc);
  return had_err;
}

GtUword gt_hcr_decoder_num_of_reads(GtHcrDecoder *hcr_dec)
{
  gt_assert(hcr_dec);
  return hcr_dec->seq_dec->num_of_reads;
}

GtUword gt_hcr_decoder_readlength(GtHcrDecoder *hcr_dec,
                                        GtUword filenum)
{
  gt_assert(hcr_dec);
  gt_assert(hcr_dec->seq_dec->num_of_files > filenum);
  return hcr_dec->seq_dec->fileinfos[filenum].readlength;
}

GtHcrEncoder *gt_hcr_encoder_new(GtStrArray *files, GtAlphabet *alpha,
                                 bool descs, GtQualRange qrange, GtTimer *timer,
                                 GtError *err)
{
  GtBaseQualDistr *bqd;
  GtHcrEncoder *hcr_enc;
  GtSeqIterator *seqit;
  GtStrArray *file;
  int had_err = 0,
      status;
  GtUword len1,
                len2,
                i,
                num_of_reads = 0;
  const GtUchar *seq,
                *qual;
  char *desc;

  gt_error_check(err);
  gt_assert(alpha && files);

  if (timer != NULL)
    gt_timer_show_progress(timer, "get <base,qual> distr", stdout);

  if (qrange.start != GT_UNDEF_UINT)
    if (qrange.start == qrange.end) {
      gt_error_set(err, "qrange.start must unequal qrange.end");
      return NULL;
    }

  hcr_enc = gt_malloc(sizeof (GtHcrEncoder));
  hcr_enc->files = files;
  hcr_enc->num_of_files = gt_str_array_size(files);
  hcr_enc->num_of_reads = 0;
  hcr_enc->page_sampling = false;
  hcr_enc->regular_sampling = false;
  hcr_enc->sampling_rate = 0;
  hcr_enc->pagesize = gt_pagesize();
  if (descs) {
    hcr_enc->encdesc_encoder = gt_encdesc_encoder_new();
    if (timer != NULL)
      gt_encdesc_encoder_set_timer(hcr_enc->encdesc_encoder, timer);
  }
  else
    hcr_enc->encdesc_encoder = NULL;

  hcr_enc->seq_encoder = gt_malloc(sizeof (GtHcrSeqEncoder));
  hcr_enc->seq_encoder->alpha = alpha;
  hcr_enc->seq_encoder->sampling = NULL;
  hcr_enc->seq_encoder->fileinfos = gt_calloc((size_t) hcr_enc->num_of_files,
                                   sizeof (*(hcr_enc->seq_encoder->fileinfos)));
  hcr_enc->seq_encoder->qrange = qrange;
  bqd = hcr_base_qual_distr_new(alpha, qrange);

  /* check if reads in the same file are of same length and get
     <base, quality> pair distribution */
  for (i = 0; i < hcr_enc->num_of_files; i++) {
    file = gt_str_array_new();
    gt_str_array_add(file, gt_str_array_get_str(files, i));
    seqit = gt_seq_iterator_fastq_new(file, err);
    if (!seqit) {
      gt_error_set(err, "cannot initialize GtSeqIteratorFastQ object");
      had_err = -1;
    }
    if (!had_err) {
      gt_seq_iterator_set_symbolmap(seqit, gt_alphabet_symbolmap(alpha));
      gt_seq_iterator_set_quality_buffer(seqit, &qual);
      status = gt_seq_iterator_next(seqit, &seq, &len1, &desc, err);

      if (status == 1) {
        num_of_reads = 1UL;
        while (!had_err) {
          status = gt_seq_iterator_next(seqit, &seq, &len2, &desc, err);
          if (status == -1)
            had_err = -1;
          if (status != 1)
            break;
          if (len2 != len1) {
            gt_error_set(err, "reads have to be of equal length");
            had_err = -1;
            break;
          }
          if (hcr_base_qual_distr_add(bqd, qual, seq, len1) != 0)
            had_err = -1;
          len1 = len2;
          num_of_reads++;
        }
      }
      else if (status == -1)
        had_err = -1;

      if (!had_err) {
        if (i == 0)
          hcr_enc->seq_encoder->fileinfos[i].readnum = num_of_reads;
        else
          hcr_enc->seq_encoder->fileinfos[i].readnum =
            hcr_enc->seq_encoder->fileinfos[i - 1].readnum + num_of_reads;
        hcr_enc->seq_encoder->fileinfos[i].readlength = len1;
      }
    }
    hcr_enc->num_of_reads += num_of_reads;
    gt_str_array_delete(file);
    gt_seq_iterator_delete(seqit);
  }
  if (!had_err)
    hcr_base_qual_distr_trim(bqd);

  if (!had_err) {
    if (timer != NULL)
      gt_timer_show_progress(timer, "build huffman tree for sequences and"
                             " qualities", stdout);
    hcr_enc->seq_encoder->huffman =
      gt_huffman_new(bqd,
                     hcr_base_qual_distr_func,
                     (GtUword) bqd->ncols * bqd->nrows);
  }
  if (!had_err) {
    hcr_enc->seq_encoder->qual_offset = bqd->qual_offset;
    hcr_base_qual_distr_delete(bqd);
    return hcr_enc;
  }
  return NULL;
}

bool gt_hcr_encoder_has_descs_support(GtHcrEncoder *hcr_enc)
{
  gt_assert(hcr_enc);
  if (!hcr_enc->encdesc_encoder)
    return false;
  else
    return true;
}

void gt_hcr_encoder_set_sampling_none(GtHcrEncoder *hcr_enc)
{
  gt_assert(hcr_enc);
  hcr_enc->page_sampling = hcr_enc->regular_sampling = false;
  if (hcr_enc->encdesc_encoder != NULL)
    gt_encdesc_encoder_set_sampling_none(hcr_enc->encdesc_encoder);
}

void gt_hcr_encoder_set_sampling_regular(GtHcrEncoder *hcr_enc)
{
  gt_assert(hcr_enc);
  hcr_enc->page_sampling = false;
  hcr_enc->regular_sampling = true;
  if (hcr_enc->encdesc_encoder != NULL)
    gt_encdesc_encoder_set_sampling_regular(hcr_enc->encdesc_encoder);
}

void gt_hcr_encoder_set_sampling_page(GtHcrEncoder *hcr_enc)
{
  gt_assert(hcr_enc);
  hcr_enc->page_sampling = true;
  hcr_enc->regular_sampling = false;
  if (hcr_enc->encdesc_encoder != NULL)
    gt_encdesc_encoder_set_sampling_page(hcr_enc->encdesc_encoder);
}

void gt_hcr_encoder_set_sampling_rate(GtHcrEncoder *hcr_enc,
                                      GtUword srate)
{
  gt_assert(hcr_enc);
  hcr_enc->sampling_rate = srate;

  if (hcr_enc->encdesc_encoder != NULL)
    gt_encdesc_encoder_set_sampling_rate(hcr_enc->encdesc_encoder, srate);
}

GtUword gt_hcr_encoder_get_sampling_rate(GtHcrEncoder *hcr_enc)
{
  gt_assert(hcr_enc);
  return hcr_enc->sampling_rate;
}

int gt_hcr_encoder_encode(GtHcrEncoder *hcr_enc, const char *name,
                          GtTimer *timer, GtError *err)
{
  int had_err = 0;
  GtStr *name1;
  gt_error_check(err);
  if (timer != NULL)
    gt_timer_show_progress(timer, "write encoding", stdout);
  if (hcr_enc->encdesc_encoder != NULL) {
    GtCstrIterator *cstr_iterator = gt_fasta_header_iterator_new(hcr_enc->files,
                                                                 err);
    had_err = gt_encdesc_encoder_encode(hcr_enc->encdesc_encoder,
                                        cstr_iterator, name, err);
    gt_cstr_iterator_delete(cstr_iterator);
  }

  if (!had_err)
    had_err = hcr_write_seq_qual_data(name, hcr_enc, timer, err);

  if (!had_err && gt_log_enabled()) {
    name1 = gt_str_new_cstr(name);
    gt_str_append_cstr(name1, HCRFILESUFFIX);
    gt_log_log("sequences with qualities encoding overview:");
    gt_log_log("**>");
    if (hcr_enc->page_sampling)
        gt_log_log("applied sampling technique: sampling every "GT_WU"th page",
                   hcr_enc->sampling_rate);
    else if (hcr_enc->regular_sampling)
        gt_log_log("applied sampling technique: sampling every "GT_WU"th read",
                   hcr_enc->sampling_rate);
    else
        gt_log_log("applied sampling technique: none");

    gt_log_log("total number of encoded nucleotide sequences with qualities: "
               ""GT_WU"", hcr_enc->num_of_reads);
    gt_log_log("total number of encoded nucleotides: "GT_LLU"",
               hcr_enc->seq_encoder->total_num_of_symbols);
    gt_log_log("bits per nucleotide encoding: %f",
               (gt_file_estimate_size(gt_str_get(name1)) * 8.0) /
                 hcr_enc->seq_encoder->total_num_of_symbols);
    gt_log_log("<**");
    gt_str_delete(name1);
  }
  return had_err;
}

void gt_hcr_encoder_delete(GtHcrEncoder *hcr_enc)
{
  if (!hcr_enc)
    return;
  gt_huffman_delete(hcr_enc->seq_encoder->huffman);
  gt_sampling_delete(hcr_enc->seq_encoder->sampling);
  gt_free(hcr_enc->seq_encoder->fileinfos);
  gt_free(hcr_enc->seq_encoder);
  gt_encdesc_encoder_delete(hcr_enc->encdesc_encoder);
  gt_free(hcr_enc);
}

void gt_hcr_decoder_delete(GtHcrDecoder *hcr_dec)
{
  if (hcr_dec != NULL) {
    hcr_seq_decoder_delete(hcr_dec->seq_dec);
    gt_encdesc_delete(hcr_dec->encdesc);
    gt_free(hcr_dec);
  }
}
