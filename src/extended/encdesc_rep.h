/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCDESC_REP_H
#define ENCDESC_REP_H

#include "core/arraydef.h"
#include "core/bittab_api.h"
#include "core/disc_distri_api.h"
#include "core/hashmap-generic.h"
#include "core/hashtable.h"
#include "core/intbits.h"
#include "extended/bitinstream.h"
#include "extended/bitoutstream.h"
#include "extended/encdesc.h"
#include "extended/huffcode.h"
#include "extended/sampling.h"

#define GT_ENCDESC_MAX_NUM_VAL_HUF 512UL

DECLARE_HASHMAP(long int, li, unsigned long long, ull, static, inline)
DEFINE_HASHMAP(long int, li, unsigned long long, ull, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

typedef struct {
  GtHashtable *li_ull_hashmap;
  long correction_base;
} EncdescHuffDist;

typedef struct {
  char *data,
       sep;
  bool is_numeric,
       is_cons,
       fieldlen_is_const,
       use_delta_coding,
       is_delta_cons,
       is_value_cons,
       use_hc,
       has_zero_padding;
  unsigned bits_per_num,
           bits_per_value,
           bits_per_len,
           max_zero;
  long max_delta,
       min_delta,
       global_delta,
       min_value,
       max_value,
       global_value,
       prev_value;
  unsigned long max_len,
                min_len,
                num_values_size,
                len,
                delta_values_size;
  GtBittab *bittab;
  GtHashtable *num_values,
              *delta_values,
              **chars;
  /* TODO: test if this can be a GtHashtable, too */
  GtDiscDistri *zero_count;
  GtHuffman **huffman_chars,
            *huffman_num,
            *huffman_zero_count;
} DescField;

struct GtEncdesc {
  bool num_of_fields_is_cons;
  unsigned bits_per_field;
  long start_of_encoding,
       start_of_samplingtab;
  unsigned long num_of_descs,
                num_of_fields,
                cur_desc;
  unsigned long long total_num_of_chars;
  DescField *fields;
  GtBitInStream *bitinstream;
  GtArrayGtUlong num_of_fields_tab;
  GtSampling *sampling;
  long pagesize;
};

struct GtEncdescEncoder {
  GtTimer *timer;
  GtEncdesc *encdesc;
  bool regular_sampling,
       page_sampling;
  unsigned long sampling_rate;
};

typedef struct {
  GtBitsequence code;
  unsigned length;
} EncdescCode;

GT_DECLAREARRAYSTRUCT(EncdescCode);

typedef struct {
  GtArrayEncdescCode *codes;
  unsigned long total_bits_prepared,
                bits_to_write,
                cur_field_start_pos,
                cur_field_num,
                cur_desc;
  char *descbuffer;
  bool sample;
} EncdescWriteInfo;

#endif
