/*
  Copyright (c) 2012 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
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
#include <errno.h>
#ifndef S_SPLINT_S
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/bittab_api.h"
#include "core/compat.h"
#include "core/cstr_api.h"
#include "core/disc_distri_api.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/hashmap-generic.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/parseutils.h"
#include "core/undef_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "extended/bitinstream.h"
#include "extended/bitoutstream.h"
#include "extended/cstr_iterator.h"
#include "extended/encdesc.h"
#include "extended/encdesc_header_io.h"
#include "extended/encdesc_rep.h"
#include "extended/huffcode.h"
#include "extended/rbtree.h"
#include "extended/sampling.h"

#define GT_ENCDESC_ARRAY_RESIZE 50
#define GT_ENCDESC_FILESUFFIX ".ede"
#define GT_ENCDESC_NUMOFSEPS 10UL
#define GT_ENCDESC_SEPS '.', '_', ',', '=', ':', '/' , '-', '|', ' ', '\0'

static inline void encdesc_prepare_num_of_fields(GtEncdesc *encdesc,
                                                 EncdescWriteInfo *info);

static void parse_number_out_of_current_field(EncdescWriteInfo *info,
                                              GtWord *retval);

static void numeric_field_prepare_zero_padding(GtEncdesc *encdesc,
                                               EncdescWriteInfo *info);

static void numeric_field_sample_prepare_verbose_value(GtEncdesc *encdesc,
                                                       EncdescWriteInfo *info,
                                                       GtUword value);

static void numeric_field_prepare_huffman_value(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword value);

static void numeric_field_prepare_verbose_value(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword value);

static inline void prepare_numeric_field(GtEncdesc *encdesc,
                                         EncdescWriteInfo *info);

static void inline regular_field_prepare_length(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword endpos);

static void inline regular_field_prepare_char(GtEncdesc *encdesc,
                                              EncdescWriteInfo *info,
                                              GtUword char_pos);

static inline void prepare_generic(EncdescWriteInfo *info,
                                   unsigned length,
                                   GtBitsequence code);

static void prepare_write_data_and_count_bits(GtEncdesc *encdesc,
                                              EncdescWriteInfo *info);

static inline void append_data_to_bitstream_and_reset_info(
                                                     EncdescWriteInfo *info,
                                                     GtBitOutStream *bitstream);

static inline void reset_info(EncdescWriteInfo *info);

static inline unsigned count_leading_zeros(const char *number);

static unsigned encdesc_digits_per_value(GtUword value,
                                         GtUword base)
{
  if (value > 0) {
    double log_val, log_base;
    log_val = log((double) value);
    log_base = log((double) base);
    return (unsigned) (floor(log_val / log_base )) + 1U;
  }
  else
    return 1U;
}

static int encdesc_hashmap_distr_add(GtHashtable *hm_distri, GtWord key)
{
  GtUint64 *valueptr;
  gt_assert(hm_distri != NULL);
  valueptr = li_ull_gt_hashmap_get(hm_distri, key);
  if (!valueptr) {
    li_ull_gt_hashmap_add(hm_distri, key, (GtUint64) 1);
    return 1;
  }
  else {
    (*valueptr) += 1;
    return 0;
  }
}

static int encdesc_analyze_descs(GtEncdesc *encdesc,
                                 GtCstrIterator *cstr_iterator,
                                 GtError *err)
{
  int status, had_err = 0;
  unsigned zero_count;
  GtUword cur_desc = 0,
                j_idx,
                k_idx,
                desc_char_idx,
                sep_idx,
                desclength,
                start_pos = 0,
                longest_desc_field_idx = 0,
                cur_field_num,
                tmp_numoffields,
                chars_len,
                len_diff;
  GtWord value,
       value_diff,
       out,
       value_delta;
  char *longest_desc = NULL,
       *mutable_desc = NULL,
       sep[GT_ENCDESC_NUMOFSEPS] = {GT_ENCDESC_SEPS};
  const char *descbuffer = NULL;
  bool /* equal, */
       found;
  DescField *cur_field;

  encdesc->num_of_fields = 0;
  encdesc->num_of_fields_is_cons = true;
  GT_INITARRAY(&encdesc->num_of_fields_tab, GtUlong);

  /* find description with maximum number of fields */
  /* TODO check constant number of fields here */
  /* TODO count totalnumofchars here */
  gt_error_check(err);
  status = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err);
  while (status != 0) {
    if (status < 0)
      return status;
    gt_assert(descbuffer != NULL);

    desclength = (GtUword) strlen(descbuffer);

    tmp_numoffields = 0;
    for (desc_char_idx = 0; desc_char_idx <= desclength; desc_char_idx++) {
      for (sep_idx = 0; sep_idx < GT_ENCDESC_NUMOFSEPS; sep_idx++) {
        if (descbuffer[desc_char_idx] == sep[sep_idx]) {
          if (desc_char_idx - start_pos > 0) {
            tmp_numoffields++;
            start_pos = desc_char_idx + 1;
          }
          break;
        }
      }
    }
    if (tmp_numoffields > encdesc->num_of_fields) {
      encdesc->num_of_fields = tmp_numoffields;
      if (longest_desc != NULL) {
        gt_free(longest_desc);
      }
      longest_desc = gt_cstr_dup(descbuffer);
    }
    cur_desc++;
    status = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err);
  }
  if (encdesc->num_of_fields == 0) {
    gt_error_set(err, "The file given seems to have no descriptions, there is "
                      "nothing to compress, aborting.");
    had_err = 1;
  }

  if (!had_err) {
    gt_assert(longest_desc != NULL);
    encdesc->fields = gt_calloc((size_t) encdesc->num_of_fields,
                                sizeof (*encdesc->fields));

    /* analyze description with maximum number of fields and initialize fields
     */
    desclength = (GtUword) strlen(longest_desc);
    start_pos = 0;

    for (desc_char_idx = 0; desc_char_idx <= desclength; desc_char_idx++) {
      found = false;
      for (sep_idx = 0; sep_idx < GT_ENCDESC_NUMOFSEPS; sep_idx++)
        if (longest_desc[desc_char_idx] == sep[sep_idx]) {
          if (desc_char_idx - start_pos > 0)
            found = true;
        }

      if (found) {
        cur_field = &encdesc->fields[longest_desc_field_idx];
        cur_field->sep = longest_desc[desc_char_idx];
        longest_desc[desc_char_idx] = '\0';
        cur_field->len = desc_char_idx - start_pos;
        cur_field->data = gt_cstr_dup(longest_desc + start_pos);

        cur_field->max_len = cur_field->len;
        cur_field->min_len = cur_field->len;
        cur_field->is_cons = true;
        cur_field->fieldlen_is_const = true;
        cur_field->has_zero_padding = false;
        cur_field->chars =
          gt_calloc((size_t) cur_field->max_len, sizeof (cur_field->chars));

        for (j_idx = 0; j_idx < cur_field->max_len; j_idx++)
          cur_field->chars[j_idx] = li_ull_gt_hashmap_new();

        cur_field->num_values = li_ull_gt_hashmap_new();
        cur_field->delta_values = li_ull_gt_hashmap_new();
        cur_field->zero_count = gt_disc_distri_new();
        cur_field->max_zero = 0;
        cur_field->num_values_size = 0;
        cur_field->delta_values_size = 0;
        if (gt_parse_long(&out, cur_field->data) == 0) {
          cur_field->is_numeric = true;
          cur_field->max_value = out;
          cur_field->min_value = out;
        }
        cur_field->bittab = gt_bittab_new(cur_field->len);
        gt_assert(cur_field->bittab != NULL);
        /* Set all bits to 1*/
        gt_bittab_complement(cur_field->bittab,
                             cur_field->bittab);

        start_pos = desc_char_idx + 1;
        longest_desc_field_idx++;
      }
    }
    gt_assert(encdesc->num_of_fields == longest_desc_field_idx);

    gt_free(longest_desc);

    if (gt_cstr_iterator_reset(cstr_iterator, err) != 0)
      had_err = 1;
  }

  if (!had_err) {
    cur_desc = 0;

    /* analyze all descriptions */
    status = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err);
    while (status != 0) {
      if (status == -1)
        return -1;
      gt_assert(descbuffer != NULL);

      gt_free(mutable_desc);
      desclength = (GtUword) strlen(descbuffer);
      encdesc->total_num_of_chars += desclength;
      mutable_desc = gt_cstr_dup(descbuffer);

      cur_field_num = 0;
      start_pos = 0;

      for (desc_char_idx = 0;
           desc_char_idx <= desclength &&
             cur_field_num < longest_desc_field_idx;
           desc_char_idx++) {
        cur_field = &encdesc->fields[cur_field_num];
        /* check for end of string, if there are less fields then in longest */
        if (mutable_desc[desc_char_idx] ==  cur_field->sep ||
            mutable_desc[desc_char_idx] == '\0') {
          mutable_desc[desc_char_idx] = '\0';

          chars_len = desc_char_idx - start_pos;
          if (chars_len > cur_field->max_len) {
            cur_field->is_cons = false;
            cur_field->fieldlen_is_const = false;
            cur_field->chars =
              gt_realloc(cur_field->chars,
                         (size_t) chars_len * sizeof (GtHashtable*));
            for (k_idx = cur_field->max_len; k_idx < chars_len; k_idx++)
              cur_field->chars[k_idx] = li_ull_gt_hashmap_new();

            cur_field->max_len = chars_len;
          }
          else if (chars_len < cur_field->min_len) {
            cur_field->is_cons = false;
            cur_field->fieldlen_is_const = false;
            cur_field->min_len = chars_len;
          }

          for (k_idx = 0; k_idx < chars_len; k_idx++) {
            (void) encdesc_hashmap_distr_add(cur_field->chars[k_idx],
                                      (GtWord) mutable_desc[start_pos + k_idx]);
          }

          if (cur_field->is_cons) {
            if (chars_len != cur_field->len) {
              cur_field->is_cons = false;
              cur_field->fieldlen_is_const = false;
            }
            else {
              if (strcmp(cur_field->data, mutable_desc + start_pos) != 0)
                cur_field->is_cons = false;
            }
          }

          if (cur_field->is_numeric) {
            if (gt_parse_long(&out, (mutable_desc + start_pos)) != 0)
              cur_field->is_numeric = false;
            else {

              value = out;

              zero_count = count_leading_zeros((mutable_desc + start_pos));
              if (zero_count > 0)
                cur_field->has_zero_padding = true;
              if (zero_count > cur_field->max_zero)
                cur_field->max_zero = zero_count;
              gt_disc_distri_add(cur_field->zero_count,
                                 (GtUword) zero_count);

              if (cur_desc == 0) {
                cur_field->global_value =
                cur_field->min_value =
                cur_field->max_value = value;
                cur_field->is_value_cons = true;
              }
              else {
                value_delta = value - cur_field->prev_value;
                if (value < cur_field->min_value) {
                  cur_field->is_value_cons = false;
                  cur_field->min_value = value;
                }
                else if (value > cur_field->max_value) {
                  cur_field->is_value_cons = false;
                  cur_field->max_value = value;
                }

                if (cur_desc == 1UL) {
                  cur_field->max_delta =
                  cur_field->min_delta =
                  cur_field->global_delta = value_delta;
                  cur_field->is_delta_cons = true;
                }
                else if (value_delta > cur_field->max_delta) {
                    cur_field->is_delta_cons = false;
                    cur_field->max_delta = value_delta;
                }
                else if (value_delta < cur_field->min_delta) {
                  cur_field->is_delta_cons = false;
                  cur_field->min_delta = value_delta;
                }

                if (cur_field->delta_values_size > 0 || cur_desc == 1UL) {
                  if (encdesc_hashmap_distr_add(cur_field->delta_values,
                                                value_delta) == 1)
                    cur_field->delta_values_size++;

                  if (cur_field->delta_values_size >
                        GT_ENCDESC_MAX_NUM_VAL_HUF) {
                    cur_field->delta_values_size = 0;
                  }
                }
              }
              if (cur_field->num_values_size > 0 || cur_desc == 0) {
                if (encdesc_hashmap_distr_add(cur_field->num_values,
                                              value) == 1)
                  cur_field->num_values_size++;
                if (cur_field->num_values_size > GT_ENCDESC_MAX_NUM_VAL_HUF) {
                  cur_field->num_values_size = 0;
                }
              }
              cur_field->prev_value = value;
            }
          }

          /* unmark non constant positions */
          if (!cur_field->is_cons) {
            for (k_idx = 0; k_idx < cur_field->len; k_idx++) {
              if (k_idx < chars_len) {
                if (cur_field->data[k_idx] != mutable_desc[k_idx + start_pos])
                  gt_bittab_unset_bit(cur_field->bittab, k_idx);
              }
              else
                gt_bittab_unset_bit(cur_field->bittab, k_idx);
            }

          }
          start_pos = desc_char_idx + 1;
          cur_field_num++;
        }
      }

      if ((cur_field_num != longest_desc_field_idx) ||
          (desc_char_idx != desclength))
        encdesc->num_of_fields_is_cons = false;

      GT_STOREINARRAY(&encdesc->num_of_fields_tab, GtUlong, 128 ,cur_field_num);

      /* TODO this heuristic is bad, if the 2nd field is missing, all following
         fields are missing too! */
      /* change this so a field can be absent in a single description */
      /* absent fields are non constant */
      for (;cur_field_num < encdesc->num_of_fields;cur_field_num++) {
        cur_field = &encdesc->fields[cur_field_num];
        cur_field->is_cons = false;
        cur_field->fieldlen_is_const = false;
        cur_field->is_numeric = false;
      }
      cur_desc++;
      status = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err);
    }
    gt_free(mutable_desc);
    mutable_desc = NULL;

    encdesc->num_of_descs = cur_desc;
    if (encdesc->num_of_fields_is_cons) {
      GT_FREEARRAY(&encdesc->num_of_fields_tab, GtUlong);
    }

    for (cur_field_num = 0;
         cur_field_num < encdesc->num_of_fields;
         cur_field_num++) {
      cur_field = &encdesc->fields[cur_field_num];
      gt_assert(cur_field->bittab != NULL);
      if (!cur_field->is_numeric) {
        if (!cur_field->is_cons) {
          len_diff = cur_field->max_len - cur_field->min_len;
          cur_field->bits_per_len =
            encdesc_digits_per_value((GtUword) len_diff, 2UL);
        }
      }
      else {
        /* TODO range can be large, but the size of the dist is more important
           for huffman coding */
        GtWord value_range, delta_range;
        value_range = labs(cur_field->max_value - cur_field->min_value);
        delta_range = labs(cur_field->max_delta - cur_field->min_delta);
        if (value_range < delta_range || cur_field->delta_values_size == 0) {
          cur_field->use_delta_coding = false;
          value_diff = value_range;
        }
        else {
          cur_field->use_delta_coding = true;
          value_diff = delta_range;
        }
        cur_field->bits_per_num =
          encdesc_digits_per_value((GtUword) value_diff, 2UL);
        cur_field->bits_per_value =
          encdesc_digits_per_value((GtUword) value_range, 2UL);
      }
    }
    encdesc->bits_per_field =
      encdesc_digits_per_value(encdesc->num_of_fields, 2UL);
  }
  return had_err;
}

static int encdesc_write_encoding(GtEncdesc *encdesc,
                                  GtCstrIterator *cstr_iterator,
                                  FILE *fp, GtError *err)
{
  int had_err = 0, striter_err;
  const char *descbuffer;
  GtUword desc_counter = 0,
                page_counter = 0,
                idx = 0,
                bits_left_in_page;
  /* const long pagesize = sysconf((int) _SC_PAGESIZE); */
  EncdescWriteInfo *info = gt_calloc((size_t) 1, sizeof (*info));
  GtBitOutStream *bitstream;

  bits_left_in_page = (GtUword) encdesc->pagesize * 8UL;

  info->codes = gt_malloc(sizeof (*info->codes));
  info->cur_desc = 0;
  GT_INITARRAY(info->codes, EncdescCode);

  /* gt_xfseek(fp, encdesc->start_of_encoding, SEEK_SET); */
  bitstream = gt_bitoutstream_new(fp);

  had_err = gt_cstr_iterator_reset(cstr_iterator, err);
  if (!had_err) {

    for (idx = 0; idx < encdesc->num_of_fields; idx++)
      encdesc->fields[idx].prev_value = 0;
    for (striter_err = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err);
         striter_err > 0;
         striter_err = gt_cstr_iterator_next(cstr_iterator, &descbuffer, err)) {
      gt_free(info->descbuffer);
      info->descbuffer = gt_cstr_dup(descbuffer);
      info->sample = false;

      info->total_bits_prepared = 0;

      prepare_write_data_and_count_bits(encdesc,
                                        info);

      /* check if a new sample has to be added */
      if (encdesc->sampling != NULL) {
        info->sample =
          gt_sampling_is_next_element_sample(encdesc->sampling,
                                             page_counter,
                                             desc_counter,
                                             info->total_bits_prepared,
                                             bits_left_in_page);
        if (info->sample) {
          /* sampled size and type of codes is different from unsampled */
          reset_info(info);
          prepare_write_data_and_count_bits(encdesc,
                                            info);
          gt_bitoutstream_flush_advance(bitstream);

          gt_sampling_add_sample(encdesc->sampling,
                                 (size_t) gt_bitoutstream_pos(bitstream),
                                 info->cur_desc);

          desc_counter = 0;
          page_counter = 0;
          bits_left_in_page = (GtUword) encdesc->pagesize * 8;
        }
      }

      while (bits_left_in_page < info->total_bits_prepared) {
        page_counter++;
        info->total_bits_prepared -= bits_left_in_page;
        bits_left_in_page = (GtUword) encdesc->pagesize * 8;
      }
      bits_left_in_page -= info->total_bits_prepared;
      /* always set first page as written */
      if (page_counter == 0)
        page_counter++;
      desc_counter++;
      info->cur_desc++;

      append_data_to_bitstream_and_reset_info(info, bitstream);
    }
    if (striter_err != 0) {
      had_err = striter_err;
      gt_assert(gt_error_is_set(err));
    }
  }

  if (!had_err) {
    gt_bitoutstream_flush(bitstream);
    encdesc->start_of_samplingtab = ftell(fp);

    if (encdesc->sampling != NULL)
      gt_sampling_write(encdesc->sampling, fp);
  }
  gt_bitoutstream_delete(bitstream);
  GT_FREEARRAY(info->codes, EncdescCode);
  gt_free(info->codes);
  gt_free(info->descbuffer);
  gt_free(info);
  return had_err;
}

static void prepare_write_data_and_count_bits(GtEncdesc *encdesc,
                                              EncdescWriteInfo *info)
{
  size_t desc_char_idx,
         field_char_idx,
         len;
  DescField *cur_field;

  len = strlen(info->descbuffer);
  info->cur_field_num = 0;
  info->cur_field_start_pos = 0;
  info->total_bits_prepared = 0;
  gt_assert(info->codes->nextfreeEncdescCode == 0);

  if (!encdesc->num_of_fields_is_cons) {
    encdesc_prepare_num_of_fields(encdesc,
                                  info);
  }

  cur_field = &encdesc->fields[info->cur_field_num];
  for (desc_char_idx = 0; desc_char_idx <= len; desc_char_idx++) {
    if (info->descbuffer[desc_char_idx] == cur_field->sep ||
        desc_char_idx == len) {

      info->descbuffer[desc_char_idx] = '\0';
      if (!cur_field->is_cons) {
        if (cur_field->is_numeric) {
          prepare_numeric_field(encdesc, info);
        }
        else {
          if (!cur_field->fieldlen_is_const) {
            regular_field_prepare_length(encdesc, info,
                                         (GtUword) desc_char_idx);
          }
          for (field_char_idx = 0;
               field_char_idx < desc_char_idx - info->cur_field_start_pos;
               field_char_idx++) {
            if (field_char_idx >= (size_t) cur_field->len ||
                !gt_bittab_bit_is_set(cur_field->bittab,
                                     (GtUword) field_char_idx)) {
              regular_field_prepare_char(encdesc, info,
                                         (GtUword) field_char_idx);
            }
          }

        }
      }
      info->cur_field_start_pos = (GtUword) desc_char_idx + 1;
      info->cur_field_num++;
      cur_field = &encdesc->fields[info->cur_field_num];
    }
  }
}

static inline void prepare_numeric_field(GtEncdesc *encdesc,
                                         EncdescWriteInfo *info)
{
  GtWord value = 0;
  GtUword to_store;
  DescField *cur_field = &encdesc->fields[info->cur_field_num];

  parse_number_out_of_current_field(info, &value);

  if (cur_field->has_zero_padding && !cur_field->fieldlen_is_const) {
    numeric_field_prepare_zero_padding(encdesc, info);
  }

  if (info->cur_desc == 0 || info->sample) {
    to_store = (GtUword) (value - cur_field->min_value);
    numeric_field_sample_prepare_verbose_value(encdesc, info, to_store);
  }
  else if (!cur_field->is_value_cons || !cur_field->is_delta_cons) {
    if (cur_field->use_delta_coding) {
      to_store = (GtUword)
        ((value - cur_field->prev_value) - cur_field->min_delta);
      gt_assert(to_store <= (GtUword) cur_field->max_delta -
                              cur_field->min_delta);
    }
    else {
      to_store = (GtUword) (value - cur_field->min_value);
      gt_assert(to_store <= (GtUword) cur_field->max_value -
                              cur_field->min_value);
    }

    if (cur_field->use_hc) {
      numeric_field_prepare_huffman_value(encdesc, info, to_store);
    }
    else {
      numeric_field_prepare_verbose_value(encdesc, info, to_store);
    }
  }
  cur_field->prev_value = value;
}

static inline void encdesc_prepare_num_of_fields(GtEncdesc *encdesc,
                                                 EncdescWriteInfo *info)
{
  GtBitsequence code;
  code =
    (GtBitsequence) encdesc->num_of_fields_tab.spaceGtUlong[info->cur_desc];
  prepare_generic(info,
                  (unsigned) encdesc->bits_per_field,
                  code);

}

static void parse_number_out_of_current_field(EncdescWriteInfo *info,
                                              GtWord *retval)
{
  char *cur_field_num = info->descbuffer + info->cur_field_start_pos;
  GT_UNUSED int had_err = gt_parse_long(retval, cur_field_num);
  gt_assert(!had_err);
}

static void numeric_field_prepare_zero_padding(GtEncdesc *encdesc,
                                               EncdescWriteInfo *info)
{
  GtBitsequence code;
  unsigned length,
           zero_count;
  GtHuffman *huffman;
  DescField *cur_field = &encdesc->fields[info->cur_field_num];

  zero_count = count_leading_zeros(info->descbuffer +
                                     info->cur_field_start_pos);

  huffman = cur_field->huffman_zero_count;
  gt_assert(huffman != NULL);
  gt_huffman_encode(huffman,
                    (GtUword) zero_count,
                    &code,
                    &length);
  prepare_generic(info, length, code);

}

static inline unsigned count_leading_zeros(const char *number)
{
  int idx;
  unsigned char count = 0;
  for (idx = 0; number[idx] == '0'; idx++) {
    count++;
  }
  if (count != 0 && number[idx] == '\0')
    count--;
  return count;
}

static void numeric_field_sample_prepare_verbose_value(GtEncdesc *encdesc,
                                                       EncdescWriteInfo *info,
                                                       GtUword value)
{
  DescField *cur_field = &encdesc->fields[info->cur_field_num];

  prepare_generic(info,
                  cur_field->bits_per_value,
                  (GtBitsequence) value);
}

static void numeric_field_prepare_huffman_value(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword value)
{
  GtBitsequence code;
  unsigned length;
  DescField *cur_field = &encdesc->fields[info->cur_field_num];
  GtHuffman *huffman;

  huffman = cur_field->huffman_num;
  gt_assert(huffman != NULL);
  gt_huffman_encode(huffman,
                    value,
                    &code,
                    &length);

  prepare_generic(info, length, code);
}

static void numeric_field_prepare_verbose_value(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword value)
{
  DescField *cur_field = &encdesc->fields[info->cur_field_num];

  prepare_generic(info,
                  cur_field->bits_per_num,
                  (GtBitsequence) value);

}

static void inline regular_field_prepare_length(GtEncdesc *encdesc,
                                                EncdescWriteInfo *info,
                                                GtUword endpos)
{
  DescField *field = &encdesc->fields[info->cur_field_num];
  GtBitsequence code = (GtBitsequence)
                       (endpos - info->cur_field_start_pos - field->min_len);
  prepare_generic(info, field->bits_per_len, code);
}

static void inline regular_field_prepare_char(GtEncdesc *encdesc,
                                              EncdescWriteInfo *info,
                                              GtUword char_pos)
{
  DescField *field = &encdesc->fields[info->cur_field_num];
  GtBitsequence code;
  unsigned length;
  GtHuffman *huffman = field->huffman_chars[char_pos];
  char cur_char = info->descbuffer[info->cur_field_start_pos + char_pos];

  gt_assert(huffman != NULL);
  gt_huffman_encode(huffman,
                    (GtUword) cur_char,
                    &code,
                    &length);
  prepare_generic(info, length, code);
}

static inline void prepare_generic(EncdescWriteInfo *info,
                                   unsigned length,
                                   GtBitsequence code)
{
  EncdescCode *this_code;

  GT_GETNEXTFREEINARRAY(this_code, info->codes, EncdescCode,
                        GT_ENCDESC_ARRAY_RESIZE);

  this_code->length = length;
  this_code->code = code;

  info->total_bits_prepared += this_code->length;
}

static inline void
append_data_to_bitstream_and_reset_info(EncdescWriteInfo *info,
                                        GtBitOutStream *bitstream)
{
  GtUword number_of_codes, idx;
  EncdescCode *codes = info->codes->spaceEncdescCode;

  number_of_codes = info->codes->nextfreeEncdescCode;

  for (idx = 0; idx < number_of_codes; idx++)
  {
    gt_bitoutstream_append(bitstream, codes[idx].code, codes[idx].length);
  }
  reset_info(info);
}

static inline void reset_info(EncdescWriteInfo *info)
{
  info->codes->nextfreeEncdescCode = 0;
  info->total_bits_prepared = 0;
}

static GtEncdesc *encdesc_new(void)
{
  GtEncdesc *encdesc;
  encdesc = gt_calloc((size_t) 1, sizeof (GtEncdesc));
  encdesc->bitinstream = NULL;
  encdesc->num_of_descs = 0;
  encdesc->total_num_of_chars = 0;
  encdesc->sampling = NULL;
  encdesc->pagesize = gt_pagesize();
  return encdesc;
}

GtEncdescEncoder* gt_encdesc_encoder_new(void)
{
  GtEncdescEncoder *ee;

  ee = gt_calloc((size_t) 1, sizeof (GtEncdescEncoder));
  ee->timer = NULL;
  ee->page_sampling = false;
  ee->regular_sampling = false;
  ee->sampling_rate = 0;

  ee->encdesc = encdesc_new();
  return ee;
}

void gt_encdesc_encoder_set_timer(GtEncdescEncoder *ee, GtTimer *timer)
{
  gt_assert(ee);
  ee->timer = timer;
}

GtTimer* gt_encdesc_encoder_get_timer(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  return ee->timer;
}

void gt_encdesc_encoder_set_sampling_none(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  ee->page_sampling = ee->regular_sampling = false;
}

void gt_encdesc_encoder_set_sampling_regular(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  ee->regular_sampling = true;
  ee->page_sampling = false;
}

void gt_encdesc_encoder_set_sampling_page(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  ee->page_sampling = true;
  ee->regular_sampling = false;
}

bool gt_encdesc_sampling_is_regular(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  return ee->regular_sampling;
}

bool gt_encdesc_sampling_is_page(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  return ee->page_sampling;
}

void gt_encdesc_encoder_set_sampling_rate(GtEncdescEncoder *ee,
                                          GtUword sampling_rate)
{
  gt_assert(ee);
  ee->sampling_rate = sampling_rate;
}

GtUword gt_encdesc_encoder_get_sampling_rate(GtEncdescEncoder *ee)
{
  gt_assert(ee);
  return ee->sampling_rate;
}

GtUint64 encdesc_hashmap_distr_get_corrected(const void *data,
                                                       GtUword key)
{
  EncdescHuffDist *dist = (EncdescHuffDist*) data;
  GtUint64 *valueptr;
  GtWord corrected_key = (GtWord) key + dist->correction_base;
  gt_assert(dist->li_ull_hashmap);
  if (!(valueptr = li_ull_gt_hashmap_get(dist->li_ull_hashmap,
                                         corrected_key)))
    return 0;
  return *valueptr;
}

static GtUint64 encdesc_distri_get_symbol_freq(const void *distri,
                                                         GtUword symbol)
{
  GtDiscDistri *distr = (GtDiscDistri*) distri;
  return gt_disc_distri_get(distr, symbol);
}

GtUint64 encdesc_hashmap_distr_get(const void *hm_distri,
                                             GtUword key)
{
  GtHashtable *hashmap = (GtHashtable*) hm_distri;
  GtUint64 *valueptr;
  gt_assert(hashmap != NULL);
  if (!(valueptr = li_ull_gt_hashmap_get(hashmap, (GtWord) key)))
    return 0;
  return *valueptr;
}

static void encdesc_init_huffman(GtEncdesc *encdesc)
{
  DescField *field;
  GtUword field_idx,
                alphabet_size = 0,
                char_idx;
  GtUword const char_alphabet_size = 256UL;
  EncdescHuffDist huffdist;

  for (field_idx = 0;
       field_idx < encdesc->num_of_fields;
       field_idx++) {
    field = &encdesc->fields[field_idx];
    if (!field->is_cons) {
      if (field->is_numeric) {
        if (field->use_delta_coding && field->use_hc) {
          huffdist.correction_base = field->min_delta;
          gt_assert(field->delta_values != NULL);
          huffdist.li_ull_hashmap = field->delta_values;
          alphabet_size =
            (GtUword) labs(field->max_delta - field->min_delta) + 1;
        }
        else if (!field->is_value_cons && field->use_hc) {
          huffdist.correction_base = field->min_value;
          gt_assert(field->num_values != NULL);
          huffdist.li_ull_hashmap = field->num_values;
          alphabet_size =
            (GtUword) labs(field->max_value - field->min_value) + 1;
        }

        if (field->use_hc) {
          field->huffman_num =
            gt_huffman_new(&huffdist, encdesc_hashmap_distr_get_corrected,
                           alphabet_size);
        }

        if (field->has_zero_padding && !field->fieldlen_is_const) {
          field->huffman_zero_count =
            gt_huffman_new(field->zero_count,
                           encdesc_distri_get_symbol_freq,
                           (GtUword) field->max_zero + 1);
        }
      }
      else {
        field->huffman_chars = gt_calloc((size_t) (field->max_len + 1),
                                         sizeof (field->huffman_chars));
        for (char_idx = 0; char_idx < field->max_len; char_idx++) {
          if (char_idx >= field->len ||
              !gt_bittab_bit_is_set(field->bittab, char_idx)) {
            field->huffman_chars[char_idx] =
              gt_huffman_new(field->chars[char_idx],
                             encdesc_hashmap_distr_get,
                             char_alphabet_size);
          }
        }
      }
    }
  }
}

int gt_encdesc_encoder_encode(GtEncdescEncoder *ee,
                              GtCstrIterator *cstr_iterator,
                              const char *name, GtError *err)
{
  int had_err = 0;
  bool is_not_at_pageborder;
  FILE *fp = NULL;
  GtStr *name1;
  GtWord pos = 0,
       start_of_encoding;
  GtUword dummy = 0,
                pagesize;

  gt_assert(ee != NULL);
  gt_assert(cstr_iterator != NULL);
  gt_assert(name != NULL);
  gt_error_check(err);
  if (ee->timer != NULL) {
    gt_timer_show_progress(ee->timer, "analyze descriptions", stdout);
  }
  had_err = encdesc_analyze_descs(ee->encdesc, cstr_iterator, err);

  if (!had_err) {
    if (ee->timer != NULL) {
      gt_timer_show_progress(ee->timer, "write encoding header", stdout);
    }

    fp = gt_fa_fopen_with_suffix(name, GT_ENCDESC_FILESUFFIX, "wb", err);
    if (fp == NULL)
      had_err = -1;
  }

  if (!had_err)
    had_err = encdesc_write_header(ee->encdesc, fp, err);
  if (!had_err) {
    if (ee->timer != NULL) {
      gt_timer_show_progress(ee->timer, "calculate huffmans", stdout);
    }
    encdesc_init_huffman(ee->encdesc);
    if (ee->timer != NULL) {
      gt_timer_show_progress(ee->timer, "write encoding", stdout);
    }
    gt_error_check(err);
    pos = ftell(fp);
    gt_xfwrite_one(&dummy, fp);

    pagesize = ee->encdesc->pagesize;
    is_not_at_pageborder = (ftell(fp) % pagesize) != 0;
    if (is_not_at_pageborder)
      start_of_encoding = (ftell(fp) / pagesize + 1) * pagesize;
    else
      start_of_encoding = ftell(fp);

    gt_xfseek(fp, start_of_encoding, SEEK_SET);

    if (ee->page_sampling)
      ee->encdesc->sampling = gt_sampling_new_page(ee->sampling_rate,
                                                   (off_t) start_of_encoding);
    else if (ee->regular_sampling)
      ee->encdesc->sampling = gt_sampling_new_regular(ee->sampling_rate,
                                                     (off_t) start_of_encoding);
    had_err = encdesc_write_encoding(ee->encdesc, cstr_iterator, fp, err);
  }
  if (!had_err) {
    gt_xfseek(fp, pos, SEEK_SET);
    if (ee->encdesc->sampling != NULL)
      gt_xfwrite(&ee->encdesc->start_of_samplingtab,
                 sizeof (ee->encdesc->start_of_samplingtab), (size_t) 1, fp);
    else {
      const GtWord null = 0;
      gt_xfwrite_one(&null, fp);
    }
  }

  gt_fa_xfclose(fp);
  if (!had_err) {
    if (ee->timer != NULL) {
      gt_timer_show_progress(ee->timer,
                             "description encoding finished", stdout);
    }
    if (gt_log_enabled()) {
      GtUword rate;
      name1 = gt_str_new_cstr(name);
      gt_str_append_cstr(name1, GT_ENCDESC_FILESUFFIX);
      gt_log_log("description encoding overview:");
      gt_log_log("==>");
      if (ee->encdesc->sampling != NULL) {
        rate = gt_sampling_get_rate(ee->encdesc->sampling);
        if (gt_sampling_is_regular(ee->encdesc->sampling)) {
          gt_log_log("applied sampling technique:"
                     " sampling every "GT_WU"th description",
                     rate);
        }
        else {
          gt_log_log("applied sampling technique:"
                     " sampling every "GT_WU"th page",
                     rate);
        }
      }
      else {
        gt_log_log("applied sampling technique: none");
      }

      gt_log_log("total number of encoded descriptions: "GT_WU"",
                 ee->encdesc->num_of_descs);
      gt_log_log("total number of encoded characters: "GT_LLU"",
                 ee->encdesc->total_num_of_chars);
      gt_log_log("bits per character encoding: %f",
                 (gt_file_estimate_size(gt_str_get(name1)) * 8.0) /
                 ee->encdesc->total_num_of_chars);
      gt_log_log("<==");
      gt_str_delete(name1);
    }
  }
  return had_err;
}

static void encdesc_read_samplingtab(GtEncdesc *encdesc,
                                     FILE *fp)
{
  bool is_not_at_pageborder;
  GT_UNUSED size_t read;

  read = gt_xfread_one(&encdesc->start_of_encoding, fp);
  gt_assert(read == (size_t) 1);

  is_not_at_pageborder = (ftell(fp) % encdesc->pagesize) != 0;
  if (is_not_at_pageborder)
    encdesc->start_of_encoding = (ftell(fp) / encdesc->pagesize + 1) *
                                 encdesc->pagesize;
  else
    encdesc->start_of_encoding = ftell(fp);

  if (encdesc->start_of_samplingtab != 0) {
    gt_xfseek(fp, encdesc->start_of_samplingtab, SEEK_SET);
    encdesc->sampling = gt_sampling_read(fp);
  }
}

GtEncdesc* gt_encdesc_load(const char *name,
                           GtError *err)
{
  int had_err = 0;
  GtEncdesc *encdesc = NULL;
  FILE *fp;
  GtStr *filename;
  int fd;
  GtUword const pages_to_map = 5UL;

  gt_assert(name);
  encdesc = encdesc_new();

  filename = gt_str_new_cstr(name);
  gt_str_append_cstr(filename, GT_ENCDESC_FILESUFFIX);
  fp = gt_fa_fopen_with_suffix(name, GT_ENCDESC_FILESUFFIX, "rb", err);
  if (fp == NULL) {
    gt_assert(gt_error_is_set(err));
    had_err = 1;
  }

  if (!had_err) {
    fd = open(gt_str_get(filename), O_RDONLY);
    if (fd == -1) {
      gt_error_set(err, "open(): cannot read file %s, error: %s",
                   gt_str_get(filename), strerror(errno));
      had_err = 1;
    }
  }
  if (!had_err)
    had_err = encdesc_read_header(encdesc, fp, err);

  if (!had_err) {
    encdesc_init_huffman(encdesc);

    encdesc_read_samplingtab(encdesc, fp);
    gt_fa_fclose(fp);

    encdesc->bitinstream =
      gt_bitinstream_new(gt_str_get(filename),
                         (size_t) encdesc->start_of_encoding, pages_to_map);
    gt_str_delete(filename);

  }
  else {
    gt_encdesc_delete(encdesc);
    encdesc = NULL;
  }
  return encdesc;
}

static int encdesc_next_desc(GtEncdesc *encdesc, GtStr *desc, GtError *err)
{
  int stat, had_err = 0;
  const int success = 1;
  bool bit,
       sampled = false;
  unsigned readbits,
           bits_to_read;
  GtWord tmp = 0;
  GtUword cur_field_num,
                fieldlen = 0,
                idx,
                numoffields,
                nearestsample,
                zero_count = 0,
                tmp_symbol = 0;
  size_t startofnearestsample;
  GtBitsequence bitseq;
  GtHuffmanBitwiseDecoder *huff_bitwise_decoder;

  gt_log_log("current_desc: "GT_WU"", encdesc->cur_desc);
  if (encdesc->cur_desc == encdesc->num_of_descs) {
    gt_error_set(err,"nothing done, eof?");
    return had_err;
  }

  if (encdesc->sampling != NULL &&
      encdesc->cur_desc == gt_sampling_get_next_elementnum(encdesc->sampling)) {
    int sample_status;
    sample_status = gt_sampling_get_next_sample(encdesc->sampling,
                                    &nearestsample,
                                    &startofnearestsample);
    if (sample_status == 1) {
      gt_bitinstream_reinit(encdesc->bitinstream,
                            startofnearestsample);
      sampled = true;
    }
    else
      had_err = sample_status;
  }

  if (had_err)
    gt_error_set(err, "sampling did not work");

  if (desc != NULL)
    gt_str_reset(desc);
  if (!had_err && !encdesc->num_of_fields_is_cons) {
    bits_to_read = encdesc->bits_per_field;
    bitseq = 0;
    for (readbits = 0; !had_err && readbits < bits_to_read; readbits++) {
      bitseq = bitseq << 1;
      if (gt_bitinstream_get_next_bit(encdesc->bitinstream, &bit) != 1) {
        gt_error_set(err, "could not get next bit");
        had_err = -1;
      }
      else if (bit)
        bitseq = bitseq | (GtBitsequence) 1;
    }
    numoffields = (GtUword) bitseq;
  }
  else
    numoffields = encdesc->num_of_fields;

  for (cur_field_num = 0;
       !had_err && cur_field_num < numoffields;
       cur_field_num++) {
    DescField *cur_field = &encdesc->fields[cur_field_num];
    if (cur_field->is_cons) {
      if (desc != NULL) {
        gt_str_append_cstr(desc, cur_field->data);
        gt_str_append_char(desc, cur_field->sep);
      }
      continue;
    }
    if (cur_field->is_numeric) {
      if (cur_field->has_zero_padding && !cur_field->fieldlen_is_const) {
        huff_bitwise_decoder = gt_huffman_bitwise_decoder_new(
                                           cur_field->huffman_zero_count, err);
        stat = -1;
        while (!had_err && stat != 0) {
          if (gt_bitinstream_get_next_bit(encdesc->bitinstream, &bit) != 1) {
            gt_error_set(err, "could not get next bit");
            had_err = -1;
          }
          else {
            stat = gt_huffman_bitwise_decoder_next(huff_bitwise_decoder, bit,
                                                   &zero_count, err);
            if (stat == -1) {
              had_err = stat;
              gt_assert(gt_error_is_set(err));
            }
          }
        }
        gt_huffman_bitwise_decoder_delete(huff_bitwise_decoder);
        for (idx = 0;
             !had_err && desc != NULL && idx < zero_count;
             idx++)
          gt_str_append_char(desc, '0');
      }
      /* read absolute value if description is first or sampled */
      if (!had_err && (encdesc->cur_desc == 0 || sampled)) {
        bits_to_read = cur_field->bits_per_value;
        bitseq = 0;
        for (readbits = 0;
             !had_err && readbits < bits_to_read;
             readbits++) {
          bitseq = bitseq << 1;
          if (gt_bitinstream_get_next_bit(encdesc->bitinstream, &bit) != 1) {
            gt_error_set(err, "could not get next bit");
            had_err = -1;
          }
          else if (bit)
            bitseq = bitseq | (GtBitsequence) 1;

        }
        if (!had_err) {
          tmp = (GtWord) bitseq + cur_field->min_value;
          cur_field->prev_value = tmp;
        }
      }
      else if (!had_err) {
        if (!cur_field->is_value_cons || !cur_field->is_delta_cons) {
          if (cur_field->bits_per_num) {
            if (cur_field->use_hc) {
              huff_bitwise_decoder = gt_huffman_bitwise_decoder_new(
                                                  cur_field->huffman_num, err);
              stat = 1;
              while (!had_err && stat != 0) {
                if (gt_bitinstream_get_next_bit(encdesc->bitinstream,
                                                &bit) != 1) {
                  gt_error_set(err, "could not get next bit");
                  had_err = -1;
                }
                else {
                  stat = gt_huffman_bitwise_decoder_next(huff_bitwise_decoder,
                                                         bit, &tmp_symbol, err);
                  if (stat == -1) {
                    had_err = stat;
                    gt_assert(gt_error_is_set(err));
                  }
                  else
                    tmp = (GtWord) tmp_symbol;
                }
              }
              gt_huffman_bitwise_decoder_delete(huff_bitwise_decoder);
            }
            else {
              bits_to_read = cur_field->bits_per_num;
              bitseq = 0;
              for (readbits = 0;
                   !had_err && readbits < bits_to_read;
                   readbits++) {
                bitseq = bitseq << 1;
                if (gt_bitinstream_get_next_bit(encdesc->bitinstream,
                                                &bit) != 1) {
                  gt_error_set(err, "could not get next bit");
                  had_err = -1;
                }
                else if (bit)
                  bitseq = bitseq | (GtBitsequence) 1;
              }
              if (!had_err)
                tmp = (GtWord) bitseq;
            }
          }
          else
            tmp = 0;
        }
        else {
          if (cur_field->use_delta_coding)
            tmp = 0;
          else
            tmp = cur_field->prev_value - cur_field->min_value;
        }
        if (cur_field->use_delta_coding)
          tmp += cur_field->prev_value + cur_field->min_delta;
        else
          tmp += cur_field->min_value;
      }
      if (!had_err) {
        cur_field->prev_value = tmp;
        if (cur_field->has_zero_padding && cur_field->fieldlen_is_const) {
          zero_count = cur_field->len -
                       encdesc_digits_per_value((GtUword) tmp, 10UL);
          for (idx = 0;
               desc != NULL && idx < zero_count;
               idx++)
            gt_str_append_char(desc, '0');
        }
        if (desc != NULL) {
          gt_str_append_ulong(desc, (GtUword) tmp);
          gt_str_append_char(desc, cur_field->sep);
        }
        continue;
      }
    }
    /* variable cur_field len */
    if (!cur_field->fieldlen_is_const) {
      bits_to_read = cur_field->bits_per_len;
      bitseq = 0;
      for (readbits = 0;
           !had_err && readbits < bits_to_read;
           readbits++) {
        bitseq = bitseq << 1;
        if (gt_bitinstream_get_next_bit(encdesc->bitinstream, &bit) != 1) {
          gt_error_set(err, "could not get next bit");
          had_err =  -1;
        }
        else if (bit)
          bitseq = bitseq | (GtBitsequence) 1;

      }
      if (!had_err)
        fieldlen = (GtUword) bitseq + cur_field->min_len;
    }
    else
      fieldlen = cur_field->len;

    for (idx = 0; !had_err && idx < fieldlen; idx++) {
      if (idx < cur_field->len &&
          gt_bittab_bit_is_set(cur_field->bittab, idx)) {
        if (desc != NULL)
          gt_str_append_char(desc, cur_field->data[idx]);
      }
      else {
        huff_bitwise_decoder = gt_huffman_bitwise_decoder_new(
                                           cur_field->huffman_chars[idx], err);
        stat = -1;
        while (!had_err && stat != 0) {
          if (gt_bitinstream_get_next_bit(encdesc->bitinstream, &bit) != 1) {
            gt_error_set(err, "could not get next bit");
            had_err = -1;
          }
          else {
            stat = gt_huffman_bitwise_decoder_next(huff_bitwise_decoder,
                                                   bit, &tmp_symbol,err);
            if (stat == -1) {
              gt_assert(gt_error_is_set(err));
              had_err = stat;
            }
            else
              tmp = (GtWord) tmp_symbol;
          }
        }
        if (!had_err && desc != NULL) {
          gt_assert(tmp < 256L);
          gt_str_append_char(desc, (char) tmp);
        }
        gt_huffman_bitwise_decoder_delete(huff_bitwise_decoder);
      }
    }
    if (!had_err && desc != NULL)
      gt_str_append_char(desc, cur_field->sep);
  }
  if (desc != NULL)
    gt_str_set_length(desc, gt_str_length(desc) - 1);
  if (!had_err) {
    encdesc->cur_desc++;
  }
  if (had_err)
    return had_err;
  else
    return success;
}

GtUword gt_encdesc_num_of_descriptions(GtEncdesc *encdesc)
{
  gt_assert(encdesc);
  return encdesc->num_of_descs;
}

int gt_encdesc_decode(GtEncdesc *encdesc,
                      GtUword num,
                      GtStr *desc,
                      GtError *err)
{
  int had_err = 0;
  GtUword descs2read = 0,
                nearestsample = 0,
                idx;
  size_t startofnearestsample = 0;

  gt_assert(encdesc);
  gt_assert(desc);
  gt_assert(num < encdesc->num_of_descs);

  if (encdesc->cur_desc == num) {
    return encdesc_next_desc(encdesc, desc, err);
  }

  if (encdesc->sampling != NULL) {
    (void) gt_sampling_get_page(encdesc->sampling,
                                num,
                                &nearestsample,
                                &startofnearestsample);
    /* nearestsample <= cur_read < readnum: current sample is the right one */
    if (nearestsample <= encdesc->cur_desc && encdesc->cur_desc <= num)
      descs2read = num - encdesc->cur_desc;
    else { /* reset decoder to new sample */
      gt_bitinstream_reinit(encdesc->bitinstream,
                            startofnearestsample);
      encdesc->cur_desc = nearestsample;
      descs2read = num - nearestsample;
    }
  }
  else {
    if (encdesc->cur_desc <= num)
      descs2read = num - encdesc->cur_desc;
    else {
      gt_bitinstream_reinit(encdesc->bitinstream,
                            (size_t) encdesc->start_of_encoding);
      descs2read = num;
      encdesc->cur_desc = 0;
    }
  }

  /* decode all description until the requested */
  for (idx = 0; !had_err && idx < descs2read; idx++) {
    if (encdesc_next_desc(encdesc, NULL, err) != 1)
      had_err = -1;
  }

  /* decode the requested description */
  if (!had_err) {
    if (encdesc_next_desc(encdesc, desc, err) != 1)
      had_err = -1;
  }
  else
    gt_assert(gt_error_is_set(err));
  return had_err;
}

static void encdesc_delete_desc_fields(DescField *fields,
                                      GtUword numoffields)
{
  GtUword idx,
                j_idx;

  if (!fields)
    return;
  for (idx = 0; idx < numoffields; idx++) {
    gt_huffman_delete(fields[idx].huffman_num);
    gt_huffman_delete(fields[idx].huffman_zero_count);
    gt_free(fields[idx].data);
    if (fields[idx].chars != NULL) {
      for (j_idx = 0; j_idx < fields[idx].max_len; j_idx++)
        gt_hashtable_delete(fields[idx].chars[j_idx]);
      gt_free(fields[idx].chars);
    }
    gt_hashtable_delete(fields[idx].delta_values);
    gt_hashtable_delete(fields[idx].num_values);
    gt_disc_distri_delete(fields[idx].zero_count);

    if (fields[idx].huffman_chars != NULL) {
      for (j_idx = 0; j_idx < fields[idx].max_len; j_idx++) {
        if (j_idx >= fields[idx].len ||
            !gt_bittab_bit_is_set(fields[idx].bittab, j_idx))
          gt_huffman_delete(fields[idx].huffman_chars[j_idx]);
      }
      gt_free(fields[idx].huffman_chars);
    }
    gt_bittab_delete(fields[idx].bittab);
  }
  gt_free(fields);
}

void gt_encdesc_delete(GtEncdesc *encdesc)
{
  if (!encdesc) return;
  gt_bitinstream_delete(encdesc->bitinstream);
  GT_FREEARRAY(&encdesc->num_of_fields_tab, GtUlong);
  encdesc_delete_desc_fields(encdesc->fields, encdesc->num_of_fields);
  gt_sampling_delete(encdesc->sampling);
  gt_free(encdesc);
}

void gt_encdesc_encoder_delete(GtEncdescEncoder *ee)
{
  if (ee != NULL) {
    gt_encdesc_delete(ee->encdesc);
    gt_free(ee);
  }
}

int gt_encdesc_unit_test(GtError *err)
{
  int had_err = 0;
  EncdescWriteInfo *info = gt_calloc((size_t) 1, sizeof (*info));
  info->codes = gt_malloc(sizeof (*info->codes));

  GT_INITARRAY(info->codes, EncdescCode);

  /* test count_leading_zeros */
  if (!had_err) {
    char *test = "000156";
    gt_ensure(count_leading_zeros(test) == 3U);

    test = "x";
    gt_ensure(count_leading_zeros(test) == 0);

    test = "0000";
    gt_ensure(count_leading_zeros(test) == 3U);

    test = "";
    gt_ensure(count_leading_zeros(test) == 0U);
  }

  /* test parse_number_out_of_current_field */
  if (!had_err) {
    GtWord retval = 0;
    info->descbuffer = "abc00666";
    info->cur_field_start_pos = 3UL;
    parse_number_out_of_current_field(info, &retval);
    gt_ensure(retval == 666L);

    info->descbuffer = "abc\000666\000abc";
    info->cur_field_start_pos = 4UL;
    parse_number_out_of_current_field(info, &retval);
    gt_ensure(retval == 666L);

    info->descbuffer = "aaa-123\000aaa";
    info->cur_field_start_pos = 3UL;
    parse_number_out_of_current_field(info, &retval);
    gt_ensure(retval == -123L);
  }

  GT_FREEARRAY(info->codes, EncdescCode);
  gt_free(info->codes);
  gt_free(info);
  return had_err;
}
