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

#include <math.h>

#include "core/ma_api.h"
#include "core/xansi_api.h"
#include "core/log_api.h"
#include "extended/encdesc_header_io.h"
#include "extended/encdesc_rep.h"

static inline void encdesc_gt_xfwrite(void *ptr,
                                      size_t size,
                                      size_t nmemb,
                                      FILE *stream)
{
  gt_xfwrite((const void*) ptr, size, nmemb, stream);
}

static inline void encdesc_gt_xfread(void *ptr,
                                      size_t size,
                                      size_t nmemb,
                                      FILE *stream)
{
  GT_UNUSED size_t written;
  written = gt_xfread(ptr, size, nmemb, stream);
  gt_assert(written == nmemb);
}

typedef void (*EncdescIOFunc)(void *ptr,
                              size_t size,
                              size_t nmemb,
                              FILE *stream);

static inline void encdesc_header_io_basics(GtEncdesc *encdesc,
                                            FILE *fp,
                                            EncdescIOFunc io_func)
{
  io_func(&(encdesc->num_of_descs), sizeof (encdesc->num_of_descs),
          (size_t) 1, fp);
  io_func(&(encdesc->num_of_fields), sizeof (encdesc->num_of_fields),
          (size_t) 1, fp);
  io_func(&(encdesc->num_of_fields_is_cons), sizeof (bool),
          (size_t) 1, fp);
  io_func(&(encdesc->bits_per_field), sizeof (encdesc->bits_per_field),
          (size_t) 1, fp);
}

static void io_field_sep_and_is_cons(DescField *field,
                           FILE *fp,
                           EncdescIOFunc io_func)
{
  io_func(&field->sep, sizeof (char), (size_t) 1, fp);
  io_func(&field->is_cons, sizeof (bool), (size_t) 1, fp);
}

static void write_cons_field_header(DescField *field,
                             FILE *fp)
{
    gt_xfwrite(&field->len, sizeof (unsigned int), (size_t) 1, fp);

    gt_xfwrite(field->data, sizeof (char), (size_t) field->len, fp);
}

static void read_cons_field_header(DescField *field,
                            FILE *fp)
{
  (void) gt_xfread(&field->len, sizeof (unsigned int), (size_t) 1, fp);
  field->data = gt_calloc((size_t) (field->len + 1), sizeof (char));
  (void) gt_xfread(field->data, sizeof (char), (size_t) field->len, fp);
}

static bool io_zero_padding_needs_dist(DescField *field,
                                       FILE *fp,
                                       EncdescIOFunc io_func)
{
  io_func(&field->has_zero_padding, sizeof (bool), (size_t) 1, fp);
  if (field->has_zero_padding) {
    io_func(&field->fieldlen_is_const, sizeof (bool), (size_t) 1, fp);
    if (field->fieldlen_is_const) {
      io_func(&field->len, sizeof (field->len), (size_t) 1, fp);
    }
    else {
      io_func(&field->max_zero, sizeof (field->max_zero), (size_t) 1, fp);
      return true;
    }
  }
  return false;
}

typedef struct data {
  FILE *fp;
  unsigned long written_elems;
  long minimum_element;
} EncsecDistriData;

static void encdesc_distri_iter_count(GT_UNUSED unsigned long symbol,
                                      GT_UNUSED unsigned long long freq,
                                      void *data)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
}

static void encdesc_distri_iter_write(unsigned long symbol,
                                      unsigned long long freq,
                                      void *data)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
  gt_xfwrite(&symbol, sizeof (symbol), (size_t) 1, this_data->fp);
  gt_xfwrite(&freq, sizeof (freq), (size_t) 1, this_data->fp);

}

static void io_min_max_calc_ranges(DescField *field,
                       FILE *fp,
                       EncdescIOFunc io_func)
{
  io_func(&field->min_value, sizeof (field->min_value), (size_t) 1, fp);
  io_func(&field->max_value, sizeof (field->min_value), (size_t) 1, fp);
  gt_assert(field->min_value <= field->max_value);
  io_func(&field->min_delta, sizeof (field->min_delta), (size_t) 1, fp);
  io_func(&field->max_delta, sizeof (field->max_delta), (size_t) 1, fp);
  gt_assert(field->min_delta <= field->max_delta);

  io_func(&field->is_value_cons, sizeof (bool), (size_t) 1, fp);
  if (field->is_value_cons)
  io_func(&field->global_value, sizeof (field->global_value), (size_t) 1, fp);
  if (field->is_value_cons) {
    gt_assert(field->min_value == field->global_value);
    gt_assert(field->max_value == field->global_value);
  }
  io_func(&field->is_delta_cons, sizeof (bool), (size_t) 1, fp);
  if (field->is_delta_cons)
  io_func(&field->global_delta, sizeof (field->global_delta), (size_t) 1, fp);
  if (field->is_delta_cons) {
    gt_assert(field->min_delta == field->global_delta);
    gt_assert(field->max_delta == field->global_delta);
  }

  io_func(&field->use_delta_coding, sizeof (field->use_delta_coding),
          (size_t) 1, fp);
  io_func(&field->delta_values_size, sizeof (field->delta_values_size),
          (size_t) 1, fp);
  io_func(&field->num_values_size, sizeof (field->num_values_size),
          (size_t) 1, fp);
  io_func(&field->bits_per_num, sizeof (field->bits_per_num), (size_t) 1, fp);
  io_func(&field->bits_per_value, sizeof (field->bits_per_value),
          (size_t) 1, fp);
}

static void numeric_field_check_distri_dependence(DescField *field,
                                                  bool *needs_delta_dist,
                                                  bool *needs_value_dist)
{
  field->use_hc = false;
  if (field->use_delta_coding) {
    if (!field->is_delta_cons &&
        field->delta_values_size <= GT_ENCDESC_MAX_NUM_VAL_HUF) {
      *needs_delta_dist = true;
      field->use_hc = true;
      gt_log_log("delta_values_size: %lu", field->delta_values_size);
    }
  }
  else {
    if (!field->is_value_cons && field->num_values_size > 0 &&
        field->num_values_size <= GT_ENCDESC_MAX_NUM_VAL_HUF) {
      *needs_value_dist = true;
      field->use_hc = true;
      gt_log_log("num_values_size: %lu", field->num_values_size);
    }
  }
}

static enum iterator_op encdesc_li_ull_hashmap_iter_count(
                                             GT_UNUSED long key,
                                             GT_UNUSED unsigned long long value,
                                             void *data,
                                             GT_UNUSED GtError *err)
{
  unsigned long *count = (unsigned long*) data;
  (void) (*count)++;
  /* always continue since we cannot fail */
  return CONTINUE_ITERATION;
}

static unsigned long get_hashmap_distri_size(GtHashtable *h_table)
{
  unsigned long count = 0;
  /* error is NULL as encdesc_li_ull_hashmap_iter_count() is sane */
  (void) li_ull_gt_hashmap_foreach(h_table,
                                   encdesc_li_ull_hashmap_iter_count,
                                   &count,
                                   NULL);
  return count;
}

static enum iterator_op encdesc_li_ull_hashmap_iter_write(
                                                       long key,
                                                       unsigned long long value,
                                                       void *data,
                                                       GT_UNUSED GtError *err)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
  gt_xfwrite(&key, sizeof (key), (size_t) 1, this_data->fp);
  gt_xfwrite(&value, sizeof (value), (size_t) 1, this_data->fp);
  /* always continue since we use gt_xfwrite() anyway */
  return CONTINUE_ITERATION;
}

static void write_hashmap_distri(EncsecDistriData *data,
                                 GtHashtable *h_table,
                                 GT_UNUSED unsigned long size)
{
  data->written_elems = 0;
  /* error is NULL as encdesc_li_ull_hashmap_iter_write() is sane */
  (void) li_ull_gt_hashmap_foreach(h_table,
                                   encdesc_li_ull_hashmap_iter_write,
                                   data,
                                   NULL);
  if (data->written_elems != size)
    gt_log_log("%lu != %lu", size, data->written_elems);
  gt_assert(data->written_elems == size);
}

static void read_hashmap_distri(unsigned long size,
                                GtHashtable *h_table,
                                FILE *fp)
{
  unsigned long idx;
  long symbol;
  unsigned long long freq;

  gt_assert(h_table != NULL);
  for (idx = 0; idx < size; idx++) {
    (void) gt_xfread(&symbol, sizeof (symbol), (size_t) 1, fp);
    (void) gt_xfread(&freq, sizeof (freq), (size_t) 1, fp);
    gt_assert(li_ull_gt_hashmap_get(h_table, symbol) == 0);
    (void) li_ull_gt_hashmap_add(h_table, symbol, freq);
  }
}

static void read_numeric_field_header(DescField *field,
                                      FILE *fp)
{
  bool needs_zero_dist = false,
       needs_delta_dist = false,
       needs_value_dist = false;
  unsigned long num_of_zero_leaves;
  /* EncdescHuffDist huffdist; */

  needs_zero_dist = io_zero_padding_needs_dist(field,
                                               fp,
                                               encdesc_gt_xfread);
  io_min_max_calc_ranges(field,
                         fp,
                         encdesc_gt_xfread);

  numeric_field_check_distri_dependence(field,
                                        &needs_delta_dist,
                                        &needs_value_dist);

  if (needs_delta_dist) {
    field->delta_values = li_ull_gt_hashmap_new();
    read_hashmap_distri(field->delta_values_size,
                        field->delta_values,
                        fp);
    /* huffdist.li_ull_hashmap = field->delta_values;
    huffdist.correction_base = field->min_delta;
    field->huffman_num = gt_huffman_new(&huffdist,
                                        encdesc_hashmap_distr_get_corrected,
                                        labs(field->max_delta -
                                          field->min_delta) + 1); */
  }

  if (needs_value_dist) {
    field->num_values = li_ull_gt_hashmap_new();
    read_hashmap_distri(field->num_values_size,
                        field->num_values,
                        fp);
    /* huffdist.li_ull_hashmap = field->num_values;
    huffdist.correction_base = field->min_value;
    field->huffman_num = gt_huffman_new(huffdist,
                                        encdesc_hashmap_distr_get_corrected,
                                        labs(field->max_value -
                                          field->min_value) + 1); */
  }

  if (needs_zero_dist) {
    unsigned long idx,
                  symbol;
    unsigned long long freq;
    field->zero_count = gt_disc_distri_new();
    (void) gt_xfread(&num_of_zero_leaves, sizeof (num_of_zero_leaves),
                     (size_t) 1, fp);
    for (idx = 0; idx < num_of_zero_leaves; idx++) {
      (void) gt_xfread(&symbol, sizeof (symbol), (size_t) 1, fp);
      (void) gt_xfread(&freq, sizeof (freq), (size_t) 1, fp);
      gt_disc_distri_add_multi(field->zero_count, symbol, freq);
    }
    /* field->huffman_zero_count =
      gt_huffman_new(field->zero_count,
                     encdesc_distri_get_symbol_freq,
                     &field->max_zero + 1); */
  }
}

static void write_numeric_field_header(DescField *field,
                                       FILE *fp)
{
  bool needs_zero_dist = false,
       needs_delta_dist = false,
       needs_value_dist = false;
  unsigned long num_of_zero_leaves;
  EncsecDistriData data;

  data.fp = fp;

  needs_zero_dist = io_zero_padding_needs_dist(field,
                                               fp,
                                               encdesc_gt_xfwrite);
  io_min_max_calc_ranges(field,
                         fp,
                         encdesc_gt_xfwrite);

  numeric_field_check_distri_dependence(field,
                                        &needs_delta_dist,
                                        &needs_value_dist);

  if (needs_delta_dist) {
    write_hashmap_distri(&data,
                         field->delta_values,
                         field->delta_values_size);
  }

  if (needs_value_dist) {
    write_hashmap_distri(&data,
                         field->num_values,
                         field->num_values_size);
  }

  if (needs_zero_dist) {
    /* TODO change this distribution, to hashtable like the others, this reduces
     * functions */
    data.written_elems = 0;
    gt_disc_distri_foreach(field->zero_count,
                           encdesc_distri_iter_count,
                           &data);
    num_of_zero_leaves = data.written_elems;
    gt_xfwrite(&num_of_zero_leaves, sizeof (num_of_zero_leaves),
               (size_t) 1, fp);
    data.written_elems = 0;
    gt_disc_distri_foreach(field->zero_count,
                           encdesc_distri_iter_write,
                           &data);
    gt_assert(data.written_elems == num_of_zero_leaves);
  }
}

static void io_field_len_header(DescField *field,
                                FILE *fp,
                                EncdescIOFunc io_func)
{
  io_func(&field->fieldlen_is_const, sizeof (bool), (size_t) 1, fp);
  io_func(&field->len, sizeof (field->len), (size_t) 1, fp);
  io_func(&field->max_len, sizeof (field->max_len), (size_t) 1, fp);
  io_func(&field->min_len, sizeof (field->min_len), (size_t) 1, fp);
  io_func(&field->bits_per_len, sizeof (field->bits_per_len), (size_t) 1, fp);
}

static void read_field_header_bittab(DescField *field,
                                     FILE *fp)
{
  unsigned long char_idx,
                num_of_chars = field->len / sizeof (char);
  char cc;
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; char_idx < num_of_chars; char_idx++) {
    (void) gt_xfread(&cc, sizeof (cc), (size_t) 1, fp);
    for (bit_idx = 0; bit_idx < sizeof (char); bit_idx++) {
      if (cc & (1 << bit_idx))
        gt_bittab_set_bit(field->bittab,
                          (unsigned long) ((sizeof (char) * char_idx) +
                                           bit_idx));
    }
  }
}

static void write_field_header_bittab(DescField *field, FILE *fp)
{
  unsigned long char_idx,
                num_of_chars = field->len / sizeof (char);
  char cc = 0;
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; char_idx < num_of_chars; char_idx++) {
    for (bit_idx = 0; bit_idx < sizeof (char); bit_idx++) {
      if (gt_bittab_bit_is_set(field->bittab,
                               (unsigned long) ((sizeof (char) * char_idx) +
                                                bit_idx)))
        cc |= 1 << bit_idx;
    }
    gt_xfwrite(&cc, sizeof (cc), (size_t) 1, fp);
    cc = 0;
  }
}

/* TODO combine field_char_dist */
static void write_field_char_dists(DescField *field,
                                   FILE *fp)
{
  unsigned long char_idx,
                distr_len;
  EncsecDistriData data;
  data.fp = fp;
  data.minimum_element = 0;
  data.written_elems = 0;

  for (char_idx = 0; char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      distr_len = get_hashmap_distri_size(field->chars[char_idx]);
      gt_xfwrite(&distr_len, sizeof (distr_len), (size_t) 1, fp);
      write_hashmap_distri(&data,
                           field->chars[char_idx],
                           distr_len);
    }
  }
}

static void read_field_char_dists(DescField *field,
                                  FILE *fp)
{
  unsigned long char_idx,
                distr_len;

  for (char_idx = 0; char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      (void) gt_xfread(&distr_len, sizeof (distr_len), (size_t) 1, fp);
      field->chars[char_idx] = li_ull_gt_hashmap_new();
      read_hashmap_distri(distr_len,
                          field->chars[char_idx],
                          fp);
    }
  }
}

void encdesc_write_header(GtEncdesc *encdesc, FILE *fp)
{
  unsigned long cur_field_num;
  DescField *cur_field;

  encdesc_header_io_basics(encdesc, fp, encdesc_gt_xfwrite);

  for (cur_field_num = 0;
       cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    io_field_sep_and_is_cons(cur_field, fp, encdesc_gt_xfwrite);

    if (cur_field->is_cons) {
      write_cons_field_header(cur_field, fp);
    }
    else {
      gt_xfwrite(&cur_field->is_numeric, sizeof (bool), (size_t) 1, fp);
      if (cur_field->is_numeric) {
        write_numeric_field_header(cur_field, fp);
      }
      else {

        io_field_len_header(cur_field, fp, encdesc_gt_xfwrite);

        gt_xfwrite(cur_field->data, sizeof (char), (size_t) cur_field->len, fp);

        write_field_header_bittab(cur_field, fp);

        write_field_char_dists(cur_field, fp);
      }
    }
  }
}

void encdesc_read_header(GtEncdesc *encdesc, FILE *fp)
{
  unsigned long cur_field_num;
  DescField *cur_field;

  encdesc_header_io_basics(encdesc, fp, encdesc_gt_xfread);

  encdesc->fields = gt_calloc((size_t) encdesc->num_of_fields,
                              sizeof (DescField));

  for (cur_field_num = 0;
       cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    cur_field->bits_per_num = 0;
    cur_field->prev_value = 0;
    io_field_sep_and_is_cons(cur_field, fp, encdesc_gt_xfread);

    if (cur_field->is_cons) {
      read_cons_field_header(cur_field, fp);
    }
    else {
      (void) gt_xfread(&cur_field->is_numeric, sizeof (bool), (size_t) 1, fp);
      if (cur_field->is_numeric) {
        read_numeric_field_header(cur_field, fp);
      }
      else {
        io_field_len_header(cur_field, fp, encdesc_gt_xfread);

        cur_field->bittab = gt_bittab_new(cur_field->len);

        cur_field->data = gt_calloc((size_t) (cur_field->len + 1),
                                    sizeof (char));
        (void) gt_xfread(cur_field->data, sizeof (char),
                         (size_t) cur_field->len, fp);

        read_field_header_bittab(cur_field, fp);

        cur_field->chars = gt_calloc((size_t) (cur_field->max_len + 1),
                                    sizeof (cur_field->chars));

        read_field_char_dists(cur_field, fp);
      }
    }
  }
}
