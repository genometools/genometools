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

#include <stdio.h>

#include "core/ma_api.h"
#include "core/log_api.h"
#include "extended/encdesc_header_io.h"
#include "extended/encdesc_rep.h"

static inline void encdesc_xfwrite(void *ptr,
                                   size_t size,
                                   size_t nmemb,
                                   FILE *stream)
{
  if (nmemb != fwrite((const void*) ptr, size, nmemb, stream)) {
    perror("gt_sampling_xfwrite could not write to file");
    exit(EXIT_FAILURE);
  }
}

static inline void encdesc_xfread(void *ptr,
                                  size_t size,
                                  size_t nmemb,
                                  FILE *stream)
{
  if (nmemb != fread(ptr, size, nmemb, stream)) {
    gt_assert(feof(stream) == 0);
    if (ferror(stream) != 0)
      perror("gt_sampling_xfread could not read from file");
    exit(EXIT_FAILURE);
  };
}

#define GT_ENCDESC_IO_ONE(element, fp)                  \
    io_func(&element, sizeof (element), (size_t) 1, fp)

#define GT_ENCDESC_WRITE_ONE(element, fp)               \
    encdesc_xfwrite(&element, sizeof (element), (size_t) 1, fp)

#define GT_ENCDESC_READ_ONE(element, fp)                \
    encdesc_xfread(&element, sizeof (element), (size_t) 1, fp)

typedef void (*EncdescIOFunc)(void *ptr,
                              size_t size,
                              size_t nmemb,
                              FILE *stream);

static inline void encdesc_header_io_basics(GtEncdesc *encdesc,
                                            FILE *fp,
                                            EncdescIOFunc io_func)
{
  GT_ENCDESC_IO_ONE((encdesc->num_of_descs), fp);
  GT_ENCDESC_IO_ONE((encdesc->num_of_fields), fp);
  GT_ENCDESC_IO_ONE((encdesc->num_of_fields_is_cons), fp);
  GT_ENCDESC_IO_ONE((encdesc->bits_per_field), fp);
}

static void io_field_sep_and_is_cons(DescField *field,
                                     FILE *fp,
                                     EncdescIOFunc io_func)
{
  GT_ENCDESC_IO_ONE(field->sep, fp);
  GT_ENCDESC_IO_ONE(field->is_cons, fp);
}

static void write_cons_field_header(DescField *field,
                                    FILE *fp)
{
  GT_ENCDESC_WRITE_ONE(field->len, fp);

  encdesc_xfwrite(field->data, sizeof (char), (size_t) field->len, fp);
}

static void read_cons_field_header(DescField *field,
                                   FILE *fp)
{
  GT_ENCDESC_READ_ONE(field->len, fp);
  field->data = gt_calloc((size_t) (field->len + 1), sizeof (char));
  encdesc_xfread(field->data, sizeof (char), (size_t) field->len, fp);
}

static bool io_zero_padding_needs_dist(DescField *field,
                                       FILE *fp,
                                       EncdescIOFunc io_func)
{
  GT_ENCDESC_IO_ONE(field->has_zero_padding, fp);
  if (field->has_zero_padding) {
    GT_ENCDESC_IO_ONE(field->fieldlen_is_const, fp);
    if (field->fieldlen_is_const)
      GT_ENCDESC_IO_ONE(field->len, fp);
    else {
      GT_ENCDESC_IO_ONE(field->max_zero, fp);
      return true;
    }
  }
  return false;
}

typedef struct data {
  FILE         *fp;
  GtUword written_elems;
  GtWord          minimum_element;
} EncsecDistriData;

static void encdesc_distri_iter_count(GT_UNUSED GtUword symbol,
                                      GT_UNUSED GtUint64 freq,
                                      void *data)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
}

static void encdesc_distri_iter_write(GtUword symbol,
                                      GtUint64 freq,
                                      void *data)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
  GT_ENCDESC_WRITE_ONE(symbol, this_data->fp);
  GT_ENCDESC_WRITE_ONE(freq, this_data->fp);

}

static void io_min_max_calc_ranges(DescField *field,
                       FILE *fp,
                       EncdescIOFunc io_func)
{
  GT_ENCDESC_IO_ONE(field->min_value, fp);
  GT_ENCDESC_IO_ONE(field->max_value, fp);
  gt_assert(field->min_value <= field->max_value);
  GT_ENCDESC_IO_ONE(field->min_delta, fp);
  GT_ENCDESC_IO_ONE(field->max_delta, fp);
  gt_assert(field->min_delta <= field->max_delta);

  GT_ENCDESC_IO_ONE(field->is_value_cons, fp);
  if (field->is_value_cons)
  GT_ENCDESC_IO_ONE(field->global_value, fp);
  if (field->is_value_cons) {
    gt_assert(field->min_value == field->global_value);
    gt_assert(field->max_value == field->global_value);
  }
  GT_ENCDESC_IO_ONE(field->is_delta_cons, fp);
  if (field->is_delta_cons)
  GT_ENCDESC_IO_ONE(field->global_delta, fp);
  if (field->is_delta_cons) {
    gt_assert(field->min_delta == field->global_delta);
    gt_assert(field->max_delta == field->global_delta);
  }

  GT_ENCDESC_IO_ONE(field->use_delta_coding, fp);
  GT_ENCDESC_IO_ONE(field->delta_values_size, fp);
  GT_ENCDESC_IO_ONE(field->num_values_size, fp);
  GT_ENCDESC_IO_ONE(field->bits_per_num, fp);
  GT_ENCDESC_IO_ONE(field->bits_per_value, fp);
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
      gt_log_log("delta_values_size: "GT_WU"", field->delta_values_size);
    }
  }
  else {
    if (!field->is_value_cons && field->num_values_size > 0 &&
        field->num_values_size <= GT_ENCDESC_MAX_NUM_VAL_HUF) {
      *needs_value_dist = true;
      field->use_hc = true;
      gt_log_log("num_values_size: "GT_WU"", field->num_values_size);
    }
  }
}

static enum iterator_op encdesc_li_ull_hashmap_iter_count(
                                             GT_UNUSED GtWord key,
                                             GT_UNUSED GtUint64 value,
                                             void *data,
                                             GT_UNUSED GtError *err)
{
  GtUword *count = (GtUword*) data;
  (void) (*count)++;
  /* always continue since we cannot fail */
  return CONTINUE_ITERATION;
}

static GtUword get_hashmap_distri_size(GtHashtable *h_table)
{
  GtUword count = 0;
  /* error is NULL as encdesc_li_ull_hashmap_iter_count() is sane */
  (void) li_ull_gt_hashmap_foreach(h_table,
                                   encdesc_li_ull_hashmap_iter_count,
                                   &count,
                                   NULL);
  return count;
}

static enum iterator_op encdesc_li_ull_hashmap_iter_write(
                                                       GtWord key,
                                                       GtUint64 value,
                                                       void *data,
                                                       GT_UNUSED GtError *err)
{
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
  GT_ENCDESC_WRITE_ONE(key, this_data->fp);
  GT_ENCDESC_WRITE_ONE(value, this_data->fp);
  /* always continue since we exit on write error anyway */
  return CONTINUE_ITERATION;
}

static void write_hashmap_distri(EncsecDistriData *data,
                                 GtHashtable *h_table,
                                 GT_UNUSED GtUword size)
{
  data->written_elems = 0;
  /* error is NULL as encdesc_li_ull_hashmap_iter_write() is sane */
  (void) li_ull_gt_hashmap_foreach(h_table,
                                   encdesc_li_ull_hashmap_iter_write,
                                   data,
                                   NULL);
  if (data->written_elems != size)
    gt_log_log(""GT_WU" != "GT_WU"", size, data->written_elems);
  gt_assert(data->written_elems == size);
}

static void read_hashmap_distri(GtUword size,
                                GtHashtable *h_table,
                                FILE *fp)
{
  GtUword idx;
  GtWord symbol;
  GtUint64 freq;

  gt_assert(h_table != NULL);
  for (idx = 0; idx < size; idx++) {
    GT_ENCDESC_READ_ONE(symbol, fp);
    GT_ENCDESC_READ_ONE(freq, fp);
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
  GtUword num_of_zero_leaves;
  /* EncdescHuffDist huffdist; */

  needs_zero_dist = io_zero_padding_needs_dist(field,
                                               fp,
                                               encdesc_xfread);
  io_min_max_calc_ranges(field,
                         fp,
                         encdesc_xfread);

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
    GtUword idx,
                  symbol;
    GtUint64 freq;
    field->zero_count = gt_disc_distri_new();
    GT_ENCDESC_READ_ONE(num_of_zero_leaves, fp);
    for (idx = 0; idx < num_of_zero_leaves; idx++) {
      GT_ENCDESC_READ_ONE(symbol, fp);
      GT_ENCDESC_READ_ONE(freq, fp);
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
  GtUword num_of_zero_leaves;
  EncsecDistriData data;

  data.fp = fp;

  needs_zero_dist = io_zero_padding_needs_dist(field,
                                               fp,
                                               encdesc_xfwrite);
  io_min_max_calc_ranges(field,
                         fp,
                         encdesc_xfwrite);

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
    GT_ENCDESC_WRITE_ONE(num_of_zero_leaves, fp);
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
  GT_ENCDESC_IO_ONE(field->fieldlen_is_const, fp);
  GT_ENCDESC_IO_ONE(field->len, fp);
  GT_ENCDESC_IO_ONE(field->max_len, fp);
  GT_ENCDESC_IO_ONE(field->min_len, fp);
  GT_ENCDESC_IO_ONE(field->bits_per_len, fp);
}

static void read_field_header_bittab(DescField *field,
                                     FILE *fp)
{
  GtUword char_idx,
                num_of_chars = field->len / sizeof (char);
  char cc;
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; char_idx < num_of_chars; char_idx++) {
    GT_ENCDESC_READ_ONE(cc, fp);
    for (bit_idx = 0; bit_idx < sizeof (char); bit_idx++) {
      if (cc & (1 << bit_idx))
        gt_bittab_set_bit(field->bittab,
                          (GtUword) ((sizeof (char) * char_idx) +
                                           bit_idx));
    }
  }
}

static void write_field_header_bittab(DescField *field, FILE *fp)
{
  GtUword char_idx,
                num_of_chars = field->len / sizeof (char);
  char cc = 0;
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; char_idx < num_of_chars; char_idx++) {
    for (bit_idx = 0; bit_idx < sizeof (char); bit_idx++) {
      if (gt_bittab_bit_is_set(field->bittab,
                               (GtUword) ((sizeof (char) * char_idx) +
                                                bit_idx)))
        cc |= 1 << bit_idx;
    }
    GT_ENCDESC_WRITE_ONE(cc, fp);
    cc = 0;
  }
}

/* TODO combine field_char_dist */
static void write_field_char_dists(DescField *field,
                                   FILE *fp)
{
  GtUword char_idx,
                distr_len;
  EncsecDistriData data;
  data.fp = fp;
  data.minimum_element = 0;
  data.written_elems = 0;

  for (char_idx = 0; char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      distr_len = get_hashmap_distri_size(field->chars[char_idx]);
      GT_ENCDESC_WRITE_ONE(distr_len, fp);
      write_hashmap_distri(&data,
                           field->chars[char_idx],
                           distr_len);
    }
  }
}

static void read_field_char_dists(DescField *field,
                                  FILE *fp)
{
  GtUword char_idx,
                distr_len;

  for (char_idx = 0; char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      GT_ENCDESC_READ_ONE(distr_len, fp);
      field->chars[char_idx] = li_ull_gt_hashmap_new();
      read_hashmap_distri(distr_len,
                          field->chars[char_idx],
                          fp);
    }
  }
}

void encdesc_write_header(GtEncdesc *encdesc, FILE *fp)
{
  GtUword cur_field_num;
  DescField *cur_field;

  encdesc_header_io_basics(encdesc, fp, encdesc_xfwrite);

  for (cur_field_num = 0;
       cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    io_field_sep_and_is_cons(cur_field, fp, encdesc_xfwrite);

    if (cur_field->is_cons) {
      write_cons_field_header(cur_field, fp);
    }
    else {
      GT_ENCDESC_WRITE_ONE(cur_field->is_numeric, fp);
      if (cur_field->is_numeric) {
        write_numeric_field_header(cur_field, fp);
      }
      else {

        io_field_len_header(cur_field, fp, encdesc_xfwrite);

        encdesc_xfwrite(cur_field->data, sizeof (char),
                        (size_t) cur_field->len, fp);

        write_field_header_bittab(cur_field, fp);

        write_field_char_dists(cur_field, fp);
      }
    }
  }
}

void encdesc_read_header(GtEncdesc *encdesc, FILE *fp)
{
  GtUword cur_field_num;
  DescField *cur_field;

  encdesc_header_io_basics(encdesc, fp, encdesc_xfread);

  encdesc->fields = gt_calloc((size_t) encdesc->num_of_fields,
                              sizeof (DescField));

  for (cur_field_num = 0;
       cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    cur_field->bits_per_num = 0;
    cur_field->prev_value = 0;
    io_field_sep_and_is_cons(cur_field, fp, encdesc_xfread);

    if (cur_field->is_cons) {
      read_cons_field_header(cur_field, fp);
    }
    else {
      GT_ENCDESC_READ_ONE(cur_field->is_numeric, fp);
      if (cur_field->is_numeric) {
        read_numeric_field_header(cur_field, fp);
      }
      else {
        io_field_len_header(cur_field, fp, encdesc_xfread);

        cur_field->bittab = gt_bittab_new(cur_field->len);

        cur_field->data = gt_calloc((size_t) (cur_field->len + 1),
                                    sizeof (char));
        encdesc_xfread(cur_field->data, sizeof (char),
                       (size_t) cur_field->len, fp);

        read_field_header_bittab(cur_field, fp);

        cur_field->chars = gt_calloc((size_t) (cur_field->max_len + 1),
                                    sizeof (cur_field->chars));

        read_field_char_dists(cur_field, fp);
      }
    }
  }
}
