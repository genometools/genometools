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
#include "extended/io_function_pointers.h"

#define GT_ENCDESC_IO_ONE(element) \
    io_func(&element, sizeof (element), (size_t) 1, fp, err)

static inline int encdesc_header_io_basics(GtEncdesc *encdesc,
                                           FILE *fp,
                                           GtIOFunc io_func,
                                           GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE((encdesc->num_of_descs));
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE((encdesc->num_of_fields));
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE((encdesc->num_of_fields_is_cons));
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE((encdesc->bits_per_field));
  return had_err;
}

static int io_field_sep_and_is_cons(DescField *field,
                                    FILE *fp,
                                    GtIOFunc io_func,
                                    GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE(field->sep);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->is_cons);
  return had_err;
}

static int io_cons_field_header(DescField *field,
                                FILE *fp,
                                GtIOFunc io_func,
                                GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE(field->len);

  if (!had_err) {
    field->data = gt_realloc(field->data,
                           (size_t) (field->len + 1) * sizeof (char));
    had_err = io_func(field->data, sizeof (char), (size_t) field->len, fp, err);
  }
  return had_err;
}

static bool zero_padding_needs_dist(DescField *field)
{
  if (field->has_zero_padding) {
    if (!field->fieldlen_is_const)
      return true;
  }
  return false;
}

static int io_zero_padding(DescField *field, FILE *fp, GtIOFunc io_func,
                           GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE(field->has_zero_padding);
  if (!had_err && field->has_zero_padding) {
    had_err = GT_ENCDESC_IO_ONE(field->fieldlen_is_const);
    if (!had_err) {
      if (field->fieldlen_is_const)
        had_err = GT_ENCDESC_IO_ONE(field->len);
      else {
        had_err = GT_ENCDESC_IO_ONE(field->max_zero);
      }
    }
  }
  return had_err;
}

typedef struct data {
  FILE    *fp;
  GtUword  written_elems;
  GtWord   minimum_element;
  GtError *err;
  int had_err;
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
  if (!this_data->had_err) {
    this_data->written_elems++;
    this_data->had_err = gt_io_error_fwrite_one(symbol, this_data->fp,
                                                this_data->err);
    if (!this_data->had_err)
      this_data->had_err = gt_io_error_fwrite_one(freq, this_data->fp,
                                                  this_data->err);
  }
}

static int io_min_max_calc_ranges(DescField *field,
                                  FILE *fp,
                                  GtIOFunc io_func,
                                  GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE(field->min_value);
  if (!had_err) {
    had_err = GT_ENCDESC_IO_ONE(field->max_value);
    gt_assert(field->min_value <= field->max_value);
  }
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->min_delta);
  if (!had_err) {
    had_err = GT_ENCDESC_IO_ONE(field->max_delta);
    gt_assert(field->min_delta <= field->max_delta);
  }

  if (!had_err) {
    had_err = GT_ENCDESC_IO_ONE(field->is_value_cons);
    if (field->is_value_cons)
      had_err = GT_ENCDESC_IO_ONE(field->global_value);
  }
  if (!had_err) {
    if (field->is_value_cons) {
      gt_assert(field->min_value == field->global_value);
      gt_assert(field->max_value == field->global_value);
    }
    had_err = GT_ENCDESC_IO_ONE(field->is_delta_cons);
  }
  if (!had_err && field->is_delta_cons) {
    had_err = GT_ENCDESC_IO_ONE(field->global_delta);
    gt_assert(field->min_delta == field->global_delta);
    gt_assert(field->max_delta == field->global_delta);
  }

  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->use_delta_coding);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->delta_values_size);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->num_values_size);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->bits_per_num);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->bits_per_value);
  return had_err;
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

static enum iterator_op
encdesc_li_ull_hashmap_iter_count(GT_UNUSED GtWord key,
                                  GT_UNUSED GtUint64 value,
                                  void *data, GT_UNUSED GtError *err)
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
                                                          GtError *err)
{
  int had_err = 0;
  EncsecDistriData *this_data = (EncsecDistriData*) data;
  this_data->written_elems++;
  had_err = gt_io_error_fwrite_one(key, this_data->fp, err);
  if (!had_err)
    had_err = gt_io_error_fwrite_one(value, this_data->fp, err);
  if (!had_err)
    return CONTINUE_ITERATION;
  return STOP_ITERATION;
}

static inline GtHashtable *read_hashmap_distri(GtHashtable *h_table,
                                               GtUword size,
                                               FILE *fp,
                                               GtError *err)
{
  int had_err = 0;
  GtUword idx;
  GtWord symbol;
  GtUint64 freq;

  h_table = li_ull_gt_hashmap_new();
  for (idx = 0; !had_err && idx < size; idx++) {
    had_err = gt_io_error_fread_one(symbol, fp, err);
    if (!had_err)
      had_err = gt_io_error_fread_one(freq, fp, err);
    if (!had_err) {
      gt_assert(li_ull_gt_hashmap_get(h_table, symbol) == 0);
      (void) li_ull_gt_hashmap_add(h_table, symbol, freq);
    }
  }
  if (had_err) {
    li_ull_gt_hashmap_delete(h_table);
    h_table = NULL;
  }
  return h_table;
}

static inline GtHashtable *write_hashmap_distri(GtHashtable *h_table,
                                                GT_UNUSED GtUword size,
                                                EncsecDistriData *data,
                                                GtError *err)
{
  int had_err = 0;
  data->written_elems = 0;
  had_err = li_ull_gt_hashmap_foreach(h_table,
                                      encdesc_li_ull_hashmap_iter_write,
                                      data, err);
  if (!had_err && data->written_elems != size)
    gt_log_log(GT_WU " != " GT_WU, size, data->written_elems);
  if (!had_err)
    gt_assert(data->written_elems == size);
  if (had_err) {
    li_ull_gt_hashmap_delete(h_table);
    h_table = NULL;
  }
  return h_table;
}

static inline GtHashtable *io_hashmap_distri(GtHashtable *h_table,
                                             GtUword size,
                                             EncsecDistriData *data,
                                             GtError *err)
{
  if (h_table != NULL) {
    h_table = write_hashmap_distri(h_table, size, data, err);
  }
  else {
    h_table = read_hashmap_distri(h_table, size, data->fp, err);
  }
  return h_table;
}

static inline GtDiscDistri *read_zero_disc_distri(GtDiscDistri *dist,
                                                  FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword idx, symbol;
  GtUint64 freq;
  GtUword num_of_zero_leaves;
  dist = gt_disc_distri_new();
  had_err = gt_io_error_fread_one(num_of_zero_leaves, fp, err);
  for (idx = 0; !had_err && idx < num_of_zero_leaves; idx++) {
    had_err = gt_io_error_fread_one(symbol, fp, err);
    if (!had_err)
      had_err = gt_io_error_fread_one(freq, fp, err);
    if (!had_err)
      gt_disc_distri_add_multi(dist, symbol, freq);
  }
  if (had_err) {
    gt_disc_distri_delete(dist);
    dist = NULL;
  }
  return dist;
}

static inline GtDiscDistri *write_zero_disc_distri(GtDiscDistri *dist,
                                                   EncsecDistriData *data,
                                                   GtError *err)
{
  GtUword num_of_zero_leaves;
  data->written_elems = 0;
  data->err = err;
  data->had_err = 0;
  gt_disc_distri_foreach(dist, encdesc_distri_iter_count, data);
  if (!data->had_err) {
    num_of_zero_leaves = data->written_elems;
  }
  data->had_err = gt_io_error_fwrite_one(num_of_zero_leaves, data->fp, err);
  if (!data->had_err) {
    data->written_elems = 0;
    gt_disc_distri_foreach(dist, encdesc_distri_iter_write, data);
    gt_assert(data->written_elems == num_of_zero_leaves);
  }
  if (data->had_err) {
    gt_disc_distri_delete(dist);
    dist = NULL;
  }
  return dist;
}

static inline GtDiscDistri *io_zero_disc_distri(GtDiscDistri *dist,
                                                EncsecDistriData *data,
                                                GtError *err)
{
  if (dist == NULL)
    dist = read_zero_disc_distri(dist, data->fp, err);
  else
    dist = write_zero_disc_distri(dist, data, err);
  return dist;
}

static int io_numeric_field_header(DescField *field, FILE *fp, GtIOFunc io_func,
                                   GtError *err)
{
  int had_err = 0;
  bool needs_zero_dist = false,
       needs_delta_dist = false,
       needs_value_dist = false;
  EncsecDistriData data;

  data.fp = fp;

  had_err = io_zero_padding(field, fp, io_func, err);
  if (!had_err) {
    needs_zero_dist = zero_padding_needs_dist(field);
    had_err = io_min_max_calc_ranges(field, fp, io_func, err);
  }
  if (!had_err) {
    numeric_field_check_distri_dependence(field, &needs_delta_dist,
                                          &needs_value_dist);
  }

  if (!had_err && needs_delta_dist) {
    field->delta_values = io_hashmap_distri(field->delta_values,
                                            field->delta_values_size,
                                            &data, err);
    if (field->delta_values == NULL)
      had_err = -1;
  }

  if (!had_err && needs_value_dist) {
    field->num_values = io_hashmap_distri(field->num_values,
                                          field->num_values_size,
                                          &data, err);
  }

  /* TODO change this distribution, to hashtable like the others, this reduces
   * functions */
  if (!had_err && needs_zero_dist) {
    field->zero_count = io_zero_disc_distri(field->zero_count, &data, err);
    if (field->zero_count == NULL)
      had_err = -1;
  }
  return had_err;
}

/* static int read_numeric_field_header(DescField *field,
   FILE *fp, GtError *err)
   {
   int had_err = 0;
   bool needs_zero_dist = false,
   needs_delta_dist = false,
   needs_value_dist = false;
   GtUword num_of_zero_leaves;
   EncsecDistriData data;

   data.fp = fp;
   data.err = err;
   data.had_err = 0;

   had_err = io_zero_padding(field, fp, gt_io_error_fread, err);
   if (!had_err) {
   needs_zero_dist = zero_padding_needs_dist(field);
   had_err = io_min_max_calc_ranges(field, fp, gt_io_error_fread, err);
   }
   if (!had_err) {
   numeric_field_check_distri_dependence(field,
   &needs_delta_dist,
   &needs_value_dist);
   }

   if (!had_err && needs_delta_dist) {
   field->delta_values = io_hashmap_distri(field->delta_values,
   field->delta_values_size,
   &data, err);
   if (field->delta_values)
   had_err = 1;
   }

   if (!had_err && needs_value_dist) {
   field->num_values = io_hashmap_distri(field->num_values,
   field->num_values_size,
   &data, err);
   if (field->num_values)
   had_err = 1;
   }

   if (!had_err && needs_zero_dist) {
   GtUword idx, symbol;
   GtUint64 freq;
   field->zero_count = gt_disc_distri_new();
   had_err = gt_io_error_fread_one(num_of_zero_leaves, fp, err);
   for (idx = 0; !had_err && idx < num_of_zero_leaves; idx++) {
   had_err = gt_io_error_fread_one(symbol, fp, err);
   if (!had_err)
   had_err = gt_io_error_fread_one(freq, fp, err);
   if (!had_err)
   gt_disc_distri_add_multi(field->zero_count, symbol, freq);
   }
   }
   return had_err;
   }

   static int write_numeric_field_header(DescField *field,
   FILE *fp, GtError *err)
   {
   int had_err = 0;
   bool needs_zero_dist = false,
   needs_delta_dist = false,
   needs_value_dist = false;
   GtUword num_of_zero_leaves;
   EncsecDistriData data;

   data.err = err;
   data.fp = fp;
   data.had_err = 0;

had_err = io_zero_padding(field, fp, gt_io_error_fwrite, err);
if (!had_err) {
  needs_zero_dist = zero_padding_needs_dist(field);
  had_err = io_min_max_calc_ranges(field, fp, gt_io_error_fwrite, err);
}

if (!had_err) {
  numeric_field_check_distri_dependence(field,
                                        &needs_delta_dist,
                                        &needs_value_dist);

  if (needs_delta_dist) {
    field->delta_values = io_hashmap_distri(field->delta_values,
                                            field->delta_values_size,
                                            &data, err);
    if (field->delta_values)
      had_err = 1;
  }
}

if (!had_err && needs_value_dist) {
  field->num_values = io_hashmap_distri(field->num_values,
                                        field->num_values_size,
                                        &data, err);
  if (field->num_values)
    had_err = 1;
}

if (!had_err && needs_zero_dist) {
  [>TODO change this distribution, to hashtable like the others, this reduces
    * functions<]
    data.written_elems = 0;
  gt_disc_distri_foreach(field->zero_count,
                         encdesc_distri_iter_count,
                         &data);
  num_of_zero_leaves = data.written_elems;
  had_err = gt_io_error_fwrite_one(num_of_zero_leaves, fp, err);
  data.written_elems = 0;
  gt_disc_distri_foreach(field->zero_count,
                         encdesc_distri_iter_write,
                         &data);
  had_err = data.had_err;
  gt_assert(data.written_elems == num_of_zero_leaves);
}
return had_err;
} */

static int io_field_len_header(DescField *field,
                               FILE *fp,
                               GtIOFunc io_func,
                               GtError *err)
{
  int had_err = 0;
  had_err = GT_ENCDESC_IO_ONE(field->fieldlen_is_const);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->len);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->max_len);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->min_len);
  if (!had_err)
    had_err = GT_ENCDESC_IO_ONE(field->bits_per_len);
  return had_err;
}

static int read_field_header_bittab(DescField *field,
                                    FILE *fp, GtError *err)
{
  int had_err = 0;
  char cc;
  GtUword char_idx,
          num_of_chars = field->len / sizeof (char);
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; char_idx < num_of_chars; char_idx++) {
    had_err = gt_io_error_fread_one(cc, fp, err);
    for (bit_idx = 0; !had_err && bit_idx < sizeof (char); bit_idx++) {
      if (cc & (1 << bit_idx))
        gt_bittab_set_bit(field->bittab,
                          (GtUword) ((sizeof (char) * char_idx) +
                                     bit_idx));
    }
  }
  return had_err;
}

static int write_field_header_bittab(DescField *field, FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword char_idx,
          num_of_chars = field->len / sizeof (char);
  char cc = 0;
  size_t bit_idx;

  if (field->len % sizeof (char) != 0)
    num_of_chars++;

  for (char_idx = 0; !had_err && char_idx < num_of_chars; char_idx++) {
    for (bit_idx = 0; !had_err && bit_idx < sizeof (char); bit_idx++) {
      if (gt_bittab_bit_is_set(field->bittab,
                               (GtUword) ((sizeof (char) * char_idx) +
                                          bit_idx)))
        cc |= 1 << bit_idx;
    }
    had_err = gt_io_error_fwrite_one(cc, fp, err);
    cc = 0;
  }
  return had_err;
}

/* TODO combine field_char_dist */
static int write_field_char_dists(DescField *field,
                                  FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword char_idx,
          distr_len;
  EncsecDistriData data;
  data.fp = fp;
  data.minimum_element = 0;
  data.written_elems = 0;
  data.had_err = 0;
  data.err = err;

  for (char_idx = 0; !had_err && char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      distr_len = get_hashmap_distri_size(field->chars[char_idx]);
      had_err = gt_io_error_fwrite_one(distr_len, fp, err);
      if (!had_err)
        field->chars[char_idx] = io_hashmap_distri(field->chars[char_idx],
                                                   distr_len,
                                                   &data, err);
      if (field->chars[char_idx] == NULL)
        had_err = 1;
    }
  }
  return had_err;
}

static int read_field_char_dists(DescField *field,
                                 FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword char_idx,
          distr_len;
  EncsecDistriData data;

  data.err = err;
  data.had_err = 0;
  data.fp = fp;

  for (char_idx = 0; !had_err && char_idx < field->max_len; char_idx++) {
    if (char_idx >= field->len ||
        !gt_bittab_bit_is_set(field->bittab, char_idx)) {

      had_err = gt_io_error_fread_one(distr_len, fp, err);
      field->chars[char_idx] = io_hashmap_distri(field->chars[char_idx],
                                                 distr_len,
                                                 &data, err);
      if (field->chars[char_idx] == NULL)
        had_err = 1;
    }
  }
  return had_err;
}

int encdesc_write_header(GtEncdesc *encdesc, FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword cur_field_num;
  DescField *cur_field;

  had_err = encdesc_header_io_basics(encdesc, fp, gt_io_error_fwrite, err);

  for (cur_field_num = 0;
       !had_err && cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    had_err = io_field_sep_and_is_cons(cur_field, fp, gt_io_error_fwrite, err);

    if (!had_err && cur_field->is_cons) {
      had_err = io_cons_field_header(cur_field, fp, gt_io_error_fwrite, err);
    }
    else if (!had_err) {
      had_err = gt_io_error_fwrite_one(cur_field->is_numeric, fp, err);
      if (!had_err && cur_field->is_numeric) {
        had_err = io_numeric_field_header(cur_field, fp, gt_io_error_fwrite,
                                          err);
      }
      else if (!had_err) {

        had_err = io_field_len_header(cur_field, fp, gt_io_error_fwrite, err);

        if (!had_err)
          had_err = gt_io_error_fwrite(cur_field->data, sizeof (char),
                                       (size_t) cur_field->len, fp, err);

        if (!had_err)
          had_err = write_field_header_bittab(cur_field, fp, err);

        if (!had_err)
          had_err = write_field_char_dists(cur_field, fp, err);
      }
    }
  }
  return had_err;
}

int encdesc_read_header(GtEncdesc *encdesc, FILE *fp, GtError *err)
{
  int had_err = 0;
  GtUword cur_field_num;
  DescField *cur_field;

  had_err = encdesc_header_io_basics(encdesc, fp, gt_io_error_fread, err);

  if (!had_err)
    encdesc->fields = gt_calloc((size_t) encdesc->num_of_fields,
                                sizeof (DescField));

  for (cur_field_num = 0;
       !had_err && cur_field_num < encdesc->num_of_fields;
       cur_field_num++) {
    cur_field = &encdesc->fields[cur_field_num];
    cur_field->bits_per_num = 0;
    cur_field->prev_value = 0;
    had_err = io_field_sep_and_is_cons(cur_field, fp, gt_io_error_fread, err);

    if (!had_err && cur_field->is_cons) {
      had_err = io_cons_field_header(cur_field, fp, gt_io_error_fread, err);
    }
    else if (!had_err) {
      had_err = gt_io_error_fread_one(cur_field->is_numeric, fp, err);
      if (!had_err && cur_field->is_numeric) {
        had_err = io_numeric_field_header(cur_field, fp, gt_io_error_fread,
                                          err);
      }
      else if (!had_err) {
        had_err = io_field_len_header(cur_field, fp, gt_io_error_fread, err);

        if (!had_err) {
          cur_field->bittab = gt_bittab_new(cur_field->len);

          cur_field->data = gt_calloc((size_t) (cur_field->len + 1),
                                      sizeof (char));
          had_err = gt_io_error_fread(cur_field->data, sizeof (char),
                                      (size_t) cur_field->len, fp, err);

        }
        if (!had_err)
          had_err = read_field_header_bittab(cur_field, fp, err);

        if (!had_err) {
          cur_field->chars = gt_calloc((size_t) (cur_field->max_len + 1),
                                       sizeof (cur_field->chars));

          had_err = read_field_char_dists(cur_field, fp, err);
        }
      }
    }
  }
  return had_err;
}
