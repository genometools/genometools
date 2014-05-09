/*
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
#include <unistd.h>
#endif

#include <stdio.h>

#include "core/assert_api.h"
#include "core/compat.h"
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/sampling.h"

typedef enum {
  GT_SAMPLING_REGULAR,
  GT_SAMPLING_PAGES,
} GtSamplingMethod;

struct GtSampling
{
  size_t          *samplingtab;
  GtUword    arraysize,
                   current_sample_elementnum,
                   current_sample_num,
                   numofsamples,
                   pagesize,
                   sampling_rate,
                  *page_sampling;
  GtSamplingMethod method;
};

static void gt_sampling_init_sampling(GtSampling *sampling,
                                   GtUword rate)
{
  sampling->numofsamples = 1UL;
  sampling->arraysize = 10UL;
  sampling->sampling_rate = rate;
  sampling->current_sample_elementnum = 0;
  sampling->current_sample_num = 0;
  sampling->pagesize = gt_pagesize();
}

GtSampling *gt_sampling_new_regular(GtUword rate, off_t first_offset)
{
  GtSampling *sampling = gt_malloc(sizeof (*sampling));
  sampling->method = GT_SAMPLING_REGULAR;
  gt_assert(rate != 0);

  gt_sampling_init_sampling(sampling, rate);

  gt_assert(first_offset % sampling->pagesize == 0);

  sampling->page_sampling = NULL;
  sampling->samplingtab = gt_malloc((size_t) sampling->arraysize *
                                    sizeof (*sampling->samplingtab));
  gt_safe_assign(sampling->samplingtab[0], first_offset);
  return sampling;
}

GtSampling *gt_sampling_new_page(GtUword rate, off_t first_offset)
{
  GtSampling *sampling = gt_sampling_new_regular(rate, first_offset);
  sampling->method = GT_SAMPLING_PAGES;
  gt_assert(rate != 0);

  sampling->page_sampling = gt_malloc((size_t) sampling->arraysize *
                                      sizeof (*sampling->page_sampling));
  sampling->page_sampling[0] = 0;
  return sampling;
}

static inline void gt_sampling_xfwrite(void *ptr,
                                       size_t size,
                                       size_t nmemb,
                                       FILE *stream)
{
  if (nmemb != fwrite((const void*) ptr, size, nmemb, stream)) {
    perror("gt_sampling_xfwrite could not write to file");
    exit(EXIT_FAILURE);
  }
}

static inline void gt_sampling_xfread(void *ptr,
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

typedef void (*SamplingIOFunc)(void *ptr,
                               size_t size,
                               size_t nmemb,
                               FILE *stream);

#define SAMPLING_IO_ONE(element, fp)                   \
    io_func(&element, sizeof (element), (size_t) 1, fp)\

static inline void gt_sampling_io_page_sampling(GtSampling *sampling,
                                             FILE *fp,
                                             SamplingIOFunc io_func)
{
  if (sampling->page_sampling == NULL) {
    sampling->page_sampling = gt_malloc((size_t) sampling->arraysize *
                                        sizeof (*sampling->page_sampling));
  }
  io_func(sampling->page_sampling,
          sizeof (*sampling->page_sampling),
          (size_t) sampling->numofsamples,
          fp);
}

static inline void gt_sampling_io_samplingtab(GtSampling *sampling,
                                           FILE *fp,
                                           SamplingIOFunc io_func)
{
  io_func(sampling->samplingtab,
          sizeof (*sampling->samplingtab),
          (size_t) sampling->numofsamples,
          fp);
}

static inline void gt_sampling_io_header(GtSampling *sampling,
                                      FILE *fp,
                                      SamplingIOFunc io_func)
{
  SAMPLING_IO_ONE(sampling->numofsamples, fp);
  gt_assert(sampling->numofsamples != 0);
  SAMPLING_IO_ONE(sampling->method, fp);
  SAMPLING_IO_ONE(sampling->sampling_rate, fp);
  gt_assert(sampling->sampling_rate != 0);
}

/* TODO: add checksums for data */
static void gt_sampling_io_header_samplingtab(GtSampling *sampling,
                                           FILE *fp,
                                           SamplingIOFunc io_func)
{
  gt_sampling_io_header(sampling, fp, io_func);
  gt_assert(sampling->method == GT_SAMPLING_REGULAR ||
            sampling->method == GT_SAMPLING_PAGES);

  if (sampling->samplingtab == NULL) {
    sampling->arraysize = sampling->numofsamples;
    sampling->samplingtab = gt_malloc((size_t) sampling->arraysize *
                                        sizeof (*sampling->samplingtab));
  }
  gt_sampling_io_samplingtab(sampling, fp, io_func);
}

void gt_sampling_write(GtSampling *sampling, FILE *fp)
{
  gt_assert(sampling);
  gt_assert(fp);

  gt_sampling_io_header_samplingtab(sampling, fp, gt_sampling_xfwrite);
  if (sampling->method == GT_SAMPLING_PAGES)
    gt_sampling_io_page_sampling(sampling, fp, gt_sampling_xfwrite);
}

GtSampling *gt_sampling_read(FILE *fp)
{
  GtSampling *sampling;

  gt_assert(fp);

  sampling = gt_malloc(sizeof (*sampling));
  sampling->samplingtab = NULL;
  sampling->page_sampling = NULL;
  sampling->current_sample_num =
    sampling->current_sample_elementnum = 0;
  sampling->pagesize = gt_pagesize();

  gt_sampling_io_header_samplingtab(sampling, fp, gt_sampling_xfread);
  if (sampling->method == GT_SAMPLING_PAGES)
    gt_sampling_io_page_sampling(sampling, fp, gt_sampling_xfread);
  gt_assert(sampling->arraysize == sampling->numofsamples);

  return sampling;
}

static void get_regular_page(GtSampling *sampling,
                             GtUword element_num,
                             GtUword *sampled_element,
                             size_t *position)
{
  sampling->current_sample_num = element_num/sampling->sampling_rate;

  *sampled_element =
    sampling->current_sample_elementnum =
    sampling->current_sample_num * sampling->sampling_rate;

  *position = sampling->samplingtab[sampling->current_sample_num];
}

static void get_pagewise_page(GtSampling *sampling,
                              GtUword element_num,
                              GtUword *sampled_element,
                              size_t *position)
{
  GtWord start = (GtWord) -1,
         end, middle;

  gt_assert(sampling->numofsamples != 0);
  /* should not overflow, because this is a small table indexing into a larger
     one. */
  gt_safe_assign(end, sampling->numofsamples);
  middle = GT_DIV2(end);
  while (end - start > (GtWord) 1) {
    if (element_num < sampling->page_sampling[middle]) {
      start = middle;
    }
    else {
      end = middle;
    }
    middle = start + GT_DIV2(end - start);
  }
  if (middle < 0) {
    middle = 0;
  }
  *sampled_element =
    sampling->current_sample_elementnum =
    sampling->page_sampling[middle];

  sampling->current_sample_num = (GtUword) middle;

  *position = sampling->samplingtab[middle];
}

void gt_sampling_get_page(GtSampling *sampling,
                          GtUword element_num,
                          GtUword *sampled_element,
                          size_t *position)
{
  gt_assert(sampling != NULL);
  gt_assert(sampled_element != NULL);
  gt_assert(position != NULL);

  switch (sampling->method) {
    case GT_SAMPLING_REGULAR:
      get_regular_page(sampling, element_num, sampled_element, position);
      break;

    case GT_SAMPLING_PAGES:
      get_pagewise_page(sampling, element_num, sampled_element, position);
      break;
  }
}

GtUword gt_sampling_get_current_elementnum(GtSampling *sampling)
{
  return sampling->current_sample_elementnum;
}

GtUword gt_sampling_get_next_elementnum(GtSampling *sampling)
{
  gt_assert(sampling->arraysize == sampling->numofsamples);
  gt_assert(sampling->current_sample_num < sampling->numofsamples);
  if (sampling->current_sample_num + 1 == sampling->numofsamples)
    return 0;
  gt_assert((sampling->current_sample_num + 1) < sampling->arraysize);
  switch (sampling->method) {
    case GT_SAMPLING_REGULAR:
      return sampling->current_sample_elementnum + sampling->sampling_rate;
    case GT_SAMPLING_PAGES:
      return sampling->page_sampling[sampling->current_sample_num + 1];
    default:
      return GT_UNDEF_UWORD;
  }
}

int gt_sampling_get_next_sample(GtSampling *sampling,
                                GtUword *sampled_element,
                                size_t *position)
{
  enum state {
    ERROR = -1,
    END,
    SUCCESS
  };
  enum state status = END;

  if (sampling->current_sample_num + 1 == sampling->numofsamples) {
    sampling->current_sample_num = 0;
    *sampled_element = sampling->current_sample_elementnum = 0;
  }
  else {
    status = SUCCESS;
    sampling->current_sample_num++;
    switch (sampling->method) {
      case GT_SAMPLING_REGULAR:
        *sampled_element =
          sampling->current_sample_elementnum += sampling->sampling_rate;
        break;
      case GT_SAMPLING_PAGES:
        *sampled_element = sampling->current_sample_elementnum =
          sampling->page_sampling[sampling->current_sample_num];
        break;
      default:
        status = ERROR;
    }
  }
  if (status != ERROR)
    *position = sampling->samplingtab[sampling->current_sample_num];
  return (int) status;
}

bool gt_sampling_is_regular(GtSampling *sampling)
{
  gt_assert(sampling);
  return sampling->method == GT_SAMPLING_REGULAR;
}

GtUword gt_sampling_get_rate(GtSampling *sampling)
{
  gt_assert(sampling);
  return sampling->sampling_rate;
}

void gt_sampling_add_sample(GtSampling *sampling,
                            size_t position,
                            GtUword element_num)
{
  gt_assert(sampling);
  gt_assert(sampling->samplingtab);

  sampling->numofsamples++;

  if (sampling->numofsamples == sampling->arraysize) {
    sampling->arraysize += sampling->arraysize/100 + 10;
    sampling->samplingtab = gt_realloc(sampling->samplingtab,
                                       (size_t) sampling->arraysize *
                                         sizeof (*sampling->samplingtab));
    if (sampling->method == GT_SAMPLING_PAGES) {
      gt_assert(sampling->page_sampling);
      sampling->page_sampling = gt_realloc(sampling->page_sampling,
                                           (size_t) sampling->arraysize *
                                             sizeof (*sampling->page_sampling));
    }
  }
  if (sampling->method == GT_SAMPLING_PAGES)
    sampling->page_sampling[sampling->numofsamples -1] = element_num;
  else
    gt_assert(element_num % sampling->sampling_rate == 0);
  sampling->samplingtab[sampling->numofsamples -1] = position;
}

bool gt_sampling_is_next_element_sample(GtSampling *sampling,
                                        GtUword pages_written,
                                        GtUword elements_written,
                                        GtUword elem_bit_size,
                                        GtUword free_pagespace_bitsize)
{
  if (sampling->method == GT_SAMPLING_REGULAR)
    return elements_written >= sampling->sampling_rate;
  else {
    if (pages_written >= sampling->sampling_rate) {
      return free_pagespace_bitsize < elem_bit_size;
    }
  }
  return false;
}

void gt_sampling_delete(GtSampling *sampling)
{
  if (!sampling) return;
  gt_free(sampling->samplingtab);
  gt_free(sampling->page_sampling);
  gt_free(sampling);
}
