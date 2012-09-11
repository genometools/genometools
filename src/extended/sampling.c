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

#include "core/assert_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/xansi_api.h"
#include "core/unused_api.h"
#include "extended/sampling.h"

typedef enum {
  GT_SAMPLING_REGULAR,
  GT_SAMPLING_PAGES,
} GtSamplingMethod;

struct GtSampling
{
  GtSamplingMethod method;
  unsigned long numofsamples,
                sampling_rate,
                arraysize,
                *page_sampling,
                current_sample_num,
                current_sample_elementnum;
  size_t *samplingtab;
  long pagesize;
};

static void sampling_init_sampling(GtSampling *sampling,
                                   unsigned long rate)
{
  sampling->numofsamples = 1UL;
  sampling->arraysize = 10UL;
  sampling->sampling_rate = rate;
  sampling->current_sample_elementnum = 0;
  sampling->current_sample_num = 0;
  sampling->pagesize = sysconf((int) _SC_PAGESIZE);
}

GtSampling *gt_sampling_new_regular(unsigned long rate, off_t first_offset)
{
  GtSampling *sampling = gt_malloc(sizeof (*sampling));
  sampling->method = GT_SAMPLING_REGULAR;
  gt_assert(rate != 0);

  sampling_init_sampling(sampling, rate);

  gt_assert(first_offset % sampling->pagesize == 0);

  sampling->page_sampling = NULL;
  sampling->samplingtab = gt_malloc((size_t) sampling->arraysize *
                                    sizeof (*sampling->samplingtab));
  gt_safe_assign(sampling->samplingtab[0], first_offset);
  return sampling;
}

GtSampling *gt_sampling_new_page(unsigned long rate, off_t first_offset)
{
  GtSampling *sampling = gt_sampling_new_regular(rate, first_offset);
  sampling->method = GT_SAMPLING_PAGES;
  gt_assert(rate != 0);

  sampling->page_sampling = gt_malloc((size_t) sampling->arraysize *
                                      sizeof (*sampling->page_sampling));
  sampling->page_sampling[0] = 0;
  return sampling;
}

static inline void sampling_gt_xfwrite(void *ptr,
                                       size_t size,
                                       size_t nmemb,
                                       FILE *stream)
{
  gt_xfwrite((const void*) ptr, size, nmemb, stream);
}

static inline void sampling_gt_xfread(void *ptr,
                                      size_t size,
                                      size_t nmemb,
                                      FILE *stream)
{
  GT_UNUSED size_t read;
  read = gt_xfread(ptr, size, nmemb, stream);
  gt_assert(read == nmemb);
}

typedef void (*SamplingIOFunc)(void *ptr,
                               size_t size,
                               size_t nmemb,
                               FILE *stream);

#define SAMPLING_IO_ONE(element, fp)                    \
  do {                                                  \
    io_func(&element, sizeof (element), (size_t) 1, fp);\
  } while (false)

static inline void sampling_io_page_sampling(GtSampling *sampling,
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

static inline void sampling_io_samplingtab(GtSampling *sampling,
                                    FILE *fp,
                                    SamplingIOFunc io_func)
{
  io_func(sampling->samplingtab,
          sizeof (*sampling->samplingtab),
          (size_t) sampling->numofsamples,
          fp);
}

static inline void sampling_io_header(GtSampling *sampling,
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
static void sampling_io_header_samplingtab(GtSampling *sampling,
                                           FILE *fp,
                                           SamplingIOFunc io_func)
{
  sampling_io_header(sampling, fp, io_func);
  gt_assert(sampling->method == GT_SAMPLING_REGULAR ||
            sampling->method == GT_SAMPLING_PAGES);

  if (sampling->samplingtab == NULL) {
    sampling->arraysize = sampling->numofsamples;
    sampling->samplingtab = gt_malloc((size_t) sampling->arraysize *
                                        sizeof (*sampling->samplingtab));
  }
  sampling_io_samplingtab(sampling, fp, io_func);
}

void gt_sampling_write(GtSampling *sampling, FILE *fp)
{
  gt_assert(sampling);
  gt_assert(fp);

  sampling_io_header_samplingtab(sampling, fp, sampling_gt_xfwrite);
  if (sampling->method == GT_SAMPLING_PAGES)
    sampling_io_page_sampling(sampling, fp, sampling_gt_xfwrite);
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
  sampling->pagesize = sysconf((int) _SC_PAGESIZE);

  sampling_io_header_samplingtab(sampling, fp, sampling_gt_xfread);
  if (sampling->method == GT_SAMPLING_PAGES)
    sampling_io_page_sampling(sampling, fp, sampling_gt_xfread);
  gt_assert(sampling->arraysize == sampling->numofsamples);

  return sampling;
}

static void get_regular_page(GtSampling *sampling,
                             unsigned long element_num,
                             unsigned long *sampled_element,
                             size_t *position)
{
  sampling->current_sample_num = element_num/sampling->sampling_rate;

  *sampled_element =
    sampling->current_sample_elementnum =
    sampling->current_sample_num * sampling->sampling_rate;

  *position = sampling->samplingtab[sampling->current_sample_num];
}

static void get_pagewise_page(GtSampling *sampling,
                              unsigned long element_num,
                              unsigned long *sampled_element,
                              size_t *position)
{
  unsigned long start = 0,
                end, middle;

  gt_assert(sampling->numofsamples != 0);
  end = sampling->numofsamples - 1;
  middle = (end - start) >> 1;
  while (start < end) {
    if (sampling->page_sampling[middle] == element_num)
      break;
    else {
      if (sampling->page_sampling[middle] > element_num) {
        end = middle - 1;
        middle = start + ((end - start) >> 1);
      }
      else {
        if (sampling->page_sampling[middle + 1] > element_num)
          break;
        else {
          start = middle + 1;
          middle = start + ((end - start) >> 1);
        }
      }
    }
  }

  *sampled_element =
    sampling->current_sample_elementnum =
    sampling->page_sampling[middle];

  sampling->current_sample_num = middle;

  *position = sampling->samplingtab[middle];
}

void gt_sampling_get_page(GtSampling *sampling,
                          unsigned long element_num,
                          unsigned long *sampled_element,
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

unsigned long gt_sampling_get_current_elementnum(GtSampling *sampling)
{
  return sampling->current_sample_elementnum;
}

unsigned long gt_sampling_get_next_elementnum(GtSampling *sampling)
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
      return GT_UNDEF_ULONG;
  }
}

int gt_sampling_get_next_sample(GtSampling *sampling,
                                unsigned long *sampled_element,
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

unsigned long gt_sampling_get_rate(GtSampling *sampling)
{
  gt_assert(sampling);
  return sampling->sampling_rate;
}

void gt_sampling_add_sample(GtSampling *sampling,
                            size_t position,
                            unsigned long element_num)
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

bool          gt_sampling_is_next_element_sample(
                                          GtSampling *sampling,
                                          unsigned long pages_written,
                                          unsigned long elements_written,
                                          unsigned long elem_bit_size,
                                          unsigned long free_pagespace_bitsize)
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
