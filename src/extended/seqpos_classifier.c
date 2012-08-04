/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#include <limits.h>
#include <string.h>
#include "core/log_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/gff3_in_stream_api.h"
#include "extended/seqpos_classifier.h"

struct GtSeqposClassifier
{
  GtNodeStream *annotation_stream;
  GtGenomeNode *gn;
  GtFeatureNode *fn;
  GtFeatureNodeIterator *fni;
  unsigned long nof_specified_ft_found;
  unsigned long next_inside;
  unsigned long next_not_inside;
  const char *specified_ft;
};

GtSeqposClassifier *gt_seqpos_classifier_new(const char *filename,
    const char *feature_type)
{
  GtSeqposClassifier *seqpos_classifier;
  seqpos_classifier = gt_malloc(sizeof (GtSeqposClassifier));
  seqpos_classifier->annotation_stream = gt_gff3_in_stream_new_sorted(filename);
  seqpos_classifier->fn = NULL;
  seqpos_classifier->fni = NULL;
  seqpos_classifier->gn = NULL;
  seqpos_classifier->nof_specified_ft_found = 0;
  seqpos_classifier->specified_ft = feature_type;
  return seqpos_classifier;
}

static int gt_seqpos_classifier_next_fn(GtSeqposClassifier *seqpos_classifier,
    GtError *err)
{
  int had_err = 0;
  gt_assert(seqpos_classifier != NULL);
  if (seqpos_classifier->fni != NULL)
  {
    gt_feature_node_iterator_delete(seqpos_classifier->fni);
    seqpos_classifier->fni = NULL;
  }
  while (true)
  {
    if (seqpos_classifier->gn != NULL)
    {
      gt_genome_node_delete(seqpos_classifier->gn);
    }
    had_err = gt_node_stream_next(seqpos_classifier->annotation_stream,
        &seqpos_classifier->gn, err);
    if (had_err != 0 || seqpos_classifier->gn == NULL)
    {
      seqpos_classifier->fn = NULL;
      seqpos_classifier->gn = NULL;
      return had_err;
    }
    else
    {
      if ((seqpos_classifier->fn =
            gt_feature_node_try_cast(seqpos_classifier->gn)) != NULL)
      {
        seqpos_classifier->fni =
          gt_feature_node_iterator_new(seqpos_classifier->fn);
        return had_err;
      }
    }
  }
}

static int gt_seqpos_classifier_next_specified_ft(
    GtSeqposClassifier *seqpos_classifier, GtRange *range,
    bool *end_of_annotation, GtError *err)
{
  int had_err = 0;
  GtFeatureNode *cfn;
  bool fni_exhausted = (seqpos_classifier->fni == NULL) ? true : false;
  gt_assert(seqpos_classifier != NULL);
  gt_assert(range != NULL);
  while (true)
  {
    if (fni_exhausted)
    {
      had_err = gt_seqpos_classifier_next_fn(seqpos_classifier, err);
      if (had_err != 0 || seqpos_classifier->fn == NULL)
      {
        *end_of_annotation = true;
        return had_err;
      }
      fni_exhausted = false;
    }
    gt_assert(seqpos_classifier->fni != NULL);
    cfn = gt_feature_node_iterator_next(seqpos_classifier->fni);
    if (cfn == NULL)
    {
      fni_exhausted = true;
    }
    else if (strcmp(gt_feature_node_get_type(cfn),
          seqpos_classifier->specified_ft) == 0)
    {
      seqpos_classifier->nof_specified_ft_found++;
      *range = gt_genome_node_get_range((GtGenomeNode*)cfn);
      gt_assert(range->start > 0);
      gt_assert(range->end > 0);
      range->start--;
      range->end--;
      *end_of_annotation = false;
      return had_err;
    }
  }
}

int gt_seqpos_classifier_position_is_inside_feature(
    GtSeqposClassifier *seqpos_classifier, unsigned long i, bool *inside,
    bool *end_of_annotation, GtError *err)
{
  int had_err = 0;
  GtRange next_specified_ft_range = { GT_UNDEF_ULONG, GT_UNDEF_ULONG };
  if (i == 0)
  {
    had_err = gt_seqpos_classifier_next_specified_ft(seqpos_classifier,
        &next_specified_ft_range, end_of_annotation, err);
    seqpos_classifier->next_inside = ULONG_MAX;
    seqpos_classifier->next_not_inside = ULONG_MAX;
    if (!had_err && !(*end_of_annotation))
    {
      seqpos_classifier->next_inside = next_specified_ft_range.start;
      seqpos_classifier->next_not_inside = next_specified_ft_range.end + 1UL;
      (*inside) = (seqpos_classifier->next_inside == 0);
    }
  }
  else
  {
    if (!(*inside))
    {
      gt_assert(i <= seqpos_classifier->next_inside);
      if (i == seqpos_classifier->next_inside)
      {
        (*inside) = true;
        /*printf("%lu..", i+1);*/
        gt_assert(seqpos_classifier->next_not_inside >
            seqpos_classifier->next_inside);
      }
    }
    else
    {
      if (i == seqpos_classifier->next_not_inside)
      {
        next_specified_ft_range.start = 0;
        while (!had_err)
        {
          had_err = gt_seqpos_classifier_next_specified_ft(seqpos_classifier,
              &next_specified_ft_range, end_of_annotation, err);
          if (!had_err)
          {
            if (*end_of_annotation)
            {
              (*inside) = false;
              seqpos_classifier->next_inside = ULONG_MAX;
              seqpos_classifier->next_not_inside = ULONG_MAX;
              break;
            }
            else
            {
              if (next_specified_ft_range.start <= i)
              {
                if (next_specified_ft_range.end >= i)
                {
                  (*inside) = true;
                  seqpos_classifier->next_not_inside =
                    next_specified_ft_range.end + 1UL;
                  break;
                }
              }
              else
              {
                /*printf("%lu\n", i);*/
                (*inside) = false;
                gt_assert(next_specified_ft_range.end > i);
                seqpos_classifier->next_inside = next_specified_ft_range.start;
                seqpos_classifier->next_not_inside =
                  next_specified_ft_range.end + 1UL;
                break;
              }
            }
          }
        }
      }
    }
  }
  return had_err;
}

unsigned long gt_seqpos_classifier_nof_features_found(
    GtSeqposClassifier *seqpos_classifier)
{
  gt_assert(seqpos_classifier);
  return seqpos_classifier->nof_specified_ft_found;
}

void gt_seqpos_classifier_delete(GtSeqposClassifier *seqpos_classifier)
{
  if (seqpos_classifier != NULL)
  {
    if (seqpos_classifier->fni != NULL)
    {
      gt_feature_node_iterator_delete(seqpos_classifier->fni);
    }
    if (seqpos_classifier->gn != NULL)
    {
      gt_genome_node_delete(seqpos_classifier->gn);
    }
    gt_node_stream_delete(seqpos_classifier->annotation_stream);
    gt_free(seqpos_classifier);
  }
}
