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

#include "core/hashmap.h"
#include "core/log.h"
#include "core/cstr_api.h"
#include "core/ma.h"
#include "extended/samfile_encseq_mapping.h"

struct GtSamfileEncseqMapping
{
  const GtEncseq *encseq;
  unsigned long *samfile2encseq;
  unsigned long nof_sequences;
};

GtSamfileEncseqMapping *gt_samfile_encseq_mapping_new(
    /*const*/ GtSamfileIterator *samfile_iterator, const GtEncseq *encseq,
    GtError *err)
{
  GtSamfileEncseqMapping *samfile_encseq_mapping;
  bool haserr = false;
  gt_assert(samfile_iterator != NULL);
  gt_assert(encseq != NULL);
  gt_error_check(err);
  samfile_encseq_mapping = gt_malloc(sizeof (GtSamfileEncseqMapping));
  samfile_encseq_mapping->nof_sequences = gt_encseq_num_of_sequences(encseq);
  samfile_encseq_mapping->samfile2encseq = NULL;
  samfile_encseq_mapping->encseq = encseq;
  if ((unsigned long)gt_samfile_iterator_number_of_references(samfile_iterator)
      != samfile_encseq_mapping->nof_sequences)
  {
    gt_error_set(err, "The number of sequences in the encoded sequence "
        "is not the same as the number of references in the SAM file.");
    haserr = true;
  }
  else
  {
    GtHashmap *seqid2seqnum;
    unsigned long i, cnum, dlen, *stored, *seqnum;
    int32_t refnum;
    const char *d;
    char *seqid = NULL;
    seqid2seqnum = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
    for (i = 0; i < samfile_encseq_mapping->nof_sequences; i++)
    {

      d = gt_encseq_description(encseq, &dlen, i);
      seqid = gt_cstr_dup_nt(d, dlen);
      for (cnum = 0; cnum < dlen; cnum++)
      {
        if (seqid[cnum] == ' ')
        {
          seqid[cnum] = '\0';
          break;
        }
      }
      if ((stored = gt_hashmap_get(seqid2seqnum, seqid)) != NULL)
      {
        gt_error_set(err, "The encseq identifier '%s' is not unique!", seqid);
        haserr = true;
        break;
      }
      else
      {
        seqnum = gt_malloc(sizeof (unsigned long));
        *seqnum = i;
        gt_hashmap_add(seqid2seqnum, seqid, seqnum);
      }
    }
    if (!haserr)
    {
      samfile_encseq_mapping->samfile2encseq = gt_malloc(
          sizeof (*samfile_encseq_mapping->samfile2encseq) *
          samfile_encseq_mapping->nof_sequences);
      for (refnum = 0; refnum < (int32_t)samfile_encseq_mapping->nof_sequences;
          refnum++)
      {
        d = gt_samfile_iterator_reference_name(samfile_iterator, refnum);
        if ((seqnum = gt_hashmap_get(seqid2seqnum, d)) == NULL)
        {
          gt_error_set(err,
              "The SAM identifier '%s' was not found in the encseq!", d);
          haserr = true;
          break;
        }
        samfile_encseq_mapping->samfile2encseq[refnum] = *seqnum;
      }
    }
    gt_hashmap_delete(seqid2seqnum);
  }
  if (haserr)
  {
    gt_samfile_encseq_mapping_delete(samfile_encseq_mapping);
    return NULL;
  }
  else
    return samfile_encseq_mapping;
}

static unsigned long gt_samfile_encseq_mapping_seqnum(
    GtSamfileEncseqMapping *samfile_encseq_mapping,
    int32_t reference_num)
{
  gt_assert(samfile_encseq_mapping != NULL);
  gt_assert(reference_num < (int32_t)samfile_encseq_mapping->nof_sequences);
  return samfile_encseq_mapping->samfile2encseq[reference_num];
}

unsigned long gt_samfile_encseq_mapping_seqpos(
    GtSamfileEncseqMapping *samfile_encseq_mapping, int32_t reference_num,
    unsigned long reference_seqpos)
{
  unsigned long seqnum, seqstartpos;
  gt_assert(samfile_encseq_mapping != NULL);
  seqnum = gt_samfile_encseq_mapping_seqnum(samfile_encseq_mapping,
      reference_num);
  seqstartpos = gt_encseq_seqstartpos(samfile_encseq_mapping->encseq, seqnum);
  return seqstartpos + reference_seqpos;
}

void gt_samfile_encseq_mapping_delete(
    GtSamfileEncseqMapping *samfile_encseq_mapping)
{
  gt_free(samfile_encseq_mapping->samfile2encseq);
  gt_free(samfile_encseq_mapping);
}
