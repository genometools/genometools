/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007/2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include "core/encseq_metadata.h"
#include "core/encseq.h"
#include "core/encseq_rep.h"
#include "core/fa.h"
#include "core/ma.h"

struct GtEncseqMetadata
{
  GtEncseqAccessType sat;
  unsigned long totallength;
  unsigned long numofdbsequences,
                numofdbfiles,
                lengthofdbfilenames;
  GtSpecialcharinfo specialcharinfo;
};

#define NEXTFREAD(VAL)\
        if (!had_err)\
        {\
          size_t ret;\
          ret = fread(&(VAL), sizeof (VAL), (size_t) 1, fp);\
          if (ferror(fp))\
          {\
            gt_error_set(err,"error when trying to read %s: %s",\
                              #VAL, strerror(errno));\
            had_err = true;\
          }\
        }

static int readfirstvaluesfromfile(GtEncseqMetadata *emd,
                                   const char *indexname, GtError *err)
{
  FILE *fp;
  bool had_err = false;
  unsigned long cc;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname, GT_ENCSEQFILESUFFIX, "rb", err);
  if (fp == NULL)
  {
    had_err = true;
  }
  NEXTFREAD(cc);
  if (!had_err)
  {
    if (cc >= (unsigned long) GT_ACCESS_TYPE_UNDEFINED)
    {
      gt_error_set(err, "illegal type %lu in \"%s%s\"", cc,
                   indexname, GT_ENCSEQFILESUFFIX);
      had_err = true;
    }
  }
  emd->sat = (GtEncseqAccessType) cc;
  NEXTFREAD(emd->totallength);
  NEXTFREAD(emd->numofdbsequences);
  NEXTFREAD(emd->numofdbfiles);
  NEXTFREAD(emd->lengthofdbfilenames);
  NEXTFREAD(emd->specialcharinfo);
  gt_fa_xfclose(fp);
  return had_err ? -1 : 0;
}

GtEncseqMetadata *gt_encseq_metadata_new(const char *indexname, GtError *err)
{
  int had_err = 0;
  GtEncseqMetadata *encseq_metadata;
  gt_assert(indexname);
  encseq_metadata = gt_malloc(sizeof (GtEncseqMetadata));
  had_err = readfirstvaluesfromfile(encseq_metadata, indexname, err);
  if (had_err) {
    gt_assert(gt_error_is_set(err));
    gt_free(encseq_metadata);
    encseq_metadata = NULL;
  }
  return encseq_metadata;
}

unsigned long gt_encseq_metadata_total_length(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->totallength;
}

unsigned long gt_encseq_metadata_num_of_sequences(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->numofdbsequences;
}

unsigned long gt_encseq_metadata_num_of_files(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->numofdbfiles;
}

unsigned long gt_encseq_metadata_length_of_filenames(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->lengthofdbfilenames;
}

GtEncseqAccessType gt_encseq_metadata_accesstype(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->sat;
}

GtSpecialcharinfo gt_encseq_metadata_specialcharinfo(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->specialcharinfo;
}

void gt_encseq_metadata_delete(GtEncseqMetadata *emd)
{
  gt_free(emd);
}
