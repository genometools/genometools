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
#include "core/intbits.h"
#include "core/ma.h"

struct GtEncseqMetadata
{
  GtEncseqAccessType sat;
  GtUchar is64bit;
  unsigned long version,
                totallength,
                numofdbsequences,
                numofdbfiles,
                lengthofdbfilenames,
                lengthofalphadef,
                minseqlen,
                maxseqlen;
  bool customalphabet;
  GtAlphabet *alpha;
  GtSpecialcharinfo specialcharinfo;
};

#define NEXTFREADWSIZE(VAL, SIZE)\
        if (!had_err)\
        {\
          GT_UNUSED size_t ret;\
          ret = fread(&(VAL), sizeof (VAL), (size_t) SIZE, fp);\
          if (ferror(fp))\
          {\
            gt_error_set(err,"error when trying to read %s: %s",\
                              #VAL, strerror(errno));\
            had_err = true;\
          }\
          if (!had_err) {\
            byteoffset += sizeof (VAL);\
            if (byteoffset % (unsigned long) GT_WORDSIZE_INBYTES > 0)\
            {\
              char buffer[GT_WORDSIZE_INBYTES];\
              size_t padunits\
                = GT_WORDSIZE_INBYTES - (byteoffset % GT_WORDSIZE_INBYTES);\
              byteoffset += (unsigned long) padunits;\
              ret = fread(buffer, (size_t) 1, (size_t) padunits, fp);\
            }\
            if (ferror(fp))\
            {\
              gt_error_set(err,"error when trying to read %s: %s",\
                                #VAL, strerror(errno));\
              had_err = true;\
            }\
          }\
        }

#define NEXTFREAD(VAL)\
        NEXTFREADWSIZE(VAL,1)

static int readfirstvaluesfromfile(GtEncseqMetadata *emd,
                                   const char *indexname, GtError *err)
{
  FILE *fp;
  bool had_err = false;
  unsigned long cc, byteoffset = 0, alphatype;
  char *alphadef;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname, GT_ENCSEQFILESUFFIX, "rb", err);
  if (fp == NULL)
  {
    had_err = true;
  }
  NEXTFREAD(emd->is64bit);
  if (!had_err)
  {
    if ((int) emd->is64bit > 1)
    {
      gt_error_set(err, "illegal platform code %u in \"%s%s\"", emd->is64bit,
                   indexname, GT_ENCSEQFILESUFFIX);
      had_err = true;
    }
    if (!had_err && ((emd->is64bit && sizeof (unsigned long) != (size_t) 8)
          || (!emd->is64bit && sizeof (unsigned long) == (size_t) 8)))
    {
      gt_error_set(err, "trying to load 64-bit index \"%s%s\" on a 32-bit "
                        "system or vice versa -- please use correct index "
                        "for this platform", indexname, GT_ENCSEQFILESUFFIX);
      had_err = true;
    }
  }
  NEXTFREAD(emd->version);
  if (!had_err)
  {
    if (emd->version < GT_ENCSEQ_VERSION)    {
      gt_error_set(err, "index \"%s%s\" is format version %lu, current is "
                        "%lu -- please re-encode",
                        indexname, GT_ENCSEQFILESUFFIX,
                        emd->version, GT_ENCSEQ_VERSION);
      had_err = true;
    }
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
  if (!had_err) {
    emd->sat = (GtEncseqAccessType) cc;
    NEXTFREAD(emd->totallength);
    NEXTFREAD(emd->numofdbsequences);
    NEXTFREAD(emd->numofdbfiles);
    NEXTFREAD(emd->lengthofdbfilenames);
    NEXTFREAD(emd->specialcharinfo);
    NEXTFREAD(emd->minseqlen);
    NEXTFREAD(emd->maxseqlen);
  }
  NEXTFREAD(alphatype);
  if (!had_err) {
    if (alphatype > 2UL) {
      gt_error_set(err, "illegal alphabet type %lu in \"%s%s\"", alphatype,
                   indexname, GT_ENCSEQFILESUFFIX);
      had_err = true;
    }
  }
  if (!had_err) {
    NEXTFREAD(emd->lengthofalphadef);
    switch (alphatype) {
      case 0:
        emd->alpha = gt_alphabet_new_dna();
        break;
      case 1:
        emd->alpha = gt_alphabet_new_protein();
        break;
      case 2:
        gt_assert(emd->lengthofalphadef > 0);
        emd->customalphabet = true;
        alphadef = gt_malloc(sizeof (char) * emd->lengthofalphadef);
        NEXTFREADWSIZE(*(alphadef), emd->lengthofalphadef);
        emd->alpha = gt_alphabet_new_from_string(alphadef,
                                                 emd->lengthofalphadef,
                                                 err);
        if (!emd->alpha) {
          had_err = true;
        }
        gt_free(alphadef);
        break;
    }
    gt_assert(emd->alpha != NULL);
  }
  gt_fa_xfclose(fp);
  return had_err ? -1 : 0;
}

GtEncseqMetadata* gt_encseq_metadata_new(const char *indexname, GtError *err)
{
  int had_err = 0;
  GtEncseqMetadata *encseq_metadata;
  gt_assert(indexname);
  encseq_metadata = gt_malloc(sizeof (GtEncseqMetadata));
  encseq_metadata->alpha = NULL;
  encseq_metadata->customalphabet = false;
  had_err = readfirstvaluesfromfile(encseq_metadata, indexname, err);
  if (had_err) {
    gt_assert(gt_error_is_set(err));
    gt_free(encseq_metadata);
    encseq_metadata = NULL;
  }
  return encseq_metadata;
}

GtAlphabet* gt_encseq_metadata_alphabet(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->alpha;
}

unsigned long gt_encseq_metadata_total_length(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->totallength;
}

unsigned long gt_encseq_metadata_version(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->version;
}

bool gt_encseq_metadata_is64bit(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return ((int) emd->is64bit == 1);
}

unsigned long gt_encseq_metadata_max_seq_length(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->maxseqlen;
}

unsigned long gt_encseq_metadata_min_seq_length(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->minseqlen;
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

unsigned long gt_encseq_metadata_length_of_alphadef(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->lengthofdbfilenames;
}

bool gt_encseq_metadata_has_custom_alphabet(GtEncseqMetadata *emd)
{
  gt_assert(emd != NULL);
  return emd->customalphabet;
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
  if (emd == NULL) return;
  if (emd->alpha != NULL)
    gt_alphabet_delete(emd->alpha);
  gt_free(emd);
}
