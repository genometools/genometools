/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/arraydef.h"
#include "core/fa.h"
#include "core/fasta_separator.h"
#include "core/log_api.h"
#include "core/xansi_api.h"
#include "core/str.h"
#include "extended/assembly_stats_calculator.h"
#include "match/rdj-contigs-writer.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-contigpaths.h"

static GtFile* gt_contigpaths_get_file(const char *indexname,
    const char *suffix, const char* mode, GtError *err)
{
  GtStr *path;
  GtFile *file;

  path = gt_str_new_cstr(indexname);
  gt_str_append_cstr(path, suffix);
  file = gt_file_new(gt_str_get(path), mode, err);
  gt_str_delete(path);
  return file;
}

#define GT_CONTIGPATHS_BUFFERSIZE ((size_t)1 << 16)

int gt_contigpaths_to_fasta(const char *indexname,
    const char *contigpaths_suffix, const char *fasta_suffix,
    const GtEncseq *encseq, GtUword min_contig_length, bool showpaths,
    bool astat, double coverage, bool load_copynum, size_t buffersize,
    GtLogger *logger, GtError *err)
{
  GtUword nofchars, seqnum, contig_length = 0;
  GtFile *infp = NULL, *outfp = NULL;
  FILE *depthinfo_fp = NULL;
  int nvalues, i;
  GtContigsWriter *cw = NULL;
  GtContigpathElem *buffer = NULL;
  unsigned char *rcn = NULL;
  int had_err = 0;

  if (buffersize == 0)
  {
    buffersize = GT_CONTIGPATHS_BUFFERSIZE * sizeof (GtContigpathElem);
    gt_log_log("buffersize = default ("GT_ZU" bytes)", buffersize);
  }
  else
  {
    size_t bsmod = buffersize % sizeof (GtContigpathElem);
    buffersize -= bsmod;
    gt_log_log("buffersize = "GT_ZU" bytes", buffersize);
  }
  gt_assert(buffersize > 0);
  gt_assert(buffersize % sizeof (GtContigpathElem) == 0);
  buffer = gt_malloc(buffersize);

  gt_assert(encseq != NULL);
  gt_error_check(err);
  gt_assert(buffer != NULL);

  infp = gt_contigpaths_get_file(indexname, contigpaths_suffix, "r", err);
  if (infp == NULL)
    had_err = -1;
  if (!had_err)
  {
    outfp = gt_contigpaths_get_file(indexname, fasta_suffix, "w", err);
    if (outfp == NULL)
    {
      gt_file_delete(infp);
      had_err = -1;
    }
  }
  if (!had_err && load_copynum)
  {
    GtFile *rcnfp = NULL;
    size_t nbytes;
    rcnfp = gt_contigpaths_get_file(indexname,
        GT_READJOINER_SUFFIX_READSCOPYNUM, "r", err);
    if (rcnfp == NULL)
    {
      had_err = -1;
    }
    else
    {
      gt_log_log("load reads copy number from file");
      nbytes = sizeof (*rcn) * gt_encseq_num_of_sequences(encseq);
      if (gt_encseq_is_mirrored(encseq))
        nbytes = GT_DIV2(nbytes);
      rcn = gt_malloc(nbytes);
      (void)gt_file_xread(rcnfp, rcn, nbytes);
      gt_file_delete(rcnfp);
    }
  }
  if (!had_err)
  {
    cw = gt_contigs_writer_new(encseq, outfp);
    if (showpaths)
      gt_contigs_writer_enable_complete_path_output(cw);
    if (astat)
    {
      depthinfo_fp = gt_fa_fopen_with_suffix(indexname,
          GT_READJOINER_SUFFIX_DEPTHINFO, "w", err);
      if (depthinfo_fp == NULL)
      {
        had_err = -1;
      }
      if (!had_err)
      {
        gt_contigs_writer_enable_astat_calculation(cw, coverage, rcn,
            depthinfo_fp);
      }
    }
  }
  if (!had_err)
  {
    while ((nvalues = gt_file_xread(infp, buffer, buffersize)) > 0)
    {
      gt_assert((size_t)nvalues <= buffersize);
      nvalues /= (sizeof (GtContigpathElem) << 1);
      for (i = 0; i < nvalues; i++)
      {
        nofchars = (GtUword)buffer[(i << 1)];
        seqnum = (GtUword)buffer[(i << 1) + 1];
        if (nofchars == 0)
        {
          /* end of contig */
          if (contig_length >= min_contig_length)
            gt_contigs_writer_write(cw);
          else
            gt_contigs_writer_abort(cw);
          gt_contigs_writer_start(cw, seqnum);
          contig_length = gt_encseq_seqlength(encseq, seqnum);
        }
        else
        {
          contig_length += nofchars;
          gt_contigs_writer_append(cw, seqnum, nofchars);
        }
      }
    }

    if (contig_length >= min_contig_length)
      gt_contigs_writer_write(cw);
    else
      gt_contigs_writer_abort(cw);
    gt_contigs_writer_show_stats(cw, logger);
  }
  if (depthinfo_fp != NULL)
    gt_fa_fclose(depthinfo_fp);
  gt_contigs_writer_delete(cw);
  gt_file_delete(infp);
  gt_file_delete(outfp);
  gt_free(buffer);
  gt_free(rcn);
  return 0;
}
