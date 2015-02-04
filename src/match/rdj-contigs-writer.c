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
#include "core/fasta.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/xansi_api.h"
#include "extended/assembly_stats_calculator.h"
#include "match/rdj-contig-info.h"
#include "match/rdj-contigs-writer.h"

struct GtContigsWriter
{
  const GtEncseq *reads;
  GtFile *outfp;
  bool show_paths, calculate_astat;
  GtAssemblyStatsCalculator *asc;
  GtEncseqReader *esr;
  GtArraychar contig;
  GtStr *contig_desc, *path_desc;
  GtUword contignum, lastseqnum, nofseqs, rlen;
  double arrival_rate;
  GtContigDepthInfo depthinfo;
  unsigned char *rcn;
  FILE *depthinfo_fp;
};

#define GT_CONTIGS_WRITER_CONTIG_INC 16384UL

GtContigsWriter *gt_contigs_writer_new(const GtEncseq *reads, GtFile *outfp)
{
  GtContigsWriter *contigs_writer;

  gt_assert(reads != NULL);

  contigs_writer = gt_malloc(sizeof (GtContigsWriter));
  contigs_writer->reads = reads;
  contigs_writer->outfp = outfp;
  contigs_writer->show_paths = false;
  contigs_writer->nofseqs = gt_encseq_num_of_sequences(contigs_writer->reads);
  contigs_writer->contignum = 0;
  contigs_writer->arrival_rate = (double)0;
  contigs_writer->depthinfo.firstread_copynum = (float)1;
  contigs_writer->depthinfo.internal_copynum = (float)0;
  contigs_writer->depthinfo.lastread_copynum = (float)1;
  contigs_writer->depthinfo.depth = 0;
  contigs_writer->depthinfo.length = 0;
  contigs_writer->depthinfo.astat = (float)0;
  contigs_writer->rlen = 0;
  contigs_writer->calculate_astat = false;
  GT_INITARRAY(&contigs_writer->contig, char);
  contigs_writer->contig_desc = gt_str_new();
  contigs_writer->path_desc = gt_str_new();
  contigs_writer->asc = gt_assembly_stats_calculator_new();
  contigs_writer->esr = gt_encseq_create_reader_with_readmode(reads,
      GT_READMODE_FORWARD, 0);
  contigs_writer->rcn = NULL;
  contigs_writer->depthinfo_fp = NULL;
  return contigs_writer;
}

void gt_contigs_writer_enable_complete_path_output(
    GtContigsWriter *contigs_writer)
{
  contigs_writer->show_paths = true;
}

void gt_contigs_writer_enable_astat_calculation(GtContigsWriter *contigs_writer,
    double coverage, unsigned char *rcn, FILE *depthinfo_fp)
{
  contigs_writer->calculate_astat = true;
  if (gt_encseq_accesstype_get(contigs_writer->reads) ==
      GT_ACCESS_TYPE_EQUALLENGTH)
    contigs_writer->rlen = gt_encseq_equallength(contigs_writer->reads);
  else
    contigs_writer->rlen = gt_encseq_total_length(contigs_writer->reads) /
      gt_encseq_num_of_sequences(contigs_writer->reads);
  contigs_writer->arrival_rate = coverage / (double)contigs_writer->rlen;
  gt_log_log("arrival_rate = %.4f", contigs_writer->arrival_rate);
  contigs_writer->rcn = rcn;
  contigs_writer->depthinfo_fp = depthinfo_fp;
}

void gt_contigs_writer_delete(GtContigsWriter *contigs_writer)
{
  if (contigs_writer == NULL)
    return;
  GT_FREEARRAY(&contigs_writer->contig, char);
  gt_str_delete(contigs_writer->contig_desc);
  gt_str_delete(contigs_writer->path_desc);
  gt_assembly_stats_calculator_delete(contigs_writer->asc);
  gt_encseq_reader_delete(contigs_writer->esr);
  gt_free(contigs_writer);
}

static void gt_contigs_writer_reset_buffers(GtContigsWriter *contigs_writer)
{
  gt_assert(contigs_writer != NULL);
  contigs_writer->contig.nextfreechar = 0;
  gt_str_reset(contigs_writer->contig_desc);
  gt_str_reset(contigs_writer->path_desc);
  contigs_writer->depthinfo.firstread_copynum = (float)0;
  contigs_writer->depthinfo.internal_copynum = (float)0;
  contigs_writer->depthinfo.lastread_copynum = (float)0;
  contigs_writer->depthinfo.depth = 0;
  contigs_writer->depthinfo.length = 0;
  contigs_writer->depthinfo.astat = (float)0;
}

void gt_contigs_writer_abort(GtContigsWriter *contigs_writer)
{
  gt_assert(contigs_writer != NULL);
  gt_contigs_writer_reset_buffers(contigs_writer);
}

GT_UNUSED
static void gt_contigs_writer_append_chars(GtContigsWriter *contigs_writer,
    GtUword nofchars)
{
  GtUword i;
  char *c;

  gt_assert(contigs_writer != NULL);
  for (i = 0; i < nofchars; i++)
  {
    GT_GETNEXTFREEINARRAY(c, &contigs_writer->contig, char,
        GT_CONTIGS_WRITER_CONTIG_INC);
    *c = gt_encseq_reader_next_decoded_char(contigs_writer->esr);
  }
}

#define GT_CONTIGS_WRITER_IS_DIRECT(SEQNUM, NOFSEQS)\
        ((SEQNUM) < ((NOFSEQS) >> 1))
#define GT_CONTIGS_WRITER_READNUM(SEQNUM, NOFSEQS)\
        (GT_CONTIGS_WRITER_IS_DIRECT(SEQNUM, NOFSEQS)\
            ? (SEQNUM)\
            : (NOFSEQS) - (SEQNUM) - 1)
#define GT_CONTIGS_WRITER_LETTER(SEQNUM, NOFSEQS)\
        (GT_CONTIGS_WRITER_IS_DIRECT(SEQNUM, NOFSEQS) ? 'E' : 'B')

void gt_contigs_writer_write(GtContigsWriter *contigs_writer)
{
  gt_assert(contigs_writer != NULL);

  if (contigs_writer->contig.nextfreechar > 0)
  {
    gt_assembly_stats_calculator_add(contigs_writer->asc,
        contigs_writer->contig.nextfreechar);
    /* build description */
    gt_str_append_cstr(contigs_writer->contig_desc, "contig_");
    gt_str_append_uword(contigs_writer->contig_desc,
        contigs_writer->contignum);
    contigs_writer->depthinfo.length = contigs_writer->contig.nextfreechar;
    gt_str_append_cstr(contigs_writer->contig_desc, " length=");
    gt_str_append_uword(contigs_writer->contig_desc,
        contigs_writer->depthinfo.length);
    gt_str_append_cstr(contigs_writer->contig_desc, " depth=");
    gt_str_append_uword(contigs_writer->contig_desc,
        contigs_writer->depthinfo.depth);
    if (contigs_writer->rcn != NULL)
    {
      contigs_writer->depthinfo.internal_copynum -= (float)contigs_writer->rcn[
        GT_CONTIGS_WRITER_READNUM(contigs_writer->lastseqnum,
            contigs_writer->nofseqs)];
      contigs_writer->depthinfo.lastread_copynum = (float)contigs_writer->rcn[
        GT_CONTIGS_WRITER_READNUM(contigs_writer->lastseqnum,
            contigs_writer->nofseqs)];
      gt_str_append_cstr(contigs_writer->contig_desc, " k=");
      gt_str_append_double(contigs_writer->contig_desc,
          (double)contigs_writer->depthinfo.internal_copynum, 2);
    }
    if (contigs_writer->calculate_astat)
    {
      float astat;
      float k = (float)0;
      if (contigs_writer->rcn != NULL)
        k = contigs_writer->depthinfo.internal_copynum;
      else
      {
        if (contigs_writer->depthinfo.depth > 2UL)
          k = (float)(contigs_writer->depthinfo.depth - 2UL);
      }
      astat = - (k * (float)log((double)2));
      if (contigs_writer->contig.nextfreechar + 1UL > contigs_writer->rlen)
      {
        astat += contigs_writer->arrival_rate * (
            contigs_writer->contig.nextfreechar + 1UL - contigs_writer->rlen);
      }
      gt_str_append_cstr(contigs_writer->contig_desc, " astat=");
      gt_str_append_double(contigs_writer->contig_desc, (double)astat, 6);
      contigs_writer->depthinfo.astat = astat;
    }
    gt_str_append_cstr(contigs_writer->contig_desc, " ");
    if (!contigs_writer->show_paths && contigs_writer->depthinfo.depth > 1UL)
    {
      /* show path summary only */
      gt_str_append_cstr(contigs_writer->path_desc,
          contigs_writer->depthinfo.depth > 2UL ? "-->...-->" :
          "-->");
      gt_str_append_uword(contigs_writer->path_desc,
          GT_CONTIGS_WRITER_READNUM(contigs_writer->lastseqnum,
              contigs_writer->nofseqs));
      gt_str_append_char(contigs_writer->path_desc,
          GT_CONTIGS_WRITER_LETTER(contigs_writer->lastseqnum,
              contigs_writer->nofseqs));
    }
    gt_str_append_str(contigs_writer->contig_desc, contigs_writer->path_desc);

    gt_fasta_show_entry(gt_str_get(contigs_writer->contig_desc),
        contigs_writer->contig.spacechar, contigs_writer->contig.nextfreechar,
        60UL, contigs_writer->outfp);
    if (contigs_writer->depthinfo_fp != NULL)
    {
      gt_xfwrite(&contigs_writer->depthinfo,
          sizeof (contigs_writer->depthinfo), (size_t)1,
          contigs_writer->depthinfo_fp);
    }
    contigs_writer->contignum++;
    gt_contigs_writer_reset_buffers(contigs_writer);
  }
}

void gt_contigs_writer_append(GtContigsWriter *contigs_writer,
    GtUword seqnum, GtUword nofchars)
{
  GtUword pos, i;
  gt_assert(contigs_writer != NULL);
  pos = gt_encseq_seqstartpos(contigs_writer->reads, seqnum) +
      gt_encseq_seqlength(contigs_writer->reads, seqnum) - nofchars;
  for (i = 0; i < nofchars; i++, pos++)
  {
    char *c;
    GT_GETNEXTFREEINARRAY(c, &contigs_writer->contig, char,
        GT_CONTIGS_WRITER_CONTIG_INC);
    *c = gt_encseq_get_encoded_char_nospecial(contigs_writer->reads, pos,
        GT_READMODE_FORWARD)["acgt"];
  }
  contigs_writer->depthinfo.depth++;
  if (contigs_writer->rcn != NULL)
  {
    contigs_writer->depthinfo.internal_copynum += contigs_writer->rcn[
      GT_CONTIGS_WRITER_READNUM(seqnum, contigs_writer->nofseqs)];
  }
  if (contigs_writer->show_paths)
  {
    gt_str_append_cstr(contigs_writer->path_desc, "-(");
    gt_str_append_uword(contigs_writer->path_desc, nofchars);
    gt_str_append_cstr(contigs_writer->path_desc, ")->");
    gt_str_append_uword(contigs_writer->path_desc,
        GT_CONTIGS_WRITER_READNUM(seqnum, contigs_writer->nofseqs));
    gt_str_append_char(contigs_writer->path_desc,
        GT_CONTIGS_WRITER_LETTER(seqnum, contigs_writer->nofseqs));
  }
  else
  {
    contigs_writer->lastseqnum = seqnum;
  }
}

void gt_contigs_writer_start(GtContigsWriter *contigs_writer,
    GtUword seqnum)
{
  GtUword pos, nofchars, i;
  gt_assert(contigs_writer != NULL);
  pos = gt_encseq_seqstartpos(contigs_writer->reads, seqnum);
  nofchars = gt_encseq_seqlength(contigs_writer->reads, seqnum);
  for (i = 0; i < nofchars; i++, pos++)
  {
    char *c;
    GT_GETNEXTFREEINARRAY(c, &contigs_writer->contig, char,
        GT_CONTIGS_WRITER_CONTIG_INC);
    *c = gt_encseq_get_encoded_char_nospecial(contigs_writer->reads, pos,
        GT_READMODE_FORWARD)["acgt"];
  }
  contigs_writer->depthinfo.depth++;
  gt_str_append_uword(contigs_writer->path_desc,
      GT_CONTIGS_WRITER_READNUM(seqnum, contigs_writer->nofseqs));
  gt_str_append_char(contigs_writer->path_desc,
      GT_CONTIGS_WRITER_LETTER(seqnum, contigs_writer->nofseqs));
  if (contigs_writer->rcn != NULL)
  {
    contigs_writer->depthinfo.firstread_copynum = (float)contigs_writer->rcn[
      GT_CONTIGS_WRITER_READNUM(seqnum, contigs_writer->nofseqs)];
  }
}

void gt_contigs_writer_show_stats(GtContigsWriter *contigs_writer,
    GtLogger *logger)
{
  gt_assert(contigs_writer != NULL);
  if (contigs_writer->contignum > 0)
    gt_assembly_stats_calculator_show(contigs_writer->asc, logger);
  else
    gt_logger_log(logger, "no contigs respect the given cutoff parameters");
}
