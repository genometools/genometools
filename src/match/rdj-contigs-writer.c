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
#include "core/ma.h"
#include "core/str.h"
#include "extended/assembly_stats_calculator.h"
#include "match/rdj-contigs-writer.h"

struct GtContigsWriter
{
  const GtEncseq *reads;
  GtFile *outfp;
  bool showpaths;
  GtAssemblyStatsCalculator *asc;
  GtEncseqReader *esr;
  GtArraychar contig;
  GtStr *contig_desc, *path_desc;
  unsigned long contignum, contigdepth, lastseqnum, nofseqs;
};

#define GT_CONTIGS_WRITER_CONTIG_INC 16384UL

GtContigsWriter *gt_contigs_writer_new(const GtEncseq *reads, GtFile *outfp,
    bool showpaths)
{
  GtContigsWriter *contigs_writer;

  gt_assert(reads != NULL);

  contigs_writer = gt_malloc(sizeof (GtContigsWriter));
  contigs_writer->reads = reads;
  contigs_writer->outfp = outfp;
  contigs_writer->showpaths = showpaths;
  contigs_writer->nofseqs = gt_encseq_num_of_sequences(contigs_writer->reads);
  contigs_writer->contignum = 0;
  contigs_writer->contigdepth = 0;
  GT_INITARRAY(&contigs_writer->contig, char);
  contigs_writer->contig_desc = gt_str_new();
  contigs_writer->path_desc = gt_str_new();
  contigs_writer->asc = gt_assembly_stats_calculator_new();
  contigs_writer->esr = gt_encseq_create_reader_with_readmode(reads,
      GT_READMODE_FORWARD, 0);
  return contigs_writer;
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
  contigs_writer->contigdepth = 0;
}

void gt_contigs_writer_abort(GtContigsWriter *contigs_writer)
{
  gt_assert(contigs_writer != NULL);
  gt_contigs_writer_reset_buffers(contigs_writer);
}

GT_UNUSED
static void gt_contigs_writer_append_chars(GtContigsWriter *contigs_writer,
    unsigned long nofchars)
{
  unsigned long i;
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
    gt_str_append_ulong(contigs_writer->contig_desc,
        contigs_writer->contignum);
    gt_str_append_cstr(contigs_writer->contig_desc, " length=");
    gt_str_append_ulong(contigs_writer->contig_desc,
        contigs_writer->contig.nextfreechar);
    gt_str_append_cstr(contigs_writer->contig_desc, " depth=");
    gt_str_append_ulong(contigs_writer->contig_desc,
        contigs_writer->contigdepth);
    gt_str_append_cstr(contigs_writer->contig_desc, " ");
    if (!contigs_writer->showpaths)
    {
      /* show path summary only */
      gt_str_append_cstr(contigs_writer->path_desc, "-->...-->");
      gt_str_append_ulong(contigs_writer->path_desc,
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
    contigs_writer->contignum++;
    gt_contigs_writer_reset_buffers(contigs_writer);
  }
}

void gt_contigs_writer_append(GtContigsWriter *contigs_writer,
    unsigned long seqnum, unsigned long nofchars)
{
  unsigned long pos, i;
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
  contigs_writer->contigdepth++;
  if (contigs_writer->showpaths)
  {
    gt_str_append_cstr(contigs_writer->path_desc, "-(");
    gt_str_append_ulong(contigs_writer->path_desc, nofchars);
    gt_str_append_cstr(contigs_writer->path_desc, ")->");
    gt_str_append_ulong(contigs_writer->path_desc,
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
    unsigned long seqnum)
{
  unsigned long pos, nofchars, i;
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
  contigs_writer->contigdepth++;
  gt_str_append_ulong(contigs_writer->path_desc,
      GT_CONTIGS_WRITER_READNUM(seqnum, contigs_writer->nofseqs));
  gt_str_append_char(contigs_writer->path_desc,
      GT_CONTIGS_WRITER_LETTER(seqnum, contigs_writer->nofseqs));
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
