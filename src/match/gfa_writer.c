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

#include "core/unused_api.h"
#include "core/log_api.h"
#include "match/gfa_writer.h"

struct GtGfaWriter
{
  GtFile *file;
  GtEncseq *encseq;
  GtGfaVersion version;
};

GtGfaWriter* gt_gfa_writer_new(GtFile *file,
                               const GtEncseq *encseq,
                               GtGfaVersion version)
{
  GtGfaWriter *gw;
  gw = gt_malloc(sizeof (*gw));
  gw->file = gt_file_ref(file);
  gw->encseq = gt_encseq_ref((GtEncseq*)encseq);
  gw->version = version;
  return gw;
}

static inline const char *gt_gfa_writer_versionstr(GtGfaVersion version)
{
  switch (version) {
    case GT_GFA_VERSION_1_0:
      return "1.0";
    case GT_GFA_VERSION_2_0:
      return "2.0";
    default:
      gt_assert(false);
  }
  return "UNKNOWN";
}

int gt_gfa_writer_show_header(GtGfaWriter *gw,
    GtUword minlen, const char *inputfilename, bool has_containments,
    bool has_transitives, GT_UNUSED GtError *err)
{
  gt_assert(gw != NULL);
  gt_file_xprintf(gw->file,
      "H\tVN:Z:%s\n"
      "H\tpn:Z:readjoiner\n"
      "H\tol:i:"GT_WU"\n"
      "H\tin:Z:%s\n"
      "H\tcn:i:%c\n"
      "H\tte:i:%c\n",
      gt_gfa_writer_versionstr(gw->version),
      minlen,
      inputfilename,
      has_containments ? '1' : '0',
      has_transitives ? '1' : '0');
  return 0;
}

static inline void gt_gfa_writer_show_gfa1_segment(GtFile *file,
    GtUword seqnum, const char *sequence)
{
  gt_file_xprintf(file, "S\t"GT_WU"\t%s\n", seqnum, sequence);
}

static inline void gt_gfa_writer_show_gfa2_segment(GtFile *file,
    GtUword seqnum, const char *sequence)
{
  gt_file_xprintf(file, "S\t"GT_WU"\t"GT_WU"\t%s\n", seqnum,
      (GtUword)strlen(sequence), sequence);
}

static inline void gt_gfa_writer_show_segment(GtGfaVersion version,
    GtFile *file, GtUword seqnum, const char *sequence)
{
    switch (version) {
      case GT_GFA_VERSION_1_0:
        gt_gfa_writer_show_gfa1_segment(file, seqnum, sequence);
        break;
      case GT_GFA_VERSION_2_0:
        gt_gfa_writer_show_gfa2_segment(file, seqnum, sequence);
        break;
      default:
        gt_assert(false);
    }
}

int gt_gfa_writer_show_segments(GtGfaWriter *gw, GT_UNUSED GtError *err)
{
  const GtTwobitencoding *nextencoded;
  GtTwobitencoding code = 0;
  GtUword seqnum = 0, nofseqs = gt_encseq_num_of_sequences(gw->encseq),
                pos = 0, next_stop, i,
                tlen = gt_encseq_total_length(gw->encseq), charsincode = 0;
  char *seqbuffer;
  bool last_seq = false;
  const char code2char[] = "ACGT";
  gt_assert(gw != NULL);
  seqbuffer = gt_malloc(sizeof (*seqbuffer) *
      (gt_encseq_max_seq_length(gw->encseq) + 1UL));
  nextencoded = gt_encseq_twobitencoding_export(gw->encseq);
  i = 0;
  while (!last_seq)
  {
    if (seqnum + 1UL == nofseqs)
    {
      last_seq = true;
      next_stop = tlen;
    }
    else
      next_stop = gt_encseq_seqstartpos(gw->encseq, seqnum + 1UL) - 1UL;
    for (/**/; i < next_stop; i++)
    {
      if (charsincode == 0)
      {
        code = *(nextencoded++);
        charsincode = (GtUword)GT_UNITSIN2BITENC;
      }
      seqbuffer[pos++] = code2char[code >> ((--charsincode) << 1) & 3];
    }
    seqbuffer[pos] = '\0';
    gt_gfa_writer_show_segment(gw->version, gw->file, seqnum, seqbuffer);
    pos = 0;
    /* consume separator */
    i++;
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (GtUword)GT_UNITSIN2BITENC;
    }
    charsincode--;
    seqnum++;
  }
  gt_free(seqbuffer);
  return 0;
}

static inline void gt_gfa_writer_show_edge_gfa1(GtFile *file,
    GtUword seqnum1, bool seq1_direct,
    GtUword seqnum2, bool seq2_direct,
    GtUword spmlen)
{
  gt_file_xprintf(file, "L\t"GT_WU"\t%c\t"GT_WU"\t%c\t"
                        GT_WU"M\n",
                        seqnum1, seq1_direct ? '+' : '-',
                        seqnum2, seq2_direct ? '+' : '-',
                        spmlen);
}

static inline void gt_gfa_writer_show_edge_gfa2(GtFile *file,
    GtUword seqnum1, bool seq1_direct,
    GtUword seqnum2, bool seq2_direct,
    GtUword start1, GtUword end1, bool end1_last,
    GtUword start2, GtUword end2, bool end2_last,
    GtUword spmlen)
{
  gt_file_xprintf(file, "E\t*\t"GT_WU"%c\t"GT_WU"%c\t"
                        GT_WU"\t"GT_WU"%s\t"GT_WU"\t"GT_WU"%s\t"
                        GT_WU"M\n",
                        seqnum1, seq1_direct ? '+' : '-',
                        seqnum2, seq2_direct ? '+' : '-',
                        start1, end1, end1_last ? "$" : "",
                        start2, end2, end2_last ? "$" : "",
                        spmlen);
}

static inline void gt_gfa_writer_show_edge(GtGfaVersion version, GtFile *file,
    GtUword seqnum1, bool seq1_direct, GtUword seqnum2, bool seq2_direct,
    GtUword start1, GtUword end1, bool end1_last,
    GtUword start2, GtUword end2, bool end2_last,
    GtUword spmlen)
{

  switch (version) {
    case GT_GFA_VERSION_1_0:
      gt_gfa_writer_show_edge_gfa1(file, seqnum1, seq1_direct, seqnum2,
          seq2_direct, spmlen);
      break;
    case GT_GFA_VERSION_2_0:
      gt_gfa_writer_show_edge_gfa2(file, seqnum1, seq1_direct, seqnum2,
          seq2_direct, start1, end1, end1_last, start2, end2, end2_last,
          spmlen);
      break;
    default:
      gt_assert(false);
  }
}

void gt_spmproc_show_gfa(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *gfa_writer)
{
  GtGfaWriter *gw = gfa_writer;
  const GtUword sl1 = gt_encseq_seqlength(gw->encseq, suffix_readnum),
                sl2 = gt_encseq_seqlength(gw->encseq, prefix_readnum);
  gt_gfa_writer_show_edge(gw->version, gw->file,
      suffix_readnum, suffixseq_direct,
      prefix_readnum, prefixseq_direct,
      suffixseq_direct ? sl1 - length : 0,
      suffixseq_direct ? sl1 : length,
      suffixseq_direct,
      prefixseq_direct ? 0 : sl2 - length,
      prefixseq_direct ? length : sl2,
      !prefixseq_direct,
      length);
}

void gt_gfa_writer_delete(GtGfaWriter *gw)
{
  gt_file_delete(gw->file);
  gt_encseq_delete(gw->encseq);
  gt_free(gw);
}
