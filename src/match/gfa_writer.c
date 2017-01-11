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
  const GtEncseq *encseq;
};

GtGfaWriter* gt_gfa_writer_new(GtFile *file, const GtEncseq *encseq)
{
  GtGfaWriter *aw;
  aw = gt_malloc(sizeof (*aw));
  aw->file = file;
  aw->encseq = encseq;
  return aw;
}

int gt_gfa_writer_show_header(GtGfaWriter *aw,
    GtUword minlen, const char *inputfilename, bool has_containments,
    bool has_transitives, GT_UNUSED GtError *err)
{
  gt_assert(aw != NULL);
  gt_file_xprintf(aw->file,
      "H\tVN:Z:"GT_GFA_VERSION"\n"
      "H\tpn:Z:readjoiner\n"
      "H\tol:i:"GT_WU"\n"
      "H\tin:Z:%s\n"
      "H\tcn:i:%c\n"
      "H\tte:i:%c\n",
      minlen,
      inputfilename,
      has_containments ? '1' : '0',
      has_transitives ? '1' : '0');
  return 0;
}

static inline void gt_gfa_writer_show_vertex_line(GtFile *file,
    GtUword seqnum, const char *sequence, GT_UNUSED bool subsequence)
{
  gt_file_xprintf(file, "S\t"GT_WU"\t"GT_WU"\t%s\n", seqnum,
      strlen(sequence), sequence);
}

int gt_gfa_writer_show_vertices(GtGfaWriter *aw, GT_UNUSED GtError *err)
{
  const GtTwobitencoding *nextencoded;
  GtTwobitencoding code = 0;
  GtUword seqnum = 0, nofseqs = gt_encseq_num_of_sequences(aw->encseq),
                pos = 0, next_stop, i,
                tlen = gt_encseq_total_length(aw->encseq), charsincode = 0;
  char *seqbuffer;
  bool last_seq = false;
  const char code2char[] = "ACGT";
  gt_assert(aw != NULL);
  seqbuffer = gt_malloc(sizeof (*seqbuffer) *
      (gt_encseq_max_seq_length(aw->encseq) + 1UL));
  nextencoded = gt_encseq_twobitencoding_export(aw->encseq);
  i = 0;
  while (!last_seq)
  {
    if (seqnum + 1UL == nofseqs)
    {
      last_seq = true;
      next_stop = tlen;
    }
    else
      next_stop = gt_encseq_seqstartpos(aw->encseq, seqnum + 1UL) - 1UL;
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
    gt_gfa_writer_show_vertex_line(aw->file, seqnum, seqbuffer, false);
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

static inline void gt_gfa_writer_show_edge_line(GtFile *file,
    GtUword seqnum1, bool seq1_direct,
    GtUword seqnum2, bool seq2_direct,
    GtUword start1, GtUword end1,
    GtUword start2, GtUword end2,
    GtUword spmlen)
{
  gt_file_xprintf(file, "E\t*\t"GT_WU"%c\t"GT_WU"%c\t"
                        GT_WU"\t"GT_WU"\t"GT_WU"\t"GT_WU"\t"
                        GT_WU"M\n",
                        seqnum1, seq1_direct ? '+' : '-',
                        seqnum2, seq2_direct ? '+' : '-',
                        start1, end1, start2, end2,
                        spmlen);
}

void gt_spmproc_show_gfa(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *gfa_writer)
{
  GtGfaWriter *aw = gfa_writer;
  const GtUword sl1 = gt_encseq_seqlength(aw->encseq, suffix_readnum),
                sl2 = gt_encseq_seqlength(aw->encseq, prefix_readnum);
  gt_gfa_writer_show_edge_line(aw->file,
      suffix_readnum, suffixseq_direct,
      prefix_readnum, prefixseq_direct,
      suffixseq_direct ? sl1 - length : 0,
      suffixseq_direct ? sl1 : length,
      prefixseq_direct ? 0 : sl2 - length,
      prefixseq_direct ? length : sl2,
      length);
}

void gt_gfa_writer_delete(GtGfaWriter *aw)
{
  gt_free(aw);
}
