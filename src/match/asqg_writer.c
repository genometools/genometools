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
#include "match/asqg_writer.h"

struct GtAsqgWriter
{
  GtFile *file;
  const GtEncseq *encseq;
};

GtAsqgWriter* gt_asqg_writer_new(GtFile *file, const GtEncseq *encseq)
{
  GtAsqgWriter *aw;
  aw = gt_malloc(sizeof (*aw));
  aw->file = file;
  aw->encseq = encseq;
  return aw;
}

int gt_asqg_writer_show_header(GtAsqgWriter *aw, float erate,
    unsigned long minlen, const char *inputfilename, bool has_containments,
    bool has_transitives, GT_UNUSED GtError *err)
{
  gt_assert(aw != NULL);
  gt_file_xprintf(aw->file,
      "HT\tVN:i:%lu\tER:f:%g\tOL:i:%lu\tIN:Z:%s\tCN:i:%c\tTE:i:%c\n",
      GT_ASQG_VERSION, erate, minlen, inputfilename,
      has_containments ? '1' : '0', has_transitives ? '1' : '0');
  return 0;
}

static inline void gt_asqg_writer_show_vertex_line(GtFile *file,
    unsigned long seqnum, const char *sequence, bool subsequence)
{
  gt_file_xprintf(file, "VT\t%lu\t%s\tSS:i:%c\n", seqnum, sequence,
      subsequence ? '1' : '0');
}

int gt_asqg_writer_show_vertices(GtAsqgWriter *aw, GT_UNUSED GtError *err)
{
  const GtTwobitencoding *nextencoded;
  GtTwobitencoding code = 0;
  unsigned long seqnum = 0, nofseqs = gt_encseq_num_of_sequences(aw->encseq),
                pos = 0, next_stop, i,
                tlen = gt_encseq_total_length(aw->encseq), charsincode = 0;
  char *seqbuffer;
  bool last_seq = false;
  const char code2char[] = "acgt";
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
      next_stop = tlen - 1UL;
    }
    else
      next_stop = gt_encseq_seqstartpos(aw->encseq, seqnum + 1UL) - 1UL;
    for (/**/; i < next_stop; i++)
    {
      if (charsincode == 0)
      {
        code = *(nextencoded++);
        charsincode = (unsigned long)GT_UNITSIN2BITENC;
      }
      seqbuffer[pos++] = code2char[code >> ((--charsincode) << 1) & 3];
    }
    seqbuffer[pos] = '\0';
    gt_asqg_writer_show_vertex_line(aw->file, seqnum, seqbuffer, false);
    pos = 0;
    /* consume separator */
    i++;
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (unsigned long)GT_UNITSIN2BITENC;
    }
    charsincode--;
    seqnum++;
  }
  gt_free(seqbuffer);
  return 0;
}

static inline void gt_asqg_writer_show_edge_line(GtFile *file,
    unsigned long seqnum1, unsigned long seqnum2, unsigned long start1,
    unsigned long end1, unsigned long seqlen1, unsigned long start2,
    unsigned long end2, unsigned long seqlen2, bool revcompl,
    unsigned long edist)
{
  gt_file_xprintf(file, "ED\t%lu %lu %lu %lu %lu %lu %lu %lu %c %lu\n",
      seqnum1, seqnum2, start1, end1, seqlen1, start2, end2, seqlen2,
      revcompl ? '1' : '0', edist);
}

void gt_spmproc_show_asgq(unsigned long suffix_readnum,
    unsigned long prefix_readnum, unsigned long length,
    bool suffixseq_direct, bool prefixseq_direct, void *asqg_writer)
{
  GtAsqgWriter *aw = asqg_writer;
  const unsigned long sl1 = gt_encseq_seqlength(aw->encseq, suffix_readnum),
                      sl2 = gt_encseq_seqlength(aw->encseq, prefix_readnum);
  gt_asqg_writer_show_edge_line(aw->file, suffix_readnum, prefix_readnum,
      suffixseq_direct ? sl1 - length : 0,
      suffixseq_direct ? sl1 - 1UL : length - 1UL, sl1,
      prefixseq_direct ? 0 : sl2 - length,
      prefixseq_direct ? length - 1UL : sl2 - 1UL, sl2,
      !suffixseq_direct || !prefixseq_direct, 0);
}

void gt_asqg_writer_delete(GtAsqgWriter *aw)
{
  gt_free(aw);
}
