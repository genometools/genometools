/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/encseq_api.h"
#include "core/range_api.h"

/* These functions wrap <GtEncseq> accessors without using anything else
   than pointers and ints as parameters or return values. This is needed to
   work around the broken 64bit support in Ruby::DL.
   May be removed when we move on to something better. */

GtUchar gt_encseq_get_encoded_char_p(const GtEncseq *encseq,
                                     unsigned long *pos,
                                     GtReadmode readmode)
{
  return gt_encseq_get_encoded_char(encseq, *pos, readmode);
}

char gt_encseq_get_decoded_char_p(const GtEncseq *encseq,
                                  unsigned long *pos,
                                  GtReadmode readmode)
{
  return gt_encseq_get_decoded_char(encseq, *pos, readmode);
}

void gt_encseq_extract_encoded_p(const GtEncseq *encseq,
                                   GtUchar *buffer,
                                   GtRange *rng)
{
  gt_encseq_extract_encoded(encseq, buffer, rng->start, rng->end);
}

void gt_encseq_extract_decoded_p(const GtEncseq *encseq,
                                 char *buffer,
                                 GtRange *rng)
{
  gt_encseq_extract_decoded(encseq, buffer, rng->start, rng->end);
}

void gt_encseq_total_length_p(const GtEncseq *encseq,
                              unsigned long *totallength)
{
  *totallength = gt_encseq_total_length(encseq);
}

void gt_encseq_seqlength_p(const GtEncseq *encseq,
                           unsigned long *seqnumber,
                           unsigned long *seqlength)
{
  *seqlength = gt_encseq_seqlength(encseq, *seqnumber);
}

void gt_encseq_seqnum_p(const GtEncseq *encseq,
                        unsigned long *seqnumber,
                        unsigned long *pos)
{
  *seqnumber = gt_encseq_seqnum(encseq, *pos);
}

void gt_encseq_filenum_p(const GtEncseq *encseq,
                         unsigned long *filenumber,
                         unsigned long *pos)
{
  *filenumber = gt_encseq_filenum(encseq, *pos);
}

void gt_encseq_seqstartpos_p(const GtEncseq *encseq,
                             unsigned long *seqnumber,
                             unsigned long *startpos)
{
  *startpos = gt_encseq_seqstartpos(encseq, *seqnumber);
}

void gt_encseq_filestartpos_p(const GtEncseq *encseq,
                              unsigned long *filenumber,
                              unsigned long *startpos)
{
  *startpos = gt_encseq_seqstartpos(encseq, *filenumber);
}

void gt_encseq_num_of_files_p(const GtEncseq *encseq,
                              unsigned long *numoffiles)
{
  *numoffiles = gt_encseq_num_of_files(encseq);
}

void gt_encseq_num_of_sequences_p(const GtEncseq *encseq,
                                  unsigned long *numofseqs)
{
  *numofseqs = gt_encseq_num_of_sequences(encseq);
}

const char* gt_encseq_description_p(const GtEncseq *encseq,
                                    unsigned long *desclen,
                                    unsigned long *seqnum)
{
  return gt_encseq_description(encseq, desclen, *seqnum);
}

void gt_encseq_effective_filelength_p(const GtEncseq *encseq,
                                      uint64_t *result,
                                      unsigned long *filenum)
{
  *result = gt_encseq_effective_filelength(encseq, *filenum);
}

GtEncseqReader* gt_encseq_create_reader_with_readmode_p(const GtEncseq *encseq,
                                                        GtReadmode readmode,
                                                        unsigned long *pos)
{
  return gt_encseq_create_reader_with_readmode(encseq, readmode, *pos);
}
