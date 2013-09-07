/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/encseq_api.h"
#include "core/ma.h"
#include "seqnumrelpos.h"

struct GtSeqnumrelpos
{
  const GtEncseq *encseq;
  GtUword relposmask;
  unsigned int bitsforrelpos;
};

GtSeqnumrelpos *gt_seqnumrelpos_new(unsigned int bitsforrelpos,
                                    const GtEncseq *encseq)
{
  GtSeqnumrelpos *snrp;

  snrp = gt_malloc(sizeof (*snrp));
  snrp->bitsforrelpos = bitsforrelpos;
  snrp->relposmask = (1UL << snrp->bitsforrelpos) - 1;
  snrp->encseq = encseq;
  return snrp;
}

void gt_seqnumrelpos_delete(GtSeqnumrelpos *snrp)
{
  gt_free(snrp);
}

GtUword gt_seqnumrelpos_decode_pos(const GtSeqnumrelpos *snrp,
                                         GtUword seqnumrelpos)
{
  GtUword seqnum, relpos;

  seqnum = (GtUword) seqnumrelpos >> (GtUword) snrp->bitsforrelpos;
  relpos = seqnumrelpos & snrp->relposmask;
  return gt_encseq_seqstartpos(snrp->encseq,seqnum) + relpos;
}

GtUword gt_seqnumrelpos_decode_seqnum(const GtSeqnumrelpos *snrp,
                                            GtUword seqnumrelpos)
{
  return (GtUword) seqnumrelpos >> (GtUword) snrp->bitsforrelpos;
}

GtUword gt_seqnumrelpos_decode_relpos(const GtSeqnumrelpos *snrp,
                                            GtUword seqnumrelpos)
{
  return seqnumrelpos & snrp->relposmask;
}

GtUword gt_seqnumrelpos_encode(const GtSeqnumrelpos *snrp,
                                     GtUword seqnum,
                                        GtUword relpos)
{
  gt_assert(relpos <= snrp->relposmask);
  return ((GtUword) seqnum << (GtUword) snrp->bitsforrelpos) | relpos;
}
