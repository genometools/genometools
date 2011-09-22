#include "core/encseq_api.h"
#include "core/ma.h"
#include "seqnumrelpos.h"

struct GtSeqnumrelpos
{
  const GtEncseq *encseq;
  unsigned long relposmask;
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

unsigned long gt_seqnumrelpos_decode_pos(const GtSeqnumrelpos *snrp,
                                         unsigned long seqnumrelpos)
{
  unsigned long seqnum, relpos;

  seqnum = seqnumrelpos >> snrp->bitsforrelpos;
  relpos = seqnumrelpos & snrp->relposmask;
  return gt_encseq_seqstartpos(snrp->encseq,seqnum) + relpos;
}

unsigned long gt_seqnumrelpos_decode_seqnum(const GtSeqnumrelpos *snrp,
                                            unsigned long seqnumrelpos)
{
  return seqnumrelpos >> snrp->bitsforrelpos;
}

unsigned long gt_seqnumrelpos_decode_relpos(const GtSeqnumrelpos *snrp,
                                            unsigned long seqnumrelpos)
{
  return seqnumrelpos & snrp->relposmask;
}

unsigned long gt_seqnumrelpos_encode(const GtSeqnumrelpos *snrp,
                                     unsigned long seqnum,
                                        unsigned long relpos)
{
  gt_assert(relpos <= snrp->relposmask);
  return (seqnum << snrp->bitsforrelpos) | relpos;
}
