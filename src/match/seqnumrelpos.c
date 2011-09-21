#include "core/encseq_api.h"
#include "core/ma.h"
#include "seqnumrelpos.h"

struct GtSeqnumrelpostab
{
  const GtEncseq *encseq;
  unsigned long relposmask, *tab;
  unsigned int bitsforrelpos;
};

GtSeqnumrelpostab *gt_seqnumrelpostab_new(unsigned long maxbucketsize,
                                          unsigned int bitsforrelpos,
                                          const GtEncseq *encseq)
{
  GtSeqnumrelpostab *snrp;

  snrp = gt_malloc(sizeof (*snrp));
  snrp->tab = gt_malloc(sizeof (*snrp->tab) * maxbucketsize);
  snrp->bitsforrelpos = bitsforrelpos;
  snrp->relposmask = (1UL << snrp->bitsforrelpos) - 1;
  snrp->encseq = encseq;
  return snrp;
}

void gt_seqnumrelpostab_delete(GtSeqnumrelpostab *snrp)
{
  if (snrp != NULL)
  {
    gt_free(snrp->tab);
    gt_free(snrp);
  }
}

unsigned long gt_seqnumrelpostab_pos(const GtSeqnumrelpostab *snrp,
                                     unsigned long idx)
{
  unsigned long seqnum, relpos;

  seqnum = snrp->tab[idx] >> snrp->bitsforrelpos;
  relpos = snrp->tab[idx] & snrp->relposmask;
  return gt_encseq_seqstartpos(snrp->encseq,seqnum) + relpos;
}

void gt_seqnumrelpostab_add(GtSeqnumrelpostab *snrp,
                            unsigned long idx,
                            unsigned long value)
{
  snrp->tab[idx] = value;
}
