#include "core/encseq_api.h"
#include "core/ma.h"
#include "seqnumrelpos.h"

struct GtSeqnumrelpostab
{
  const GtEncseq *encseq;
  unsigned long numofentries, relposmask, *tab;
  unsigned int bitsforrelpos;
};

GtSeqnumrelpostab *gt_seqnumrelpostab_new(unsigned long numofentries,
                                          unsigned int bitsforrelpos,
                                          const GtEncseq *encseq)
{
  GtSeqnumrelpostab *snrp;

  snrp = gt_malloc(sizeof (*snrp));
  snrp->tab = gt_malloc(sizeof (*snrp->tab) * numofentries);
  snrp->numofentries = numofentries;
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

  gt_assert(idx < snrp->numofentries);
  seqnum = snrp->tab[idx] >> snrp->bitsforrelpos;
  relpos = snrp->tab[idx] & snrp->relposmask;
  return gt_encseq_seqstartpos(snrp->encseq,seqnum) + relpos;
}

unsigned long gt_seqnumrelpostab_seqnum(const GtSeqnumrelpostab *snrp,
                                        unsigned long idx)
{
  gt_assert(idx < snrp->numofentries);
  return snrp->tab[idx] >> snrp->bitsforrelpos;
}

unsigned long gt_seqnumrelpostab_relpos(const GtSeqnumrelpostab *snrp,
                                        unsigned long idx)
{
  gt_assert(idx < snrp->numofentries);
  return snrp->tab[idx] & snrp->relposmask;
}

void gt_seqnumrelpostab_add(GtSeqnumrelpostab *snrp,
                            unsigned long idx,
                            unsigned long value)
{
  gt_assert(idx < snrp->numofentries);
  snrp->tab[idx] = value;
}

unsigned long gt_seqnumrelpostab_encode(const GtSeqnumrelpostab *snrp,
                                        unsigned long seqnum,
                                        unsigned long relpos)
{
  gt_assert(relpos <= snrp->relposmask);
  return (seqnum << snrp->bitsforrelpos) | relpos;
}
