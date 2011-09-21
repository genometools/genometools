#ifndef SEQNUMRELPOS_H
#define SEQNUMRELPOS_H

typedef struct GtSeqnumrelpostab GtSeqnumrelpostab;

GtSeqnumrelpostab *gt_seqnumrelpostab_new(unsigned long maxbucketsize,
                                          unsigned int bitsforrelpos,
                                          const GtEncseq *encseq);

void gt_seqnumrelpostab_delete(GtSeqnumrelpostab *snrp);

unsigned long gt_seqnumrelpostab_pos(const GtSeqnumrelpostab *snrp,
                                     unsigned long idx);

void gt_seqnumrelpostab_add(GtSeqnumrelpostab *snrp,
                            unsigned long idx,
                            unsigned long value);

unsigned long gt_seqnumrelpostab_encode(const GtSeqnumrelpostab *snrp,
                                        unsigned long seqnum,
                                        unsigned long relpos);

#endif
