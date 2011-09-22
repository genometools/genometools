#ifndef SEQNUMRELPOS_H
#define SEQNUMRELPOS_H

typedef struct GtSeqnumrelpos GtSeqnumrelpos;

GtSeqnumrelpos *gt_seqnumrelpos_new(unsigned int bitsforrelpos,
                                    const GtEncseq *encseq);

void gt_seqnumrelpos_delete(GtSeqnumrelpos *snrp);

unsigned long gt_seqnumrelpos_decode_pos(const GtSeqnumrelpos *snrp,
                                         unsigned long seqnumrelpos);

unsigned long gt_seqnumrelpos_decode_seqnum(const GtSeqnumrelpos *snrp,
                                            unsigned long seqnumrelpos);

unsigned long gt_seqnumrelpos_decode_relpos(const GtSeqnumrelpos *snrp,
                                            unsigned long seqnumrelpos);

unsigned long gt_seqnumrelpos_encode(const GtSeqnumrelpos *snrp,
                                     unsigned long seqnum,
                                     unsigned long relpos);

#endif
