#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"

#define GT_SEQORENCSEQ_INIT_ENCSEQ(SORE,ENCSEQ)\
        (SORE)->encseq = ENCSEQ;\
        (SORE)->seq = NULL;\
        (SORE)->desc = NULL;\
        (SORE)->seqlength = 0

#define GT_SEQORENCSEQ_INIT_SEQ(SORE,SEQ,SEQDESC,SEQLENGTH)\
        (SORE)->encseq = NULL;\
        (SORE)->seq = SEQ;\
        (SORE)->desc = SEQDESC;\
        (SORE)->seqlength = SEQLENGTH

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
  GtWord seqlength;
} GtSeqorEncseq;

#endif
