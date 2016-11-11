#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"

#define GT_QUERYSEQORENCSEQ_INIT_ENCSEQ(QUERYES,ENCSEQ)\
        (QUERYES).encseq = ENCSEQ;\
        (QUERYES).seq = NULL;\
        (QUERYES).desc = NULL;\
        (QUERYES).seqlength = 0

#define GT_QUERYSEQORENCSEQ_INIT_SEQ(QUERYES,SEQ,SEQDESC,SEQLENGTH)\
        (QUERYES).encseq = NULL;\
        (QUERYES).seq = SEQ;\
        (QUERYES).desc = SEQDESC;\
        (QUERYES).seqlength = SEQLENGTH

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
  GtWord seqlength;
} GtSeqorEncseq;

#endif
