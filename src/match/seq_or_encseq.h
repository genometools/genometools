#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"

#define GT_QUERYSEQORENCSEQ_INIT_ENCSEQ(QUERYES,ENCSEQ,SELFMATCH)\
        (QUERYES).encseq = ENCSEQ;\
        (QUERYES).seq = NULL;\
        (QUERYES).desc = NULL;\
        (QUERYES).selfmatch = SELFMATCH

#define GT_QUERYSEQORENCSEQ_INIT_SEQ(QUERYES,SEQ,SEQDESC,SELFMATCH)\
        (QUERYES).encseq = NULL;\
        (QUERYES).seq = SEQ;\
        (QUERYES).desc = SEQDESC;\
        (QUERYES).SELFMATCH = SELFMATCH

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
  bool selfmatch;
} GtSeqorEncseq;

#endif
