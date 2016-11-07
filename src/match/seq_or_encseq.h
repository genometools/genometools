#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"

#define GT_QUERYSEQORENCSEQ_INIT_ENCSEQ(QUERYES,ENCSEQ,READMODE,SELFMATCH)\
        (QUERYES).encseq = ENCSEQ;\
        (QUERYES).seq = NULL;\
        (QUERYES).desc = NULL;\
        (QUERYES).no_query = (!GT_ISDIRREVERSE(READMODE) && (SELFMATCH)) \
                                ? true : false

#define GT_QUERYSEQORENCSEQ_INIT_SEQ(QUERYES,SEQ,SEQDESC,READMODE,SELFMATCH)\
        (QUERYES).encseq = NULL;\
        (QUERYES).seq = SEQ;\
        (QUERYES).desc = SEQDESC;\
        (QUERYES).no_query = (!GT_ISDIRREVERSE(READMODE) && (SELFMATCH)) \
                                ? true : false

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
  bool no_query;
} GtSeqorEncseq;

#endif
