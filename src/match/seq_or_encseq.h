#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"
#define GT_NO_QUERY(READMODE,SELFMATCH)\
        ((!GT_ISDIRREVERSE(READMODE) && (SELFMATCH)) ? true : false)

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
} GtSeqorEncseq;

#endif
