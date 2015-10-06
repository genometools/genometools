#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/encseq_api.h"
typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
} GtSeqorEncseq;
#endif
