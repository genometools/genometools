#ifndef SEQ_OR_ENCSEQ_H
#define SEQ_OR_ENCSEQ_H
#include "core/assert_api.h"
#include "core/encseq_api.h"

#define GT_SEQORENCSEQ_INIT_ENCSEQ(SORE,ENCSEQ)\
        (SORE)->encseq = ENCSEQ;\
        (SORE)->seq = NULL;\
        (SORE)->desc = NULL;\
        (SORE)->seqlength = GT_UWORD_MAX;\
        (SORE)->seqstartpos = GT_UWORD_MAX;\
        (SORE)->characters = NULL;\
        (SORE)->wildcardshow = 0;\
        (SORE)->haswildcards = true

#define GT_SEQORENCSEQ_ADD_SEQ_COORDS(SORE,SEQSTARTPOS,SEQLENGTH)\
        gt_assert((SORE)->encseq != NULL);\
        (SORE)->seqstartpos = SEQSTARTPOS;\
        (SORE)->seqlength = SEQLENGTH

#define GT_SEQORENCSEQ_INIT_SEQ(SORE,SEQ,SEQDESC,SEQLENGTH,CHARACTERS,\
                                WILDCARDSHOW,HASWILDCARDS)\
        (SORE)->encseq = NULL;\
        (SORE)->seq = SEQ;\
        (SORE)->desc = SEQDESC;\
        (SORE)->seqlength = SEQLENGTH;\
        (SORE)->seqstartpos = 0;\
        (SORE)->characters = CHARACTERS;\
        (SORE)->wildcardshow = WILDCARDSHOW;\
        (SORE)->haswildcards = HASWILDCARDS

typedef struct
{
  const GtUchar *seq;
  const GtEncseq *encseq;
  const char *desc; /* only used if seq != NULL and display_seq_desc */
  GtWord seqstartpos, seqlength;
  const GtUchar *characters;
  GtUchar wildcardshow;
  bool haswildcards;
} GtSeqorEncseq;

#endif
