#ifndef FT_FRONT_GENERATION_H
#define FT_FRONT_GENERATION_H
#include <stdint.h>
#include <stdbool.h>
#include "core/types_api.h"
#include "ft-eoplist.h"

#define FT_EOP_MISMATCH  1
#define FT_EOP_INSERTION (1 << 1)
#define FT_EOP_DELETION  (1 << 2)

typedef struct
{
  GtUword alignedlen, row, distance, trimleft, max_mismatches;
} Polished_point;

typedef struct GtFronttrace GtFronttrace;

GtFronttrace *front_trace_new(void);

void front_trace_delete(GtFronttrace *front_trace);

void front_trace_reset(GtFronttrace *front_trace,GtUword sumseqlen);

void front_trace_add_gen(GtFronttrace *front_trace,GtUword trimleft,
                         GtUword valid);

void front_trace_add_trace(GtFronttrace *front_trace,uint8_t backreference,
                           unsigned int lcs);

void front_trace2eoplist(bool polished,
                         GtEoplist *eoplist,
                         GtFronttrace *front_trace,
                         const Polished_point *pp,
                         GtUword pol_size,
                         GtWord match_score,
                         GtWord difference_score,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen);

#endif
