#ifndef FT_FRONT_GENERATION_H
#define FT_FRONT_GENERATION_H
#include <stdint.h>
#include <stdbool.h>
#include "match/ft-eoplist.h"

#define FT_EOP_MISMATCH  1
#define FT_EOP_INSERTION (1 << 1)
#define FT_EOP_DELETION  (1 << 2)

typedef struct
{
  GtUword alignedlen, row, distance, trimleft, max_mismatches;
} GtFtPolished_point;

typedef struct GtFrontTrace GtFrontTrace;

GtFrontTrace *front_trace_new(void);

void front_trace_delete(GtFrontTrace *front_trace);

void front_trace_reset(GtFrontTrace *front_trace,GtUword sumseqlen);

void front_trace_add_gen(GtFrontTrace *front_trace,GtUword trimleft,
                         GtUword valid);

void front_trace_add_trace(GtFrontTrace *front_trace,uint8_t backreference,
                           uint32_t localmatch_count);

void front_trace2eoplist(bool polished,
                         GtEoplist *eoplist,
                         GtFrontTrace *front_trace,
                         const GtFtPolished_point *pp,
                         GtUword pol_size,
                         GtWord match_score,
                         GtWord difference_score,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen);

void gt_front_trace2eoplist_full_front_directed(GtEoplist *eoplist,
                                                const GtFrontTrace *front_trace,
                                                GtUword distance,
                                                const GtUchar *useq,
                                                GtUword ulen,
                                                const GtUchar *vseq,
                                                GtUword vlen);

#endif
