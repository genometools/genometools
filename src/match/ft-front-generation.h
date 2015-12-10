#ifndef FT_FRONT_GENERATION_H
#define FT_FRONT_GENERATION_H
#include <stdint.h>
#include "core/types_api.h"
#ifndef OUTSIDE_OF_GT
#include "core/arraydef.h"
#else
typedef struct
{
  GtUword nextfreeuint8_t,
          allocateduint8_t;
  uint8_t *spaceuint8_t;
} GtArrayuint8_t;

#endif

#define FT_EOP_REPLACEMENT 1
#define FT_EOP_INSERTION   (1 << 1)
#define FT_EOP_DELETION    (1 << 2)

#ifndef OUTSIDE_OF_GT
#define FT_EOPCODE_MAXREPLACEMENT 254
#define FT_EOPCODE_DELETION       254
#define FT_EOPCODE_INSERTION      255
#endif

typedef struct
{
  GtUword alignedlen, row, distance, trimleft;
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
                         GtArrayuint8_t *eoplist,
                         GtFronttrace *front_trace,
                         const Polished_point *pp,
                         GtUword pol_size,
                         GtWord match_score,
                         GtWord difference_score,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen);

void front_trace_multireplacement(GtArrayuint8_t *eoplist,GtUword repnum);

#endif
