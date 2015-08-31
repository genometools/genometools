#ifndef FT_FRONT_GENERATION_H
#define FT_FRONT_GENERATION_H
#include <stdint.h>
#include "core/types_api.h"
#ifndef OUTSIDE_OF_GT
#include "core/arraydef.h"
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

typedef struct Fronttrace Fronttrace;

Fronttrace *front_trace_new(void);

void front_trace_delete(Fronttrace *front_trace);

void front_trace_reset(Fronttrace *front_trace,GtUword sumseqlen);

void front_trace_add_gen(Fronttrace *front_trace,GtUword trimleft,
                         GtUword valid);

void front_trace_add_trace(Fronttrace *front_trace,uint8_t backreference,
                           unsigned int lcs);

#ifdef OUTSIDE_OF_GT
void front_trace_verify(const Fronttrace *front_trace,
                        const Polished_point *pp,
                        const GtUchar *useq,
                        GtUword ulen,
                        const GtUchar *vseq,
                        GtUword vlen);

void front_trace_verify_all(const Fronttrace *front_trace,
                            const Polished_point *pp,
                            const GtUchar *useq,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vlen);
#else
void front_trace_multireplacement(GtArrayuint8_t *eoplist,GtUword repnum);

void front_trace2eoplist(GtArrayuint8_t *eoplist,
                         const Fronttrace *front_trace,
                         const Polished_point *pp,
                         GT_UNUSED GtUword ulen,
                         GT_UNUSED GtUword vlen);

#endif
#endif
