#ifndef RECONSTRUCTALIGNMENT_H
#define RECONSTRUCTALIGNMENT_H
#include "core/types_api.h"
#include "extended/alignment.h"

void reconstructalignment(GtAlignment *align,
                          const GtUword *Ctab,
                          const GtUword vlen);

GtUword construct_trivial_alignment(GtAlignment *align, GtUword len,
                                    const GtWord gapcost,
                                    void (*indel)(GtAlignment*));
#endif
