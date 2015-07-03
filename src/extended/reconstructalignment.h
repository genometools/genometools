#ifndef RECONSTRUCT_H
#define RECONSTRUCT_H
#include "core/types_api.h"
#include "extended/alignment.h"

void reconstructalignment(GtAlignment *align,
                          const GtUword *Ctab,
                          const GtUword vlen);

GtUword determineCtab0(GtUword *Ctab,
                            const GtUchar vseq0,
                            const GtUchar *useq);
#endif
