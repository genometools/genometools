#ifndef LINEARSPACE_H
#define LINEARSPACE_H

#include "core/unused_api.h"
#include "core/error.h"
#include "extended/alignment.h"

void gt_computelinearspace2(bool showevalue,
                            const GtUchar *useq,
                             GtUword ulen,
                             const GtUchar *vseq,
                             GtUword vlen,
                             const GtWord matchcost,
                             const GtWord mismatchcost,
                             const GtWord gapcost,
                             FILE *fp);

GtUword gt_calc_linearalign_with_costs(const GtUchar *useq, GtUword ulen,
                            const GtUchar *vseq, GtUword vlen,
                            GtAlignment *align,
                            const GtWord matchcost,
                            const GtWord mismatchcost,
                            const GtWord gapcost);
/* extreme Codecuplizierung!!! entspricht linearedist mit varibalen Kosten*/
#endif
