#ifndef LINEARESPACE_ALIGN_H
#define LINEARSPACE_ALIGN_H

#include "core/unused_api.h"
#include "core/error.h"
#include "extended/alignment.h"

void gt_computelinearspace_with_output(const GtUchar *useq,
                                       GtUword ulen,
                                       const GtUchar *vseq,
                                       GtUword vlen,
                                       GtWord matchcost,
                                       GtWord mismatchcost,
                                       GtWord gapcost);
                           
GtUword gt_computelinearspace_with_costs(const GtUchar *u, GtUword ulen,
                                       const GtUchar *v, GtUword vlen,
                                       GtAlignment *align,
                                       const GtWord matchcost,
                                       const GtWord mismatchcost,
                                       const GtWord gapcost);
                                       
/* extreme Codecuplizierung!!! entspricht linearedist mit varibalen Kosten*/                             
#endif
