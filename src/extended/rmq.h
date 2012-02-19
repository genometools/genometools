/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Simon J. Puglisi <simon.puglisi@rmit.edu.au>
*/

#ifndef RMQ_H
#define RMQ_H

#include "core/error_api.h"
#include "core/range_api.h"

typedef struct GtRMQ GtRMQ;

GtRMQ*        gt_rmq_new(long* arr, unsigned long size);
unsigned long gt_rmq_find_min_index(GtRMQ *rmq, unsigned long start,
                                    unsigned long end);
void          gt_rmq_delete(GtRMQ *rmq);

int           gt_rmq_unit_test(GtError *err);
#endif
