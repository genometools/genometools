#ifndef FIRSTCODES_SCAN_H
#define FIRSTCODES_SCAN_H

#include "core/intbits.h"
#include "core/codetype.h"

void gt_firstcodes_kmerscan(const GtBitsequence *twobitencoding,
                            unsigned long equallength,
                            unsigned long totallength,
                            unsigned long maxunitindex,
                            unsigned int kmersize,
                            void (*processcode)(GtCodetype,GtCodetype,
                                                unsigned long,void *),
                            void *data);

#endif
