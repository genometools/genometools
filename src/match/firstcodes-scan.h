#ifndef FIRSTCODES_SCAN_H
#define FIRSTCODES_SCAN_H

#include "core/intbits.h"
#include "core/codetype.h"
#include "core/encseq_api.h"

void gt_firstcodes_kmerscan(const GtBitsequence *twobitencoding,
                            bool withcheck,
                            unsigned long equallength,
                            unsigned long totallength,
                            unsigned long maxunitindex,
                            unsigned int kmersize,
                            void (*processcode)(GtCodetype,GtCodetype,
                                                unsigned long,void *),
                            void *data);

void gt_firstcode_runkmerscan(const GtEncseq *encseq,
                              bool withcheck,unsigned int kmersize);

#endif
