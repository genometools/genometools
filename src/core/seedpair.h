#ifndef SEEDPAIR_H
#define SEEDPAIR_H
#include <stdint.h>

typedef uint32_t GtSeedExtendPosition;
typedef uint32_t GtSeedExtendSeqnum;

typedef struct {
  GtSeedExtendSeqnum bseqnum; /*  2nd important sort criterion */
  GtSeedExtendSeqnum aseqnum; /* most important sort criterion */
  GtSeedExtendPosition bpos;
  GtSeedExtendPosition apos;  /*  3rd important sort criterion */
} GtSeedExtendSeedPair;

typedef struct {
  long double a;
  GtSeedExtendPosition bpos;
} GtQuasiSeedExtendSeedPair;

#endif
