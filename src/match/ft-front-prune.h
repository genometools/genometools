#ifndef FT_FRONT_PRUNE_H
#define FT_FRONT_PRUNE_H
#include "core/types_api.h"
#include "ft-polish.h"
#include "ft-front-generation.h"
#include "core/encseq_api.h"

typedef struct
{
  void *space;
  GtUword offset, allocated;
} GtAllocatedMemory;

#ifndef OUTSIDE_OF_GT
typedef enum
{
  GT_EXTEND_CHAR_ACCESS_ENCSEQ,
  GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER,
  GT_EXTEND_CHAR_ACCESS_ANY
} GtExtendCharAccess;

typedef struct
{
  const GtEncseq *encseq;
  GtAllocatedMemory *sequence_cache;
  GtEncseqReader *encseq_r;
  GtExtendCharAccess extend_char_access;
  GtUword totallength;
} FTsequenceResources;
#endif

GtUword front_prune_edist_inplace(
#ifndef OUTSIDE_OF_GT
                       bool forward,
                       GtAllocatedMemory *frontspace_reservoir,
#endif
                       Polished_point *best_polished_point,
                       Fronttrace *fronttrace,
                       const Polishing_info *pol_info,
                       GtUword history,
                       GtUword minmatchnum,
                       GtUword maxalignedlendifference,
                       FTsequenceResources *ufsr,
                       GtUword ustart,
                       GtUword uulen,
                       FTsequenceResources *vfsr,
                       GtUword vstart,
                       GtUword vlen);

#endif
