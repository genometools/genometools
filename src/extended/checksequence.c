#ifndef CHECKSEQUENCE_H
#define CHECKSEQUENCE_H
#include "extended/checksequence.h"
#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

bool gap_symbol_in_sequence(const GtUchar *seq, GtUword len)
{
  const GtUchar *sptr;

  for (sptr = seq; sptr < seq + len; sptr++)
  {
    if (*sptr == LINEAR_EDIST_GAP)
    {
      return true;
    }
  }
  return false;
}
#endif
