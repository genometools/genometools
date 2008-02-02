#include <inttypes.h>
#include "libgtcore/strarray.h"
#include "libgtcore/seqiterator.h"
#include "spacedef.h"

typedef struct
{
  uint64_t unitnum;
  SeqIterator *seqit;
  bool newseq;
  char *desc;
} Substriter;

typedef struct
{
  const Uchar *queryptr;
  unsigned long remaining;
} Substring;

Substriter *substriter_new(const StrArray *queryfilenames,
                           const Uchar *symbolmap)
{
  Substriter *substriter;
  ALLOCASSIGNSPACE(substriter,NULL,Substriter,1);
  substriter->unitnum = 0;
  substriter->seqit 
    = seqiterator_new(queryfilenames,symbolmap,true);
  substriter->newseq = true;
  substriter->desc = NULL;
  return substriter;
}

int substriter_next(Substring *substring,Substriter *substriter,Error *err)
{
  while(true)
  {
    if (substriter->newseq)
    {
      int retval;
  
      retval = seqiterator_next(substriter->seqit,
                                &substring->queryptr,
                                &substring->remaining,
                                &substriter->desc,
                                err);
      if (retval <= 0)
      {
        FREESPACE(substriter->desc);
        return retval;
      }
      assert(substring->remaining > 0);
      substriter->newseq = false;
      break;
    }
    assert(substring->remaining > 0);
    substring->remaining--;
    substring->queryptr++;
    if (substring->remaining > 0)
    {
      break;
    }
    substriter->newseq = true;
    FREESPACE(substriter->desc);
  }
  return 1;
}
