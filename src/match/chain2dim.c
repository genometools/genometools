#include "chaindef.h"

/*
  The basic information required for each fragment is stored
  in a structure of the following type. The user has to specify
  those components which a tagged `user defined'. The chaining
  algorithms computes the remaining components score and previous
  in chain.
*/

typedef struct
{
  Seqpos startpos[2],  /* start of fragments in the 2 dimensions, userdef */
         endpos[2],    /* end of fragments in the 2 dimensions, userdef */
         firstinchain,   /* first element in chain, compute */
         previousinchain;  /* previous index in chain, compute */
  GtChainscoretype
         weight, /* weight of fragment, user defined */
         initialgap, /* gap to start of sequences, user defined */
         terminalgap, /* gap to last positions of fragment, user defined */
         score; /* score of highest scoreing chain ending here, compute */
} GtFragmentinfo;

struct GtFragmentinfotable
{
   GtFragmentinfo *fragments;
   unsigned long nextfree, allocated;
};

GtFragmentinfotable *fragmentinfotable_new(unsigned long numberoffragments)
{
  GtFragmentinfotable *fragmentinfotable
    = gt_malloc(sizeof(*fragmentinfotable));
  fragmentinfotable->fragments
    = gt_malloc(sizeof(*fragmentinfotable->fragments) * numberoffragments);
  fragmentinfotable->nextfree = 0;
  fragmentinfotable->allocated = numberoffragments;
  return fragmentinfotable;
}

void fragmentinfotable_delete(GtFragmentinfotable *fragmentinfotable)
{
  gt_free(fragmentinfotable->fragments);
  gt_free(fragmentinfotable);
}

void fragmentinfotable_add(GtFragmentinfotable *fragmentinfotable,
                           Seqpos start1,Seqpos end1,
                           Seqpos start2,Seqpos end2,
                           GtChainscoretype initialgap,
                           GtChainscoretype terminalgap,
                           GtChainscoretype weight)
{
  GtFragmentinfo *frag;

  gt_assert(fragmentinfotable->nextfree < fragmentinfotable->allocated);
  frag = fragmentinfotable->fragments + fragmentinfotable->nextfree;
  frag->startpos[0] = start1;
  frag->startpos[1] = start2;
  frag->endpos[0] = end1;
  frag->endpos[1] = end2;
  frag->initialgap = initialgap;
  frag->terminalgap = terminalgap;
  frag->weight = weight;
}
