#include <stdint.h>
#include "core/arraydef_api.h"
#include "match/rectangle-store.h"

struct GtArrayGtDiagbandseedRectangle
{
  GtDiagbandseedRectangle *spaceGtDiagbandseedRectangle;
  GtUword allocatedGtDiagbandseedRectangle,
          nextfreeGtDiagbandseedRectangle;
};

static int gt_rectangle_cmp(const GtDiagbandseedRectangle *l,
                            const GtDiagbandseedRectangle *r)
{
  if (l->a_end < r->a_end)
  {
    return -1;
  }
  if (l->a_end > r->a_end)
  {
    return 1;
  }
  if (l->b_end < r->b_end)
  {
    return -1;
  }
  if (l->b_end > r->b_end)
  {
    return 1;
  }
  /*
  fprintf(stderr,"%s: cmp(%u,%u)(%u,%u) vs (%u,%u)(%u,%u)\n",__func__,
                  l->a_start,l->a_end,l->b_start,l->b_end,
                  r->a_start,r->a_end,r->b_start,r->b_end);
  */
  return 0;
}

GtArrayGtDiagbandseedRectangle *gt_rectangle_store_new(void)
{
  GtArrayGtDiagbandseedRectangle *rectangle_store
    = gt_malloc(sizeof *rectangle_store);
  GT_INITARRAY(rectangle_store,GtDiagbandseedRectangle);
  return rectangle_store;
}

void gt_rectangle_store_delete(GtArrayGtDiagbandseedRectangle *rectangle_store)
{
  if (rectangle_store != NULL)
  {
    GT_FREEARRAY(rectangle_store,GtDiagbandseedRectangle);
    gt_free(rectangle_store);
  }
}

void gt_rectangle_store_add(GtArrayGtDiagbandseedRectangle *rectangle_store,
                            const GtDiagbandseedRectangle *key)
{
  GtDiagbandseedRectangle *ptr;

  GT_CHECKARRAYSPACE(rectangle_store,GtDiagbandseedRectangle,
                     256UL +
                     rectangle_store->allocatedGtDiagbandseedRectangle * 0.2);
  for (ptr = rectangle_store->spaceGtDiagbandseedRectangle +
             rectangle_store->nextfreeGtDiagbandseedRectangle;
       ptr > rectangle_store->spaceGtDiagbandseedRectangle &&
       gt_rectangle_cmp(ptr-1,key) > 0; ptr--)
  {
    *ptr = *(ptr-1);
  }
  *ptr = *key;
  rectangle_store->nextfreeGtDiagbandseedRectangle++;
}

const GtDiagbandseedRectangle *gt_rectangle_previous_equal(
                   const GtArrayGtDiagbandseedRectangle *rectangle_store,
                   const GtDiagbandseedRectangle *key)
{
  const GtDiagbandseedRectangle
    *midptr,
    *found = NULL,
    *leftptr = rectangle_store->spaceGtDiagbandseedRectangle,
    *rightptr = leftptr + rectangle_store->nextfreeGtDiagbandseedRectangle - 1;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + (GtUword) (rightptr-leftptr)/2;
    if (gt_rectangle_cmp(key,midptr) < 0)
    {
      rightptr = midptr - 1;
    } else
    {
      /* Now key *midptr <= key. Hence *midptr satisfies the
         requirement, but it may not be the largest element <= key. */
      found = midptr;
      /* if *midptr < key, then there may be a larger element and
         therefore we search for it in the range
         midptr+1 .. rightptr. if *midptr == key, then we have found
         the largest element and we can stop the search. */
      if (gt_rectangle_cmp(key,midptr) > 0)
      {
        leftptr = midptr + 1;
      } else
      {
        break;
      }
    }
  }
  /* if found == NULL, then we have never seen an element <= key and
     so we return NULL. Otherwise found has been set in the else case
     of the loop */
  return found;
}

bool gt_rectangle_overlap(
                   const GtArrayGtDiagbandseedRectangle *rectangle_store,
                   const GtDiagbandseedRectangle *key)
{
  const GtDiagbandseedRectangle *ptr;
/*
    = gt_rectangle_previous_equal(rectangle_store,key)
*/
   ;

/*
  if (ptr == NULL)
  {
    ptr = rectangle_store->spaceGtDiagbandseedRectangle;
  }
*/
  ptr = rectangle_store->spaceGtDiagbandseedRectangle;
  while (ptr < rectangle_store->spaceGtDiagbandseedRectangle +
               rectangle_store->nextfreeGtDiagbandseedRectangle)
  {
    if (ptr->a_start <= key->a_end && ptr->a_end >= key->a_start &&
        ptr->b_start <= key->b_end && ptr->b_end >= key->b_start)
    {
      return true;
    }
    ptr++;
  }
  return false;
}

void gt_rectangle_store_show(const GtArrayGtDiagbandseedRectangle
                                *rectangle_store)
{
  const GtDiagbandseedRectangle *key;

  printf("# %s\n",__func__);
  for (key = rectangle_store->spaceGtDiagbandseedRectangle;
       key < rectangle_store->spaceGtDiagbandseedRectangle +
             rectangle_store->nextfreeGtDiagbandseedRectangle; key++)
  {
     printf("# %u %u %u %u\n",key->a_start,key->a_end,
                              key->b_start,key->b_end);
  }
}
