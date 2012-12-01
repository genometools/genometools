#include "core/ma_api.h"
#include "extended/rbtree.h"
#include "extended/ranked_list.h"

struct GtRankedlist
{
  unsigned long currentsize,           /* current size of the ranked list */
                maxsize;               /* maximal size of the ranked list */
  GtRBTree  *root;                     /* root of tree */
  GtRBTreeCompareFunc comparefunction; /* comparefunction */
  void *worstelement,                  /* reference to worst key */
       *compareinfo;                   /* info needed by compare function */
};

GtRankedlist *gt_ranked_list_new(unsigned long maxsize,
                                 GtRBTreeCompareFunc comparefunction,
                                 void *compareinfo)
{
  GtRankedlist *ranked_list;

  ranked_list = gt_malloc(sizeof (*ranked_list));
  ranked_list->currentsize = 0;
  ranked_list->maxsize = maxsize;
  ranked_list->comparefunction = comparefunction;
  ranked_list->compareinfo = compareinfo;
  ranked_list->root = gt_rbtree_new(comparefunction, NULL, NULL);
  ranked_list->worstelement = NULL;
  return ranked_list;
}

void gt_ranked_list_insert(GtRankedlist *ranked_list,void *elemin)
{
  bool nodecreated = false;

  if (ranked_list->currentsize < ranked_list->maxsize)
  {
    if (ranked_list->currentsize == 0 ||
        ranked_list->comparefunction(elemin,ranked_list->worstelement,
                                     ranked_list->compareinfo) < 0)
    {
      ranked_list->worstelement = elemin;
    }
    (void) gt_rbtree_search(ranked_list->root, elemin, &nodecreated);
    if (nodecreated)
    {
      ranked_list->currentsize++;
    }
  } else
  {
/*
  new element is not as bad as worst element, so insert it and
  and delete the worst element
*/
    if (ranked_list->comparefunction(ranked_list->worstelement,elemin,
                                     ranked_list->compareinfo) < 0)
    {
      (void) gt_rbtree_search(ranked_list->root, elemin, &nodecreated);
      if (nodecreated)
      {
        if (gt_rbtree_erase(ranked_list->root, ranked_list->worstelement) != 0)
        {
          fprintf(stderr,"%s: deletion failed\n",__func__);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        ranked_list->worstelement = gt_rbtree_minimum_key(ranked_list->root);
      }
    }
  }
}

void *gt_ranked_list_minimum_key(const GtRankedlist *ranked_list)
{
  return gt_rbtree_minimum_key(ranked_list->root);
}

void gt_ranked_list_delete(GtRankedlist *ranked_list)
{
  gt_rbtree_delete(ranked_list->root);
  gt_free(ranked_list);
}
