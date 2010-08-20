/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

/*
 This module is derived from the file tsearch.c implementing
 redblack trees. Some bugs and inconsistencies
 have been fixed by Stefan Kurtz, <kurtz@zbh.uni-hamburg.de> and
 some functions have been added.

 Copyright (C) 1995, 1996, 1997, 2000 Free Software Foundation, Inc.
 This file is part of the GNU C Library.
 Contributed by Bernd Schmidt <crux@Pool.Informatik.RWTH-Aachen.DE>, 1997.

 The GNU C Library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 The GNU C Library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with the GNU C Library; if not, write to the Free
 Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 02111-1307 USA.

 Tree search for red/black trees. The algorithm for adding nodes is taken
 from one of the many "Algorithms" books by Robert Sedgewick, although the
 implementation differs. The algorithm for deleting nodes can probably be
 found in a book named "Introduction to Algorithms" by
 Cormen/Leiserson/Rivest.  At least that's the book that my professor took
 most algorithms from during the "Data Structures" course...

 Red/black trees are binary trees in which the edges are colored either
 red through or black. They have the following properties:
 1. The number of black edges on every path from the root to
 a leaf is constant.
 2. No two red edges are adjacent. Therefore there is an upper bound
 on the length of every path, it's O(log n) where n is the number
 of nodes in the tree. No path can be longer than 1+2*P where P
 is the length of the shortest path in the tree. Useful for the
 implementation:
 3. If one of the children of a node is NULL, then the other one
 is red (if it exists).

 In the implementation, not the edges are colored, but the nodes.
 The color interpreted as the color of the edge leading to this node.
 The color is meaningless for the root node, but we color the root
 node black for convenience. All added nodes are red initially.

 Adding to a red/black tree is rather easy. The right place is searched
 with a usual binary tree search. Additionally, whenever a node N is
 reached that has two red successors, the successors are colored black
 and the node itself colored red. This moves red edges up the tree
 where they pose less of a problem once we get to really insert the
 new node.  Changing N's color to red may violate rule 2, however, so
 rotations may become necessary to restore the invariants.
 Adding a new red leaf may violate the same rule, so
 afterwards an additional check is run and the tree possibly rotated.

 Deleting is hairy. There are mainly two nodes involved: the node to be
 deleted * (n1), and another node that is to be unchained from the tree
 (n2). If n1 has a successor (the node with a smallest key that is larger
 than n1), then the
 successor becomes n2 and its contents are copied into n1, otherwise n1
 becomes n2. Unchaining a node may violate rule 1: if n2 is black, one
 subtree is missing one black edge afterwards.  The algorithm must try to
 move this error upwards towards the root, so that the subtree that
 does not have enough black edges becomes the whole tree.
 Once that happens, the error has disappeared. It may not be necessary
 to go all the way up, since it is possible that rotations and
 recoloring can fix the error before that.

 Although the deletion algorithm must walk upwards through the tree,
 we do not store parent pointers in the nodes. Instead, delete allocates
 a small array of parent pointers and fills it while GT_RBT_DESCENDING the tree.
 Since we know that the length of a path is O(log n), where n is the
 number of nodes, this is likely to use less memory.

 Tree rotations look like this:
      A                C
     / \              / \
    B   C            A   G
   / \ / \  -->     / \
   D E F G         B   F
                  / \
                 D   E

 In this case, A has been rotated left.
 This preserves the ordering of the binary tree.
*/

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/mathsupport.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/stack-inlined.h"
#include "extended/redblack.h"

#define RBT_CHECK_RETURN_CODE\
        if (retcode < 0 || retcode == 1)\
        {\
          return retcode;\
        }

struct GtRBTnode
{
  bool red;
  GtKeytype key;
  struct GtRBTnode *left, *right;
};

#ifdef GT_RBT_DEBUG

/*
 Routines to check tree invariants.
 */

#define CHECK_TREE(NODE) check_tree(NODE)

static void check_tree_recurse(const GtRBTnode *p,
                               unsigned long d_sofar,
                               unsigned long d_total)
{
  if (p == NULL)
  {
    gt_assert(d_sofar == d_total);
    return;
  }
  check_tree_recurse (p->left,
                      d_sofar + (unsigned long) (p->left
                                               && !p->left->red),
                      d_total);
  check_tree_recurse (p->right,
                      d_sofar + (unsigned long) (p->right
                                                 && !p->right->
                                                 red), d_total);
  if (p->left != NULL)
  {
    gt_assert(!(p->left->red && p->red));
  }
  if (p->right != NULL)
  {
    gt_assert(!(p->right->red && p->red));
  }
}

static void check_tree(const GtRBTnode *root)
{
  unsigned long cnt = 0;
  GtRBTnode *p;

  if (root == NULL)
  {
    return;
  }
  /* gt_assert(!root->red); */
  /* root->red = false; */
  for (p = root->left; p != NULL; p = p->left)
  {
    if (!p->red)
    {
      cnt++;
    }
  }
  check_tree_recurse (root, 0, cnt);
}

#else

#define CHECK_TREE(a)                  /* Nothing */

#endif /* NDEBUG */

/**
 * Possibly "split" a node with two red successors,
 * and/or fix up two red edges
 * in a row.  rootp is a pointer to the lowest node we visited, PARENTP and
 * GPARENTP pointers to its parent/grandparent.  P_R and GP_R contain the
 * comparison values that determined which way was taken in the tree to reach
 * rootp.  MODE is 1 if we need not do the split, but must check for two red
 * edges between GPARENTP and rootp.
 */
static void maybe_split_for_insert(GtRBTnode **rootp,
                                   GtRBTnode **parentp,
                                   GtRBTnode **gparentp,
                                   int p_r,
                                   int gp_r,
                                   unsigned long mode)
{
  GtRBTnode *root = *rootp,
          **rp,
          **lp;

  rp = &(*rootp)->right;
  lp = &(*rootp)->left;

  /*
     See if we have to split this node (both successors red).
   */
  if (mode == 1UL
      || ((*rp) != NULL && (*lp) != NULL && (*rp)->red
          && (*lp)->red))
  {
    /*
       This node becomes red, its successors black.
     */
    root->red = true;
    if (*rp != NULL)
    {
      (*rp)->red = false;
    }
    if (*lp != NULL)
    {
      (*lp)->red = false;
    }

    /*
     * If the parent of this node is also red, we have to do rotations.
     */
    if (parentp != NULL && (*parentp)->red)
    {
      GtRBTnode *gp = *gparentp,
        *p = *parentp;

      /*
       * There are two main cases: 1. The edge types (left or right) of the two
       * red edges differ. 2. Both red edges are of the same type. There exist
       * two symmetries of each case, so there is a total of 4 cases.
       */
      if ((p_r > 0) != (gp_r > 0))
      {
        /*
         * Put the child at the top of the tree, with its parent and
         * grandparent as successors.
         */
        p->red = true;
        gp->red = true;
        root->red = false;
        if (p_r < 0)
        {
          /*
             Child is left of parent.
           */
          p->left = *rp;
          *rp = p;
          gp->right = *lp;
          *lp = gp;
        }
        else
        {
          /*
             Child is right of parent.
           */
          p->right = *lp;
          *lp = p;
          gp->left = *rp;
          *rp = gp;
        }
        *gparentp = root;
      }
      else
      {
        *gparentp = *parentp;
        /*
         * Parent becomes the top of the tree, grandparent and child are its
         * successors.
         */
        p->red = false;
        gp->red = true;
        if (p_r < 0)
        {
          /*
             Left edges.
           */
          gp->left = p->right;
          p->right = gp;
        }
        else
        {
          /*
             Right edges.
           */
          gp->right = p->left;
          p->left = gp;
        }
      }
    }
  }
}

/**
 * Find or insert datum with given key into search tree. rootp
 * is the address of tree root, cmpfun the ordering function.
 */
GtKeytype gt_rbt_search(const GtKeytype key,
                        bool *nodecreated,
                        GtRBTnode **rootp,
                        GtDictcomparefunction cmpfun,
                        void *cmpinfo)
{
  GtRBTnode *newnode,
          **parentp = NULL,
          **gparentp = NULL,
          **nextp;
          int r = 0,
              p_r = 0,
              gp_r = 0;   /* No they might not, Mr Compiler. */

  if (rootp == NULL)
  {
    *nodecreated = false;
    return GtKeytypeerror;
  }

  /*
     This saves some additional tests below.
   */
  if (*rootp != NULL)
  {
    (*rootp)->red = false;
  }
  CHECK_TREE (*rootp);
  nextp = rootp;
  while (*nextp != NULL)
  {
    GtRBTnode *root = *rootp;

    r = cmpfun (key, root->key, cmpinfo);
    if (r == 0)
    {
      *nodecreated = false;
      return root->key;
    }

    maybe_split_for_insert (rootp, parentp, gparentp, p_r, gp_r, 0);
    /*
     * If that did any rotations, parentp and gparentp are now garbage. That
     * doesn't matter, because the values they contain are never used again in
     * that case.
     */

    nextp = r < 0 ? &root->left : &root->right;
    if (*nextp == NULL)
    {
      break;
    }

    gparentp = parentp;
    parentp = rootp;
    rootp = nextp;

    gp_r = p_r;
    p_r = r;
  }

  newnode = (GtRBTnode *) gt_malloc (sizeof (GtRBTnode));
  *nextp = newnode;                    /* link new node to old */
  newnode->key = key;                  /* initialize new node */
  newnode->red = true;
  newnode->left = NULL;
  newnode->right = NULL;
  if (nextp != rootp)
  {
    /*
       There may be two red edges in a row now, which we must
       avoid by rotating the tree.
     */
    maybe_split_for_insert (nextp, rootp, parentp, r, p_r, 1UL);
  }
  *nodecreated = true;
  return newnode->key;
}

/**
 * Find datum in search tree. KEY is the key to be located, rootp is the
 * address of tree root, cmpfun the ordering function.
 */
GtKeytype gt_rbt_find(const GtKeytype key,
                      const GtRBTnode *root,
                      GtDictcomparefunction cmpfun,
                      void *cmpinfo)
{
  if (root == NULL)
  {
    return GtKeytypeerror;
  }
  CHECK_TREE (root);
  while (root != NULL)
  {
    int r;

    r = cmpfun (key, root->key, cmpinfo);
    if (r == 0)
    {
      return root->key;
    }
    if (r < 0)
    {
      root = root->left;
    }
    else
    {
      root = root->right;
    }
  }
  return GtKeytypeerror;
}

typedef GtRBTnode ** GtRBTnodeptrptr;

GT_STACK_DECLARESTRUCT(GtRBTnodeptrptr,32UL);

/**
 * Delete node with given key. rootp is the
 * address of the root of tree, cmpfun the comparison function.
 */
int gt_rbt_delete(const GtKeytype key,
                  GtRBTnode **rootp,
                  GtDictcomparefunction cmpfun,
                  void *cmpinfo)
{
  GtRBTnode *p,
            *q,
            *r,
            *root,
            *unchained;
  int cmp;
  GtStackGtRBTnodeptrptr nodestack;

  GT_STACK_INIT(&nodestack,16UL);
  p = *rootp;
  if (p == NULL)
  {
    return -1;
  }
  CHECK_TREE(p);
  while ((cmp = cmpfun (key, (*rootp)->key, cmpinfo)) != 0)
  {
    GT_STACK_PUSH(&nodestack,rootp);
    p = *rootp;
    if (cmp < 0)
    {
      rootp = &(*rootp)->left;
    }
    else
    {
      rootp = &(*rootp)->right;
    }
    if (*rootp == NULL)
    {
      GT_STACK_DELETE(&nodestack);
      return -1;
    }
  }

  /*
     We don't unchain the node we want to delete.
     Instead, we overwrite it with
     its successor and unchain the successor.  If there is no successor, we
     really unchain the node to be deleted.
   */

  root = *rootp;

  r = root->right;
  q = root->left;

  if (q == NULL || r == NULL)
  {
    unchained = root;
  }
  else
  {
    GtRBTnode **parent = rootp,
      **up = &root->right;

    for (;;)
    {
      GT_STACK_PUSH(&nodestack,parent);
      parent = up;
      if ((*up)->left == NULL)
      {
        break;
      }
      up = &(*up)->left;
    }
    unchained = *up;
  }

  /*
   * We know that either the left or right successor of UNCHAINED is NULL. R
   * becomes the other one, it is chained into the parent of UNCHAINED.
   */
  r = unchained->left;
  if (r == NULL)
  {
    r = unchained->right;
  }
  if (GT_STACK_ISEMPTY(&nodestack))
  {
    *rootp = r;
  }
  else
  {
    q = *(GT_STACK_TOP(&nodestack));
    if (unchained == q->right)
    {
      q->right = r;
    }
    else
    {
      q->left = r;
    }
  }

  if (unchained != root)
  {
    root->key = unchained->key;
  }
  if (!unchained->red)
  {
    /*
       Now we lost a black edge, which means that the number of
       black edges on every path is no longer constant.
       We must balance the tree.

       nodestack now contains all parents of R.
       R is likely to be NULL in the
       first iteration.

       NULL nodes are considered black throughout - this is necessary for
       correctness.
     */
    while (!GT_STACK_ISEMPTY(&nodestack) && (r == NULL || !r->red))
    {
      GtRBTnode **pp = GT_STACK_TOP(&nodestack);

      p = *pp;
      /*
         Two symmetric cases.
       */
      if (r == p->left)
      {
        /*
         * Q is R's brother, P is R's parent.  The subtree with root R has one
         * black edge less than the subtree with root Q.
         */
        q = p->right;
        if (q != NULL && q->red)
        {
          /*
           * If Q is red, we know that P is black. We rotate P left so that Q
           * becomes the top node in the tree, with P below it.  P is colored
           * red, Q is colored black. This action does not change the black
           * edge count for any leaf in the tree, but we will be able to
           * recognize one of the following situations, which all require that
           * Q is black.
           */
          q->red = false;
          p->red = true;
          /*
             Left rotate p.
           */
          p->right = q->left;
          q->left = p;
          *pp = q;
          /*
           * Make sure pp is right if the case below tries to use it.
           */
          pp = &q->left;
          GT_STACK_PUSH(&nodestack,pp); /* this has been added by S.K. */
          q = p->right;
        }
        gt_assert(q != NULL);
        /*
         * We know that Q can't be NULL here.  We also know that Q is black.
         */
        if ((q->left == NULL || !q->left->red)
            && (q->right == NULL || !q->right->red))
        {
          /*
           * Q has two black successors.  We can simply color Q red. The whole
           * subtree with root P is now missing one black edge. Note that this
           * action can temporarily make the tree invalid (if P is red).  But
           * we will exit the loop in that case and set P black, which both
           * makes the tree valid and also makes the black edge count come out
           * right.  If P is black, we are at least one step closer to the root
           * and we'll try again the next iteration.
           */
          q->red = true;
          r = p;
        }
        else
        {
          /*
           * Q is black, one of Q's successors is red.  We can repair the tree
           * with one operation and will exit the loop afterwards.
           */
          if (q->right == NULL || !q->right->red)
          {
            /*
             * The left one is red.  We perform the same action as in
             * maybe_split_for_insert where two red edges are adjacent but
             * point in different directions: Q's left successor (let's call it
             * Q2) becomes the top of the subtree we are looking at, its parent
             * (Q) and grandparent (P) become its successors. The former
             * successors of Q2 are placed below P and Q. P becomes black, and
             * Q2 gets the color that P had. This changes the black edge count
             * only for node R and its successors.
             */
            GtRBTnode *q2 = q->left;

            q2->red = p->red;
            p->right = q2->left;

            q->left = q2->right;
            q2->right = q;
            q2->left = p;
            *pp = q2;
            p->red = false;
          }
          else
          {
            /*
             * It's the right one.  Rotate P left. P becomes black, and Q gets
             * the color that P had.  Q's right successor also becomes black.
             * This changes the black edge count only for node R and its
             * successors.
             */
            q->red = p->red;
            p->red = false;

            q->right->red = false;

            /*
               left rotate p
             */
            p->right = q->left;
            q->left = p;
            *pp = q;
          }

          /*
             We're done.
           */
          GT_STACK_MAKEALMOSTEMPTY(&nodestack);
          r = NULL;
        }
      }
      else
      {
        /*
           Comments: see above.
         */
        q = p->left;
        if (q != NULL && q->red)
        {
          q->red = false;
          p->red = true;
          p->left = q->right;
          q->right = p;
          *pp = q;
          pp = &q->right;
          GT_STACK_PUSH(&nodestack,pp);   /* this has been added by S.K. */
          q = p->left;
        }
        gt_assert(q != NULL);
        if ((q->right == NULL || !q->right->red)
            && (q->left == NULL || !q->left->red))
        {
          q->red = true;
          r = p;
        }
        else
        {
          if (q->left == NULL || !q->left->red)
          {
            GtRBTnode *q2 = q->right;

            q2->red = p->red;
            p->left = q2->right;
            q->right = q2->left;
            q2->left = q;
            q2->right = p;
            *pp = q2;
            p->red = false;
          }
          else
          {
            q->red = p->red;
            p->red = false;
            q->left->red = false;
            p->left = q->right;
            q->right = p;
            *pp = q;
          }
          GT_STACK_MAKEALMOSTEMPTY(&nodestack);
          r = NULL;
        }
      }
      GT_STACK_DECREMENTTOP(&nodestack);
    }
    if (r != NULL)
    {
      r->red = false;
    }
  }
  gt_free (unchained);
  GT_STACK_DELETE(&nodestack);
  return 0;
}

/*
 Walk the nodes of a tree. root is the root of the tree to be walked,
 action the function to be called at each node. LEVEL is the level of
 root in the whole tree.
 */

static int mytreerecurse(const GtRBTnode *root,
                         GtDictaction action,
                         unsigned long level,
                         void *actinfo)
{
  if (root->left == NULL && root->right == NULL)
  {
    if (action (root->key, GT_RBT_LEAF, level, actinfo) != 0)
    {
      return -1;
    }
  }
  else
  {
    if (action (root->key, GT_RBT_PREORDER, level, actinfo) != 0)
    {
      return -2;
    }
    if (root->left != NULL)
    {
      if (mytreerecurse (root->left, action, level + 1, actinfo) != 0)
      {
        return -3;
      }
    }
    if (action (root->key, GT_RBT_POSTORDER, level, actinfo) != 0)
    {
      return -4;
    }
    if (root->right != NULL)
    {
      if (mytreerecurse (root->right, action, level + 1, actinfo) != 0)
      {
        return -5;
      }
    }
    if (action (root->key, GT_RBT_ENDORDER, level, actinfo) != 0)
    {
      return -6;
    }
  }
  return 0;
}

static int mytreerecursewithstop (const GtRBTnode *root,
                                  GtDictaction action,
                                  unsigned long level,
                                  void *actinfo)
{
  int retcode;

  if (root->left == NULL && root->right == NULL)
  {
    retcode = action (root->key, GT_RBT_LEAF, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  } else
  {
    retcode = action (root->key, GT_RBT_PREORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->left != NULL)
    {
      retcode = mytreerecursewithstop (root->left, action,
                                       level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action (root->key, GT_RBT_POSTORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
    if (root->right != NULL)
    {
      retcode = mytreerecursewithstop (root->right, action,
                                       level + 1, actinfo);
      RBT_CHECK_RETURN_CODE;
    }
    retcode = action (root->key, GT_RBT_ENDORDER, level, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

static int mytreerecursereverseorder (const GtRBTnode *root,
                                      GtDictaction action,
                                      unsigned long level,
                                      void *actinfo)
{
  if (root->left == NULL && root->right == NULL)
  {
    if (action (root->key, GT_RBT_LEAF, level, actinfo) != 0)
    {
      return -1;
    }
  } else
  {
    if (action (root->key, GT_RBT_PREORDER, level, actinfo) != 0)
    {
      return -2;
    }
    if (root->right != NULL)
    {
      if (mytreerecursereverseorder (root->right, action, level + 1,
                                     actinfo) != 0)
      {
        return -3;
      }
    }
    if (action (root->key, GT_RBT_POSTORDER, level, actinfo) != 0)
    {
      return -4;
    }
    if (root->left != NULL)
    {
      if (mytreerecursereverseorder
          (root->left, action, level + 1,
           actinfo) != 0)
      {
        return -5;
      }
    }
    if (action (root->key, GT_RBT_ENDORDER, level, actinfo) != 0)
    {
      return -6;
    }
  }
  return 0;
}

/*
 Walk the nodes of a tree. root is the root of the tree to be walked, action
 the function to be called at each node.
 */

int gt_rbt_walk (const GtRBTnode *root,GtDictaction action,void *actinfo)
{
  CHECK_TREE(root);
  if (root != NULL && action != NULL)
  {
    if (mytreerecurse (root, action, 0, actinfo) != 0)
    {
      return -1;
    }
  }
  return 0;
}

int gt_rbt_walkwithstop (const GtRBTnode *root, GtDictaction action,
                         void *actinfo)
{
  CHECK_TREE(root);
  if (root != NULL && action != NULL)
  {
    int retcode = mytreerecursewithstop (root, action, 0, actinfo);
    RBT_CHECK_RETURN_CODE;
  }
  return 0;
}

int gt_rbt_walkreverseorder(const GtRBTnode *root, GtDictaction action,
                            void *actinfo)
{
  CHECK_TREE(root);
  if (root != NULL && action != NULL)
  {
    if (mytreerecursereverseorder (root, action, 0, actinfo) != 0)
    {
      return -1;
    }
  }
  return 0;
}

/**
 find minimum key
 */

GtKeytype gt_rbt_minimumkey (const GtRBTnode *root)
{
  if (root == NULL)
  {
    return GtKeytypeerror;
  }
  while (root->left != NULL)
  {
    root = root->left;
  }
  return root->key;
}

/*
 find maximum key
 */

GtKeytype gt_rbt_maximumkey (const GtRBTnode *root)
{
  if (root == NULL)
  {
    return GtKeytypeerror;
  }
  while (root->right != NULL)
  {
    root = root->right;
  }
  return root->key;
}

void gt_rbt_treeshape (const GtRBTnode *root,unsigned long level)
{
  printf ("visit node %lu at level %lu\n",
          (unsigned long) *((unsigned long *) root->key),
          level);
  if (root->left != NULL)
  {
    printf ("visit left branch\n");
    gt_rbt_treeshape (root->left, level + 1);
  }
  if (root->right != NULL)
  {
    printf ("visit right branch\n");
    gt_rbt_treeshape (root->right, level + 1);
  }
  printf ("bracktrack\n");
}

/**
 * find largest element strictly smaller than key
 */
GtKeytype gt_rbt_previouskey(const GtKeytype key,
                             const GtRBTnode *root,
                             GtDictcomparefunction cmpfun,
                             void *cmpinfo)
{
  int cmp;
  const GtRBTnode *current = root,
                *found = NULL;

  while (current != NULL)
  {
    cmp = cmpfun (key, current->key, cmpinfo);
    if (cmp == 0)
    {
      if (current->left == NULL)
      {
        if (found == NULL)
        {
          return GtKeytypeerror;
        }
        return found->key;
      }
      return gt_rbt_maximumkey (current->left);
    } else
    {
      if (cmp < 0)
      {
        current = current->left;
      } else
      {
        found = current;
        current = current->right;
      }
    }
  }
  if (found == NULL)
  {
    return GtKeytypeerror;
  }
  return found->key;
}

/**
 * find largest element smaller than or equal to the key
 */
GtKeytype gt_rbt_previousequalkey(const GtKeytype key,
                                  const GtRBTnode *root,
                                  GtDictcomparefunction cmpfun,
                                  void *cmpinfo)
{
  int cmp;
  const GtRBTnode *current = root,
                *found = NULL;

  while (current != NULL)
  {
    cmp = cmpfun (key, current->key, cmpinfo);
    if (cmp == 0)
    {
      return current->key;
    } else
    {
      if (cmp < 0)
      {
        current = current->left;
      } else
      {
        found = current;
        current = current->right;
      }
    }
  }
  if (found == NULL)
  {
    return GtKeytypeerror;
  }
  return found->key;
}

/**
 * find smallest element strictly larger than key
 */
GtKeytype gt_rbt_nextkey(const GtKeytype key,
                         const GtRBTnode *root,
                         GtDictcomparefunction cmpfun,
                         void *cmpinfo)
{
  int cmp;
  const GtRBTnode *current = root,
                *found = NULL;

  while (current != NULL)
  {
    cmp = cmpfun (key, current->key, cmpinfo);
    if (cmp == 0)
    {
      if (current->right == NULL)
      {
        if (found == NULL)
        {
          return GtKeytypeerror;
        }
        return found->key;
      }
      return gt_rbt_minimumkey (current->right);
    } else
    {
      if (cmp < 0)
      {
        found = current;
        current = current->left;
      } else
      {
        current = current->right;
      }
    }
  }
  if (found == NULL)
  {
    return GtKeytypeerror;
  }
  return found->key;
}

/**
 * Lorem ipsum dolor sit amet
 */
GtKeytype gt_rbt_nextequalkey(const GtKeytype key,
                              const GtRBTnode *root,
                              GtDictcomparefunction cmpfun,
                              void *cmpinfo)
{
  int cmp;
  const GtRBTnode *current = root,
                *found = NULL;

  while (current != NULL)
  {
    cmp = cmpfun (key, current->key, cmpinfo);
    if (cmp == 0)
    {
      return current->key;
    } else
    {
      if (cmp < 0)
      {
        found = current;
        current = current->left;
      }
      else
      {
        current = current->right;
      }
    }
  }
  if (found == NULL)
  {
    return GtKeytypeerror;
  }
  return found->key;
}

/**
 * The standardized functions miss an important functionality:
 * the tree cannot be removed easily. We provide a function to do this.
 * If the boolean freekey is true, then also the key field is freed.
 * This only works if key points to a malloced block.
 */
static void redblacktreedestroyrecurse(bool dofreekey,
                                       GtFreekeyfunction freekey,
                                       void *freeinfo,
                                       GtRBTnode *root)
{
  if (root->left != NULL)
  {
    redblacktreedestroyrecurse (dofreekey, freekey, freeinfo,
                                root->left);
  }
  if (root->right != NULL)
  {
    redblacktreedestroyrecurse (dofreekey, freekey, freeinfo,
                                root->right);
  }
  if (dofreekey)
  {
    if (freekey == NULL)
    {
      gt_free ((void *) root->key);
    }
    else
    {
      freekey (root->key, freeinfo);
    }
  }
  /*
     Free the node itself.
   */
  gt_free (root);
}

void gt_rbt_destroy(bool dofreekey,
                    GtFreekeyfunction freekey,
                    void *freeinfo,
                    GtRBTnode *root)
{
  CHECK_TREE(root);
  if (root != NULL)
  {
    redblacktreedestroyrecurse (dofreekey, freekey, freeinfo, root);
  }
}

GtKeytype gt_rbt_extractrootkey (const GtRBTnode *root)
{
  if (root == NULL)
  {
    return GtKeytypeerror;
  }
  return root->key;
}

static int rangetreerecurse(const GtRBTnode *root,
                            GtDictaction action,
                            unsigned long level,
                            void *actinfo,
                            GtComparewithkey greaterequalleft,
                            GtComparewithkey lowerequalright,
                            void *cmpinfo)
{
  if (root->left == NULL && root->right == NULL)
  {
    if (greaterequalleft (root->key, cmpinfo)
        && lowerequalright (root->key, cmpinfo))
    {
      if (action (root->key, GT_RBT_LEAF, level, actinfo) != 0)
      {
        return -1;
      }
    }
  }
  else
  {
    if (root->left != NULL)
    {
      if (greaterequalleft (root->key, cmpinfo))
      {
        if (rangetreerecurse (root->left, action, level + 1, actinfo,
                              greaterequalleft, lowerequalright, cmpinfo) != 0)
        {
          return -2;
        }
      }
    }
    if (greaterequalleft (root->key, cmpinfo)
        && lowerequalright (root->key, cmpinfo))
    {
      if (action (root->key, GT_RBT_POSTORDER, level, actinfo) != 0)
      {
        return -3;
      }
    }
    if (root->right != NULL)
    {
      if (lowerequalright (root->key, cmpinfo))
      {
        if (rangetreerecurse (root->right, action, level + 1, actinfo,
                              greaterequalleft, lowerequalright, cmpinfo) != 0)
        {
          return -4;
        }
      }
    }
  }
  return 0;
}

int gt_rbt_walkrange(const GtRBTnode *root,
                     GtDictaction action,
                     void *actinfo,
                     GtComparewithkey greaterequalleft,
                     GtComparewithkey lowerequalright,
                     void *cmpinfo)
{
  CHECK_TREE(root);
  if (root != NULL && action != NULL)
  {
    if (rangetreerecurse (root, action, 0, actinfo, greaterequalleft,
                          lowerequalright, cmpinfo) != 0)
    {
      return -1;
    }
  }
  return 0;
}

/**
 * This file contains unit tests for the red black tree datastructure
 */

#define SEED 0
#define PASSES 100
#define SIZE 100

typedef enum
{
  GT_RBT_ASCENDING,
  GT_RBT_DESCENDING,
  GT_RBT_RANDOMORDER
} SearchOrder;

typedef enum
{
  Build,
  Build_and_del,
  Delete,
  Find
} DoAction;

/* Set to 1 if a test is flunked.
 static signed long error = 0;
 */

/* The keys we add to the tree.  */
static unsigned long xtab[SIZE];

/*
 * Pointers into the key array, possibly permutated, to define an order for
 * insertion/removal.
 */
static unsigned long ytab[SIZE];

/* Flags set for each element visited during a tree walk.  */
static unsigned long ztab[SIZE];

/*
 * Depths for all the elements, to check that the depth is constant for all
 * three visits.
 */
static unsigned long depths[SIZE];

/* Maximum depth during a tree walk.  */
static unsigned long max_depth;

/* Compare two keys.  */
static int cmp_fn (const GtKeytype a,const GtKeytype b,GT_UNUSED void *info)
{
  unsigned long va, vb;

  va = *(unsigned long *) a;
  vb = *(unsigned long *) b;
  if (va < vb)
  {
    return -1;
  }
  if (va > vb)
  {
    return 1;
  }
  return 0;
}

/* Permute an array of integers.  */

static void permuteintarray (unsigned long *arr)
{
  unsigned long i, j, c;

  for (i = 0; i < (unsigned long) SIZE; ++i)
  {
    j = (unsigned long) (random () % SIZE);
    c = arr[i];
    arr[i] = arr[j];
    arr[j] = c;
  }
}

static int walk_action (const GtKeytype nodekey,
                        const GtRbtVisit which,
                        const unsigned long depth,
                        GT_UNUSED void *actinfo)
{
  unsigned long key = *(unsigned long *) nodekey;

  if (depth > max_depth)
  {
    max_depth = depth;
  }
  if (which == GT_RBT_LEAF || which == GT_RBT_PREORDER)
  {
    ++ztab[key];
    depths[key] = depth;
  } else
  {
    if (depths[key] != depth)
    {
      gt_xfputs("Depth for one element is not constant during tree walk.\n",
                stdout);
      return -1;
    }
  }
  return 0;
}

static int walk_tree (const void *root,unsigned long expected_count)
{
  unsigned long i;
  int error = 0;

  memset (ztab, 0, sizeof ztab);
  max_depth = 0;

  if (gt_rbt_walk (root, walk_action, NULL) != 0)
  {
    gt_xfputs("walk failed\n", stdout);
    error = 1;
  }
  for (i = 0; i < expected_count; ++i)
  {
    if (ztab[i] != 1UL)
    {
      gt_xfputs("Node was not visited.\n", stdout);
      error = 1;
    }
  }
  if (max_depth >
      (unsigned long) (log ((double) expected_count) * 2.0 +
                       2.0))
  {
    gt_xfputs("Depth too large during tree walk.\n", stdout);
    error = 1;
  }
  return error;
}

/* Perform an operation on a tree.  */
static bool mangle_tree (SearchOrder how,
                         DoAction what,
                         GtRBTnode **root,
                         unsigned long lag,
                         GtError *err)
{
  unsigned long i;
  bool nodecreated, haserr = false;

  if (how == GT_RBT_RANDOMORDER)
  {
    for (i = 0; i < (unsigned long) SIZE; ++i)
    {
      ytab[i] = i;
    }
    permuteintarray (ytab);
  }
  for (i = 0; i < SIZE + lag; ++i)
  {
    void *elem;
    unsigned long j, k;

    switch (how)
    {
      case GT_RBT_RANDOMORDER:
        if (i >= lag)
        {
          k = ytab[i - lag];
        } else
        {
          k = ytab[SIZE - i - 1 + lag];
        }
        j = ytab[i];
        break;

      case GT_RBT_ASCENDING:
        k = i - lag;
        j = i;
        break;

      case GT_RBT_DESCENDING:
        k = SIZE - i - 1 + lag;
        j = SIZE - i - 1;
        break;

      default:
        /*
         * This never should happen, but gcc isn't smart enough to recognize it.
         */
        abort ();
    }

    switch (what)
    {
      case Build_and_del:
      case Build:
        if (i < (unsigned long) SIZE)
        {
          if (gt_rbt_find (xtab + j, *root, cmp_fn, NULL) !=
              NULL)
          {
            gt_error_set(err,"Found element which is not in tree yet");
            haserr = true;
          }
          if (!haserr)
          {
            elem = gt_rbt_search (xtab + j, &nodecreated, root, cmp_fn, NULL);
            if (elem == NULL || gt_rbt_find (xtab + j, *root, cmp_fn,
                                                  NULL) == NULL)
            {
              gt_error_set(err,"Couldn't find element after it was added");
              haserr = true;
            }
          }
        }
        if (haserr || what == Build || i < lag)
        {
          break;
        }
        j = k;
        /*@fallthrough@*/
      case Delete:
        elem = gt_rbt_find (xtab + j, *root, cmp_fn, NULL);
        if (elem == NULL
            || gt_rbt_delete (xtab + j, root, cmp_fn,
                                   NULL) != 0)
        {
          gt_error_set(err,"Error deleting element");
          haserr = true;
        }
        break;

      case Find:
        if (gt_rbt_find (xtab + j, *root, cmp_fn, NULL) == NULL)
        {
          gt_error_set(err,"Couldn't find element after it was added");
          haserr = true;
        }
        break;
    }
  }
  return haserr;
}

#define MANGLECHECK(ORDER,MODE,LAG)\
        ensure (had_err,!mangle_tree (ORDER, MODE, &root, LAG,err))

#define WALKCHECK\
        ensure (had_err,!walk_tree (root, (unsigned long) SIZE))

/**
 * Unit test for the red black tree datastructure
 */

int gt_rbt_unit_test (GtError *err)
{
  int had_err = 0;
  GtRBTnode *root = NULL;
  unsigned long i, j;

  gt_error_check (err);
  for (i = 0; i < (unsigned long) SIZE; ++i)
  {
    xtab[i] = i;
  }
 /**
  * Do this loop several times to get different permutations for the random
  * case.
  */
  for (i = 0; i < (unsigned long) PASSES; ++i)
  {
    MANGLECHECK(GT_RBT_ASCENDING, Build, 0);
    MANGLECHECK(GT_RBT_ASCENDING, Find, 0);
    MANGLECHECK(GT_RBT_DESCENDING, Find, 0);
    MANGLECHECK(GT_RBT_RANDOMORDER, Find, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_ASCENDING, Delete, 0);

    MANGLECHECK (GT_RBT_ASCENDING, Build, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_DESCENDING, Delete, 0);

    MANGLECHECK (GT_RBT_ASCENDING, Build, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_RANDOMORDER, Delete, 0);

    MANGLECHECK (GT_RBT_DESCENDING, Build, 0);
    MANGLECHECK (GT_RBT_ASCENDING, Find, 0);
    MANGLECHECK (GT_RBT_DESCENDING, Find, 0);
    MANGLECHECK (GT_RBT_RANDOMORDER, Find, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_DESCENDING, Delete, 0);

    MANGLECHECK (GT_RBT_DESCENDING, Build, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_DESCENDING, Delete, 0);

    MANGLECHECK (GT_RBT_DESCENDING, Build, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_RANDOMORDER, Delete, 0);

    MANGLECHECK (GT_RBT_RANDOMORDER, Build, 0);
    MANGLECHECK (GT_RBT_ASCENDING, Find, 0);
    MANGLECHECK (GT_RBT_DESCENDING, Find, 0);
    MANGLECHECK (GT_RBT_RANDOMORDER, Find, 0);
    WALKCHECK;
    MANGLECHECK (GT_RBT_RANDOMORDER, Delete, 0);

    for (j = 1UL; j < (unsigned long) SIZE; j *= 2)
    {
      MANGLECHECK (GT_RBT_RANDOMORDER, Build_and_del, j);
    }
  }
  for (i = 1UL; i < (unsigned long) SIZE; i *= 2)
  {
    MANGLECHECK (GT_RBT_ASCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_DESCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_ASCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_DESCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_ASCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_DESCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_ASCENDING, Build_and_del, i);
    MANGLECHECK (GT_RBT_DESCENDING, Build_and_del, i);
  }
  return had_err;
}
