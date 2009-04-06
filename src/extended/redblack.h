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

#ifndef REDBLACK_H
#define REDBLACK_H

#include "core/error.h"

typedef enum
{
  GT_RBT_PREORDER,
  GT_RBT_POSTORDER,
  GT_RBT_ENDORDER,
  GT_RBT_LEAF
} GtRbtVisit;

typedef void *GtKeytype;

#define GtKeytypeerror NULL

typedef struct GtRBTnode GtRBTnode;

typedef int (*Dictcomparefunction) (const GtKeytype,const GtKeytype, void *);
typedef void (*Dictshowelem) (const GtKeytype,void *);
typedef int (*Dictaction) (const GtKeytype,GtRbtVisit,unsigned long, void *);
typedef void (*Freekeyfunction) (const GtKeytype,void *);
typedef bool (*Comparewithkey) (const GtKeytype, void *);

GtKeytype gt_rbt_search (const GtKeytype key,
                    bool *nodecreated,
                    GtRBTnode **rootp,
                    Dictcomparefunction cmpfun,
                    void *cmpinfo);

GtKeytype gt_rbt_find (const GtKeytype key,
                  const GtRBTnode *root,
                  Dictcomparefunction cmpfun,
                  void *cmpinfo);

int gt_rbt_delete (const GtKeytype key,
                GtRBTnode **rootp,
                Dictcomparefunction cmpfun,
                void *cmpinfo);

int gt_rbt_walk (const GtRBTnode *root,Dictaction action,void *actinfo);

int gt_rbt_walkwithstop (const GtRBTnode *root,Dictaction action,void *actinfo);

int gt_rbt_walkreverseorder (const GtRBTnode *root,Dictaction action,
                             void *actinfo);

GtKeytype gt_rbt_minimumkey (const GtRBTnode *root);

GtKeytype gt_rbt_maximumkey (const GtRBTnode *root);

void gt_rbt_treeshape (const GtRBTnode *root,unsigned long level);

GtKeytype gt_rbt_previouskey (const GtKeytype key,
                         const GtRBTnode *root,
                         Dictcomparefunction cmpfun,
                         void *cmpinfo);

GtKeytype gt_rbt_previousequalkey (const GtKeytype key,
                              const GtRBTnode *root,
                              Dictcomparefunction cmpfun,
                              void *cmpinfo);

GtKeytype gt_rbt_nextkey (const GtKeytype key,
                     const GtRBTnode *root,
                     Dictcomparefunction cmpfun,
                     void *cmpinfo);

GtKeytype gt_rbt_nextequalkey (const GtKeytype key,
                          const GtRBTnode *root,
                          Dictcomparefunction cmpfun,
                          void *cmpinfo);

void gt_rbt_destroy (bool dofreekey,
                  Freekeyfunction freekey,
                  void *freeinfo,
                  GtRBTnode *root);

GtKeytype gt_rbt_extractrootkey (const GtRBTnode *root);

int gt_rbt_walkrange (const GtRBTnode *root,
                   Dictaction action,
                   void *actinfo,
                   Comparewithkey greaterequalleft,
                   Comparewithkey lowerequalright,
                   void *cmpinfo);

int gt_rbt_unit_test (GtError *err);

#endif
