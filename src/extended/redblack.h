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
  preorder,
  postorder,
  endorder,
  leaf
} VISIT;

typedef void *Keytype;

#define Keytypeerror NULL

typedef struct RBTnode RBTnode;

typedef int (*Dictcomparefunction) (const Keytype,const Keytype, void *);
typedef void (*Dictshowelem) (const Keytype,void *);
typedef int (*Dictaction) (const Keytype,VISIT,unsigned long, void *);
typedef void (*Freekeyfunction) (const Keytype,void *);
typedef bool (*Comparewithkey) (const Keytype, void *);

Keytype rbt_search (const Keytype key,
                    bool *nodecreated,
                    RBTnode **rootp,
                    Dictcomparefunction cmpfun,
                    void *cmpinfo);

Keytype rbt_find (const Keytype key,
                  const RBTnode *root,
                  Dictcomparefunction cmpfun,
                  void *cmpinfo);

int rbt_delete (const Keytype key,
                RBTnode **rootp,
                Dictcomparefunction cmpfun,
                void *cmpinfo);

int rbt_walk (const RBTnode *root,Dictaction action,void *actinfo);

int rbt_walkwithstop (const RBTnode *root,Dictaction action,void *actinfo);

int rbt_walkreverseorder (const RBTnode *root,Dictaction action,void *actinfo);

Keytype rbt_minimumkey (const RBTnode *root);

Keytype rbt_maximumkey (const RBTnode *root);

void rbt_treeshape (const RBTnode *root,unsigned long level);

Keytype rbt_previouskey (const Keytype key,
                         const RBTnode *root,
                         Dictcomparefunction cmpfun,
                         void *cmpinfo);

Keytype rbt_previousequalkey (const Keytype key,
                              const RBTnode *root,
                              Dictcomparefunction cmpfun,
                              void *cmpinfo);

Keytype rbt_nextkey (const Keytype key,
                     const RBTnode *root,
                     Dictcomparefunction cmpfun,
                     void *cmpinfo);

Keytype rbt_nextequalkey (const Keytype key,
                          const RBTnode *root,
                          Dictcomparefunction cmpfun,
                          void *cmpinfo);

void rbt_destroy (bool dofreekey,
                  Freekeyfunction freekey,
                  void *freeinfo,
                  RBTnode *root);

Keytype rbt_extractrootkey (const RBTnode *root);

int rbt_walkrange (const RBTnode *root,
                   Dictaction action,
                   void *actinfo,
                   Comparewithkey greaterequalleft,
                   Comparewithkey lowerequalright,
                   void *cmpinfo);

int rbt_unit_test (GT_Error *err);

#endif
