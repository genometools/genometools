/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/ma.h"
#include "extended/maxcoordvalue.h"
#include "extended/linspace_management.h"

struct GtLinspaceManagement{
  void             *valueTabspace,
                   *rTabspace,
                   *crosspointTabspace;
  GtUword          ulen,
                   timesquarefactor;
  size_t           valueTabsize,
                   rTabsize,
                   crosspointTabsize,
                   spacepeak; /*sum of space in bytes*/
  GtMaxcoordvalue *maxscoordvaluespace;
};

GtLinspaceManagement* gt_linspace_management_new()
{
  GtLinspaceManagement *spacemanager;
  spacemanager = gt_malloc(sizeof(*spacemanager));
  spacemanager->valueTabspace = NULL;
  spacemanager->rTabspace = NULL;
  spacemanager->crosspointTabspace = NULL;
  spacemanager->maxscoordvaluespace = NULL;
  spacemanager->valueTabsize = 0;
  spacemanager->rTabsize = 0;
  spacemanager->crosspointTabsize = 0;
  spacemanager->timesquarefactor = 1;
  spacemanager->ulen = 0;
  spacemanager->spacepeak = 0;
  return spacemanager;
}

void gt_linspace_management_delete(GtLinspaceManagement *spacemanager)
{
  if (spacemanager != NULL)
  {
    gt_free(spacemanager->valueTabspace);
    gt_free(spacemanager->rTabspace);
    if (spacemanager->crosspointTabspace != NULL)
      gt_free(spacemanager->crosspointTabspace);
    gt_maxcoordvalue_delete(spacemanager->maxscoordvaluespace);
    gt_free(spacemanager);
  }
}

/* space in bytes */
size_t gt_linspace_management_get_spacepeak(const GtLinspaceManagement
                                                                  *spacemanager)
{
  gt_assert(spacemanager != NULL);
  return spacemanager->spacepeak;
}

/* resize space */
static void gt_linspace_management_check_generic(GtLinspaceManagement
                                                 *spacemanager,
                                                 GtUword ulen, GtUword vlen,
                                                 size_t valuesize,
                                                 size_t rtabsize,
                                                 size_t crosspointsize,
                                                 bool local)
{
  size_t space = 0, localspace = 0; /* space in bytes */
  /*if (spacemanager == NULL)
    spacemanager = gt_new_linspace_management();*/

  gt_assert(spacemanager != NULL);

  if (spacemanager->valueTabsize < (ulen+1)*valuesize)
  {
    spacemanager->valueTabspace = gt_realloc(spacemanager->valueTabspace,
                                            (ulen+1)*valuesize);
    spacemanager->valueTabsize = (ulen+1)*valuesize;
  }

  if (spacemanager->rTabsize < (ulen+1)*rtabsize)
  {
    spacemanager->rTabspace = gt_realloc(spacemanager->rTabspace,
                                        (ulen+1)*rtabsize);
    spacemanager->rTabsize = (ulen+1)*rtabsize;
  }
  if (spacemanager->crosspointTabsize < (vlen+1)*crosspointsize)
  {
    spacemanager->crosspointTabspace =
           gt_realloc(spacemanager->crosspointTabspace,(vlen+1)*crosspointsize);
    spacemanager->crosspointTabsize = (vlen+1)*crosspointsize;
  }
  if (local)
  {
    if (spacemanager->maxscoordvaluespace == NULL)
      spacemanager->maxscoordvaluespace = gt_maxcoordvalue_new();
    else
      gt_maxcoordvalue_reset(spacemanager->maxscoordvaluespace);
  }
  if (spacemanager->maxscoordvaluespace != NULL)
    localspace = 2 * sizeof (GtUwordPair) + sizeof (GtWord);
            /* = sizeof (GtMaxcoordvalue)*/

  /* determine space peak */
  space = spacemanager->valueTabsize + spacemanager->rTabsize +
          spacemanager->crosspointTabsize + localspace;

  if (space > spacemanager->spacepeak)
    spacemanager->spacepeak = space;
}

void gt_linspace_management_check(GtLinspaceManagement *spacemanager,
                                  GtUword ulen, GtUword vlen,
                                  size_t valuesize,
                                  size_t rtabsize,
                                  size_t crosspointsize)
{
  gt_linspace_management_check_generic(spacemanager,
                              ulen, vlen,
                              valuesize,
                              rtabsize,
                              crosspointsize,
                              false);
  spacemanager->ulen = ulen;
}

void  gt_linspace_management_check_local(GtLinspaceManagement *spacemanager,
                                         GtUword ulen, GtUword vlen,
                                         size_t valuesize, size_t rstabsize)
{
  gt_linspace_management_check_generic(spacemanager,
                                       ulen, vlen,
                                       valuesize,
                                       rstabsize,
                                       0, true);
  spacemanager->ulen = ulen;
}

static bool checksquare(GtLinspaceManagement *spacemanager,
                        GtUword ulen, GtUword vlen,
                        size_t valuesize, size_t rsize,
                        bool local)
{
  GtUword TSfactor;
  gt_assert(spacemanager);

  TSfactor = spacemanager->timesquarefactor;
  if ((ulen+1)*(vlen+1)*valuesize <= spacemanager->valueTabsize)
  {
    if (local)
      gt_maxcoordvalue_reset(spacemanager->maxscoordvaluespace);
    return true;
  }
  else if ((ulen+1)*(vlen+1) <= (spacemanager->ulen+1)*TSfactor)
  {
    if (!local)
    {
      gt_linspace_management_check_generic(spacemanager,(ulen+1)*(vlen+1)-1,
                                           vlen,valuesize,rsize,0,false);
    }
    else
    {
       gt_linspace_management_check_generic(spacemanager,(ulen+1)*(vlen+1)-1,
                                            vlen,valuesize,rsize,0,true);
    }
    return true;
  }
  return false;
}

bool gt_linspace_management_checksquare(GtLinspaceManagement *spacemanager,
                                        GtUword ulen,
                                        GtUword vlen,
                                        size_t valuesize,
                                        size_t rsize)
{
  return checksquare(spacemanager, ulen, vlen, valuesize, rsize,false);
}

bool gt_linspace_management_checksquare_local(
                                             GtLinspaceManagement *spacemanager,
                                             GtUword ulen,
                                             GtUword vlen,
                                             size_t valuesize,
                                             size_t rsize)
{
  return checksquare(spacemanager, ulen, vlen, valuesize, rsize,true);
}

void gt_linspace_management_set_ulen(GtLinspaceManagement *spacemanager,
                                     GtUword ulen)
{
  gt_assert(spacemanager != NULL);
  spacemanager->ulen = ulen;
}

void *gt_linspace_management_get_valueTabspace(const GtLinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->valueTabspace);
  return NULL;
}

void *gt_linspace_management_get_rTabspace(const GtLinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->rTabspace);
  return NULL;
}

void *gt_linspace_management_get_crosspointTabspace(const GtLinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->crosspointTabspace);
  return NULL;
}

size_t gt_linspace_management_get_valueTabsize(const GtLinspaceManagement
                                                                  *spacemanager)
{
  gt_assert(spacemanager != NULL);
  return spacemanager->valueTabsize;
}

/* space for GtMaxcoordvalue struct */
void *gt_linspace_management_get_maxspace(const GtLinspaceManagement
                                                                 *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->maxscoordvaluespace);
  return NULL;
}

static inline bool check(const GtLinspaceManagement *spacemanager,
                         GtUword ulen, GtUword vlen)
{
  return ((ulen+1)*(vlen+1)*sizeof(GtUword) <= spacemanager->valueTabsize);
}

GtUword **gt_linspace_management_change_to_square(GtLinspaceManagement
                                                 *spacemanager,
                                                 GtUword ulen, GtUword vlen)
{
  GtUword **E;
  GtUword idx;
  gt_assert(check(spacemanager, ulen, vlen));

  E = gt_linspace_management_get_rTabspace(spacemanager);
  *E = gt_linspace_management_get_valueTabspace(spacemanager);

  for (idx=1;idx<ulen+1;idx++)
    E[idx]=E[idx-1]+vlen+1;

  return E;
}

void gt_linspace_management_set_TSfactor(GtLinspaceManagement *spacemanager,
                                        GtUword timesquarefactor)
{
  gt_assert(spacemanager != NULL);
  spacemanager->timesquarefactor = timesquarefactor;
}
