/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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

#include "extended/linspaceManagement.h"
struct LinspaceManagement{
  void             *valueTabspace,
                   *rTabspace,
                   *crosspointTabspace;
  GtUword          valuesize,
                   rsize,
                   crosspointsize,
                   timesquarefactor,
                   ulen;
  size_t           spacepeak; /*space in bytes*/
  Gtmaxcoordvalue *maxscoordvaluespace;
};

LinspaceManagement* gt_linspaceManagement_new()
{
  LinspaceManagement *spacemanager;
  spacemanager = gt_malloc(sizeof(*spacemanager));
  spacemanager->valueTabspace = NULL;
  spacemanager->rTabspace = NULL;
  spacemanager->crosspointTabspace = NULL;
  spacemanager->maxscoordvaluespace = NULL;
  spacemanager->valuesize = 0;
  spacemanager->rsize = 0;
  spacemanager->crosspointsize = 0;
  spacemanager->timesquarefactor = 1;
  spacemanager->ulen = 0;
  spacemanager->spacepeak = 0;
  return spacemanager;
}

void gt_linspaceManagement_delete(LinspaceManagement *spacemanager)
{
  if (spacemanager != NULL)
  {
    gt_free(spacemanager->valueTabspace);
    gt_free(spacemanager->rTabspace);
    if (spacemanager->crosspointTabspace != NULL)
      gt_free(spacemanager->crosspointTabspace);
    gt_max_delete(spacemanager->maxscoordvaluespace);
    gt_free(spacemanager);
  }
}

/* space in bytes */
size_t gt_linspaceManagement_get_spacepeak(const LinspaceManagement
                                                                  *spacemanager)
{
  gt_assert(spacemanager != NULL);
  return spacemanager->spacepeak;
}

/* resize space */
static void gt_linspaceManagement_check_generic(LinspaceManagement *spacemanager,
                                                GtUword ulen, GtUword vlen,
                                                size_t valuesize,
                                                size_t rtabsize,
                                                size_t crosspointsize,
                                                bool local)
{
  size_t space = 0, localspace = 0; /* space in bytes */
  /*if (spacemanager == NULL)
    spacemanager = gt_new_linspaceManagement();*/

  gt_assert(spacemanager != NULL);
  gt_assert(spacemanager->valuesize == spacemanager->rsize);

  if (spacemanager->valuesize < ulen+1)
  {
    spacemanager->valueTabspace = gt_realloc(spacemanager->valueTabspace,
                                            (ulen+1)*valuesize);
    spacemanager->valuesize = ulen+1;

    spacemanager->rTabspace = gt_realloc(spacemanager->rTabspace,
                                        (ulen+1)*rtabsize);
    spacemanager->rsize = ulen+1;
  }
  if (crosspointsize && spacemanager->crosspointsize < vlen+1)
  {
    spacemanager->crosspointTabspace =
           gt_realloc(spacemanager->crosspointTabspace,(vlen+1)*crosspointsize);
    spacemanager->crosspointsize = vlen+1;
  }
  if (local)
  {
    if (spacemanager->maxscoordvaluespace == NULL)
      spacemanager->maxscoordvaluespace = gt_max_new();
    else
      gt_max_reset(spacemanager->maxscoordvaluespace);

    localspace = 2 * sizeof (GtUwordPair)+ sizeof (GtWord);
            /* = sizeof (Gtmaxcoordvalue)*/
  }

  /* determine space peak */
  space = (ulen+1) * valuesize + (ulen +1) * rtabsize
        + (vlen+1) * crosspointsize + localspace;
  if (space > spacemanager->spacepeak)
    spacemanager->spacepeak = space;
}

void gt_linspaceManagement_check(LinspaceManagement *spacemanager,
                                 GtUword ulen, GtUword vlen,
                                 size_t valuesize,
                                 size_t rtabsize,
                                 size_t crosspointsize)
{
  gt_linspaceManagement_check_generic(spacemanager,
                              ulen, vlen,
                              valuesize,
                              rtabsize,
                              crosspointsize,
                              false);
  spacemanager->ulen = ulen;
}

void  gt_linspaceManagement_check_local(LinspaceManagement *spacemanager,
                                        GtUword ulen, GtUword vlen,
                                        size_t valuesize, size_t rstabsize)
{
  gt_linspaceManagement_check_generic(spacemanager,
                              ulen, vlen,
                              valuesize,
                              rstabsize,
                              0, true);
  spacemanager->ulen = ulen;
}

static bool checksquare(LinspaceManagement *spacemanager,
                        GtUword ulen, GtUword vlen,
                        size_t valuesize, size_t rsize,
                        bool local)
{
  GtUword TSfactor;
  gt_assert(spacemanager);

  TSfactor = spacemanager->timesquarefactor;
  if ((ulen+1)*(vlen+1) <= spacemanager->valuesize)
    return true;
  else if (ulen == 1|| vlen == 1 ||
          ((ulen+1)*(vlen+1) <= (spacemanager->ulen+1)*TSfactor))
  {
    if (!local)
    {
      gt_linspaceManagement_check_generic(spacemanager,(ulen+1)*(vlen+1)-1,
                                          vlen,valuesize,rsize,0,false);
    }
    else
    {
       gt_linspaceManagement_check_generic(spacemanager,(ulen+1)*(vlen+1)-1,
                                           vlen,valuesize,rsize,0,true);
    }
    return true;
  }
  return false;
}

bool gt_linspaceManagement_checksquare(LinspaceManagement *spacemanager,
                                       GtUword ulen,
                                       GtUword vlen,
                                       size_t valuesize,
                                       size_t rsize)
{
  return checksquare(spacemanager, ulen, vlen, valuesize, rsize,false);
}

bool gt_linspaceManagement_checksquare_local(LinspaceManagement *spacemanager,
                                             GtUword ulen,
                                             GtUword vlen,
                                             size_t valuesize,
                                             size_t rsize)
{
  return checksquare(spacemanager, ulen, vlen, valuesize, rsize,true);
}

void gt_linspaceManagement_set_ulen(LinspaceManagement *spacemanager,
                                    GtUword ulen)
{
  gt_assert(spacemanager != NULL);
  spacemanager->ulen = ulen;
}

void *gt_linspaceManagement_get_valueTabspace(const LinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->valueTabspace);
  return NULL;
}

void *gt_linspaceManagement_get_rTabspace(const LinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->rTabspace);
  return NULL;
}

void *gt_linspaceManagement_get_crosspointTabspace(const LinspaceManagement
                                                                  *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->crosspointTabspace);
  return NULL;
}

GtUword gt_linspaceManagement_get_valuesize(const LinspaceManagement
                                                                  *spacemanager)
{
  gt_assert(spacemanager != NULL);
  return spacemanager->valuesize;
}

/* space for Gtmaxcoordvalue struct */
void *gt_linspaceManagement_get_maxspace(const LinspaceManagement *spacemanager)
{
  if (spacemanager != NULL)
    return (spacemanager->maxscoordvaluespace);
  return NULL;
}

static inline bool check(const LinspaceManagement *spacemanager,
                         GtUword ulen, GtUword vlen)
{
  return ((ulen+1)*(vlen+1) <= spacemanager->valuesize);
}

GtUword **gt_linspaceManagement_change_to_square(LinspaceManagement *spacemanager,
                                                 GtUword ulen, GtUword vlen)
{
  GtUword **E;
  GtUword idx;
  gt_assert(check(spacemanager, ulen, vlen));

  E=gt_linspaceManagement_get_valueTabspace(spacemanager);
  *E=gt_linspaceManagement_get_rTabspace(spacemanager);

  for (idx=1;idx<ulen+1;idx++)
    E[idx]=E[idx-1]+vlen+1;

  return E;
}

void gt_linspaceManagement_set_TSfactor(LinspaceManagement *spacemanager,
                                        GtUword timesquarefactor)
{
  gt_assert(spacemanager != NULL);
  spacemanager->timesquarefactor = timesquarefactor;
}

GtUchar* sequence_to_lower_case(const GtUchar *seq, GtUword len)
{
  GtUword i;
  GtUchar *low_seq;

  low_seq = gt_malloc(sizeof(*low_seq)*(len+1));
  for (i = 0; i < len; i++)
    low_seq[i] = tolower((int)seq[i]);
  low_seq[i] = '\0';

  return low_seq;
}

inline GtWord add_safe(GtWord val1, GtWord val2, GtWord exception)
{
  return (val1 != exception) ? val1 + val2 : exception;
}

inline GtWord add_safe_max(GtWord val1, GtWord val2)
{
  return add_safe(val1,val2,GT_WORD_MAX);
}

inline GtWord add_safe_min(GtWord val1, GtWord val2)
{
  return add_safe(val1,val2,GT_WORD_MIN);
}

inline GtUword add_usafe(GtUword val1, GtUword val2, GtUword exception)
{
  return (val1 != exception) ? val1 + val2 : exception;
}

inline GtUword add_safe_umax(GtUword val1, GtUword val2)
{
  return add_usafe(val1, val2, GT_UWORD_MAX);
}
