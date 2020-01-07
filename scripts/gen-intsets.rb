#!/usr/bin/env ruby
# Copyright (c) 2014-2016 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
# Copyright (c) 2014-2016 Center for Bioinformatics, University of Hamburg

# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.

# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

require 'optparse'
require 'ostruct'
require 'erb'

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'codegen_module'

#begin=c
CRIGHT = <<-CRIGHT
/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>\
<% if name!="Dirk Willrodt" %>
  Copyright (c) <%=year%> <%=name%> <<%=email%>><% end %>
  Copyright (c) 2014<% if not year==2014 %>-<%=year%><% end %> Center for \
Bioinformatics, University of Hamburg

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
  THIS FILE IS GENERATED by\n  #{$0}.
  DO NOT EDIT.
*/

CRIGHT
#end=c

#begin=c
IMPL = <<-IMPL
#include <inttypes.h>
#include <limits.h>

#include "core/assert_api.h"
#include "core/divmodmul_api.h"
#include "core/ensure_api.h"
#include "core/intbits.h"
#include "core/ma_api.h"
#include "core/mathsupport_api.h"
#include "core/unused_api.h"
#include "extended/intset_<%=bits%>.h"
#include "extended/io_function_pointers.h"

#define gt_intset_<%=bits%>_cast(cvar) \\
        gt_intset_cast(gt_intset_<%=bits%>_class(), cvar)

#define GT_ELEM2SECTION_M(X) GT_ELEM2SECTION(X, members->logsectionsize)

#define GT_INTSET_<%=bits%>_TYPE ((GtUword) <%=bits%>)

#define gt_intset_<%=bits%>_io_one(element) \\
        io_func(&element, sizeof (element), (size_t) 1, fp, err)

struct GtIntset<%=bits%> {
  GtIntset parent_instance;
  uint<%=bits%>_t *elements;
};

GtIntset* gt_intset_<%=bits%>_new(GtUword maxelement, GtUword num_of_elems)
{
  GtIntset *intset;
  GtIntset<%=bits%> *intset_<%=bits%>;
  GtIntsetMembers *members;
  GtUword idx;

  gt_assert(num_of_elems != 0);

  intset = gt_intset_create(gt_intset_<%=bits%>_class());
  intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  members = intset->members;

  members->currentsectionnum = 0;
  members->maxelement = maxelement;
  members->nextfree = 0;
  members->num_of_elems = num_of_elems;
  members->previouselem = ULONG_MAX;
  members->refcount = 0;

  members->logsectionsize = GT_BITS_FOR_TYPE(uint<%=bits%>_t);
  members->numofsections = GT_ELEM2SECTION_M(maxelement) + 1;

  intset_<%=bits%>->elements =
    gt_malloc(sizeof (*intset_<%=bits%>->elements) * num_of_elems);

  members->sectionstart = gt_malloc(sizeof (*members->sectionstart) *
                                    (members->numofsections + 1));

  members->sectionstart[0] = 0;
  for (idx = (GtUword) 1; idx <= members->numofsections; idx++) {
    members->sectionstart[idx] = num_of_elems;
  }
  return intset;
}

static bool gt_intset_<%=bits%>_elems_is_valid(GtIntset *intset)
{
  GtUword idx, sec_idx = 0;
  GtIntset<%=bits%> *intset_<%=bits%>;
  GtIntsetMembers *members;
  intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  members = intset->members;

  for (idx = (GtUword) 1; idx < members->num_of_elems; idx++) {
    while (idx > members->sectionstart[sec_idx])
      sec_idx++;

    if (idx != members->sectionstart[sec_idx] &&
        intset_<%=bits%>->elements[idx] <= intset_<%=bits%>->elements[idx - 1])
      return false;
  }
  return true;
}

static bool gt_intset_<%=bits%>_secstart_is_valid(GtIntset *intset)
{
  GtUword idx;
  GtIntsetMembers *members;
  members = intset->members;

  for (idx = (GtUword) 1; idx <= members->numofsections; idx++) {
    if (members->sectionstart[idx] < members->sectionstart[idx - 1])
      return false;
  }
  return true;
}

GtIntset *gt_intset_<%=bits%>_io_fp(GtIntset *intset, FILE *fp, GtError *err,
                                    GtIOFunc io_func)
{
  int had_err = 0;
  GtUword type = (GtUword) GT_INTSET_<%=bits%>_TYPE;
  GtIntset<%=bits%> *intset_<%=bits%>;
  GtIntsetMembers *members;

  gt_error_check(err);

  intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);

  had_err = gt_intset_<%=bits%>_io_one(type);
  if (!had_err && type != GT_INTSET_<%=bits%>_TYPE) {
    /* only applies to reading */
    had_err = 1;
    gt_error_set(err, "Trying to read GtIntset<%=bits%> from file,"
                 " type does not match!");
  }
  if (!had_err) {
    members = intset->members;
    had_err = gt_intset_<%=bits%>_io_one(members->currentsectionnum);
    if (!had_err)
      had_err = gt_intset_<%=bits%>_io_one(members->maxelement);
    if (!had_err)
      had_err = gt_intset_<%=bits%>_io_one(members->nextfree);
    if (!had_err)
      had_err = gt_intset_<%=bits%>_io_one(members->num_of_elems);
    if (!had_err)
      had_err = gt_intset_<%=bits%>_io_one(members->previouselem);
    if (!had_err) {
      members->logsectionsize = GT_BITS_FOR_TYPE(uint<%=bits%>_t);
      members->numofsections = GT_ELEM2SECTION_M(members->maxelement) + 1;
      members->sectionstart = gt_realloc(members->sectionstart,
                sizeof (*members->sectionstart) * (members->numofsections + 1));
    }
    had_err = io_func(members->sectionstart, sizeof (*members->sectionstart),
                      (size_t) (members->numofsections + 1), fp, err);
    if (!had_err && members->sectionstart[0] != 0) {
      had_err = 1;
      gt_error_set(err, "Unexpected value in sectionstart[0]: "
                   GT_WU " expected 0!", members->sectionstart[0]);
    }
  }
  if (!had_err) {
    intset_<%=bits%>->elements = gt_realloc(intset_<%=bits%>->elements,
                  sizeof (*intset_<%=bits%>->elements) * members->num_of_elems);
    had_err = io_func(intset_<%=bits%>->elements,
                      sizeof (*intset_<%=bits%>->elements),
                      (size_t) members->num_of_elems, fp, err);
  }
  if (had_err) {
    gt_intset_<%=bits%>_delete(intset);
    intset = NULL;
  }
  return intset;

}

GtIntset *gt_intset_<%=bits%>_io(GtIntset *intset, FILE *fp, GtError *err)
{
  GtIntset<%=bits%> *intset_<%=bits%>;
  if (intset == NULL) {
    intset = gt_intset_create(gt_intset_<%=bits%>_class());
    intset->members->sectionstart = NULL;
    intset->members->refcount = 0;
    intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
    intset_<%=bits%>->elements = NULL;
    intset = gt_intset_<%=bits%>_io_fp(intset, fp, err, gt_io_error_fread);
  }
  else {
    intset = gt_intset_<%=bits%>_io_fp(intset, fp, err, gt_io_error_fwrite);
  }
  return intset;
}

GtIntset *gt_intset_<%=bits%>_new_from_file(FILE *fp, GtError *err)
{
  gt_assert(fp != NULL);
  gt_error_check(err);
  return gt_intset_<%=bits%>_io(NULL, fp, err);
}

GtIntset *gt_intset_<%=bits%>_write(GtIntset *intset, FILE *fp, GtError *err)
{
  gt_assert(intset != NULL);
  gt_assert(fp != NULL);
  gt_error_check(err);
  return gt_intset_<%=bits%>_io(intset, fp, err);
}

void gt_intset_<%=bits%>_add(GtIntset *intset, GtUword elem)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *sectionstart = members->sectionstart;
  gt_assert(members->nextfree < members->num_of_elems);
  gt_assert(elem <= members->maxelement);
  gt_assert(members->previouselem == ULONG_MAX || members->previouselem < elem);
  while (elem >= GT_SECTIONMINELEM(members->currentsectionnum + 1)) {
    gt_assert(members->currentsectionnum < members->numofsections);
    sectionstart[members->currentsectionnum + 1] = members->nextfree;
    members->currentsectionnum++;
  }
  gt_assert(GT_SECTIONMINELEM(members->currentsectionnum) <= elem &&
            elem < GT_SECTIONMINELEM(members->currentsectionnum+1) &&
            GT_ELEM2SECTION_M(elem) == members->currentsectionnum);
  intset_<%=bits%>->elements[members->nextfree++] = (uint<%=bits%>_t) elem;
  members->previouselem = elem;
}

static GtUword gt_intset_<%=bits%>_sec_idx_largest_leq(GtUword *sectionstart,
                                               \
<% if bits != 8 %> <% end %>GtUword idx)
{
  GtUword result = 0;
  while (sectionstart[result] <= idx)
    result++;
  return result - 1;
}

static GtUword
gt_intset_<%=bits%>_binarysearch_sec_idx_largest_leq(GtUword *secstart_begin,
                                             \
<% if bits != 8 %> <% end %>GtUword *secstart_end,
                                             \
<% if bits != 8 %> <% end %>GtUword idx)
{
  GtUword *midptr = NULL, *found = NULL,
          *startorig = secstart_begin;
  if (*secstart_begin <= idx)
    found = secstart_begin;
  while (secstart_begin < secstart_end) {
    midptr = secstart_begin + (GtUword) GT_DIV2(secstart_end - secstart_begin);
    if (*midptr < idx) {
      found = midptr;
      if (*midptr == idx) {
        break;
      }
      secstart_begin = midptr + 1;
    }
    else {
      secstart_end = midptr - 1;
    }
  }
  gt_assert(found != NULL);
  while (found[1] <= idx)
    found++;
  return (GtUword) (found - startorig);
}

static GtUword gt_intset_<%=bits%>_get_test(GtIntset *intset, GtUword idx)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *sectionstart = members->sectionstart;
  gt_assert(idx < members->nextfree);

  return (gt_intset_<%=bits%>_sec_idx_largest_leq(sectionstart, idx) <<
         members->logsectionsize) + intset_<%=bits%>->elements[idx];
}

GtUword gt_intset_<%=bits%>_get(GtIntset *intset, GtUword idx)
{
  GtUword quotient;
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *sectionstart = members->sectionstart;
  gt_assert(idx < members->nextfree);

  quotient = gt_intset_<%=bits%>_binarysearch_sec_idx_largest_leq(
                                      sectionstart,
                                      sectionstart + members->numofsections - 1,
                                      idx);
  return (quotient << members->logsectionsize) +
         intset_<%=bits%>->elements[idx];
}

GtUword gt_intset_<%=bits%>_size(GtIntset *intset)
{
  GT_UNUSED GtIntset<%=bits%> *intset_<%=bits%> = \
gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;
  return members->nextfree;
}

static bool gt_intset_<%=bits%>_binarysearch_is_member(\
const uint<%=bits%>_t *leftptr,
<% if bits != 8 %> <%end%>                                               \
const uint<%=bits%>_t *rightptr,
<% if bits != 8 %> <%end%>                                               \
uint<%=bits%>_t elem)
{
  const uint<%=bits%>_t *midptr;
    while (leftptr <= rightptr) {
      midptr = leftptr + (GtUword) GT_DIV2(rightptr - leftptr);
      if (elem < *midptr) {
        rightptr = midptr - 1;
      }
      else {
        if (elem > *midptr)
          leftptr = midptr + 1;
        else
          return true;
      }
    }
  return false;
}

bool gt_intset_<%=bits%>_is_member(GtIntset *intset, GtUword elem)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;
  GtUword *sectionstart = members->sectionstart;
  if (elem <= members->maxelement)
  {
    const GtUword sectionnum = GT_ELEM2SECTION_M(elem);

    if (sectionstart[sectionnum] < sectionstart[sectionnum+1]) {
      return gt_intset_<%=bits%>_binarysearch_is_member(
<% if bits == 8 %> <% end %>                           \
intset_<%=bits%>->elements + sectionstart[sectionnum],
<% if bits == 8 %> <% end %>                           \
intset_<%=bits%>->elements + sectionstart[sectionnum+1] - 1,
<% if bits == 8 %> <% end %>                           (uint64_t) elem);
    }
  }
  return false;
}

static GtUword gt_intset_<%=bits%>_idx_sm_geq(\
const uint<%=bits%>_t *leftptr,
<% if bits != 8 %> <%end%>                                      \
const uint<%=bits%>_t *rightptr,
<% if bits != 8 %> <%end%>                                      \
uint<%=bits%>_t value)
{
  const uint<%=bits%>_t *leftorig = leftptr;
  if (value < *leftptr)
    return 0;
  if (value > *rightptr)
    return 1UL + (GtUword) (rightptr - leftptr);
  gt_assert(value <= *rightptr);
  while (*leftptr < value)
    leftptr++;
  return (GtUword) (leftptr - leftorig);
}

static GtUword gt_intset_<%=bits%>_binarysearch_idx_sm_geq(\
const uint<%=bits%>_t *leftptr,
<% if bits != 8 %> <%end%>                                                   \
const uint<%=bits%>_t *rightptr,
<% if bits != 8 %> <%end%>                                                   \
uint<%=bits%>_t value)
{
  const uint<%=bits%>_t *midptr = NULL,
        *leftorig = leftptr;

  gt_assert(leftptr <= rightptr);
  if (value <= *leftptr)
    return 0;
  if (value > *rightptr)
    return 1UL + (GtUword) (rightptr - leftptr);
  while (leftptr < rightptr) {
    midptr = leftptr + (GtUword) GT_DIV2(rightptr - leftptr);
    if (value <= *midptr)
      rightptr = midptr;
    else {
      leftptr = midptr + 1;
    }
  }
  return (GtUword) (leftptr - leftorig);
}

static GtUword gt_intset_<%=bits%>_get_idx_smallest_geq_test(GtIntset *intset,
<% if bits == 8 %> <% end %>                                                   \
  GtUword value)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;

  GtUword sectionnum = GT_ELEM2SECTION_M(value);

  if (value > members->previouselem)
    return members->num_of_elems;

  gt_assert(value <= members->maxelement);
  if (members->sectionstart[sectionnum] < members->sectionstart[sectionnum+1]) {
    return members->sectionstart[sectionnum] +
           gt_intset_<%=bits%>_idx_sm_geq(
<% if bits == 8 %> <% end %>                  intset_<%=bits%>->elements + \
members->sectionstart[sectionnum],
<% if bits == 8 %> <% end %>                  intset_<%=bits%>->elements + \
members->sectionstart[sectionnum+1] - 1,
<% if bits == 8 %> <% end %>                  (uint<%=bits%>_t) value);
  }
  return members->sectionstart[sectionnum];
}

GtUword gt_intset_<%=bits%>_get_idx_smallest_geq(GtIntset *intset, \
GtUword value)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  GtIntsetMembers *members = intset->members;

  GtUword sectionnum = GT_ELEM2SECTION_M(value);

  if (value > members->previouselem)
    return members->num_of_elems;

  gt_assert(value <= members->maxelement);

  if (members->sectionstart[sectionnum] < members->sectionstart[sectionnum+1]) {
    return members->sectionstart[sectionnum] +
           gt_intset_<%=bits%>_binarysearch_idx_sm_geq(
<% if bits == 8 %> <% end %>                  intset_<%=bits%>->elements + \
members->sectionstart[sectionnum],
<% if bits == 8 %> <% end %>                  intset_<%=bits%>->elements + \
members->sectionstart[sectionnum+1] - 1,
<% if bits == 8 %> <% end %>                  (uint<%=bits%>_t) value);
  }
  return members->sectionstart[sectionnum];
}

void gt_intset_<%=bits%>_delete(GtIntset *intset)
{
  GtIntset<%=bits%> *intset_<%=bits%> = gt_intset_<%=bits%>_cast(intset);
  if (intset_<%=bits%> != NULL) {
    gt_free(intset_<%=bits%>->elements);
  }
}

size_t gt_intset_<%=bits%>_size_of_rep(GtUword maxelement, GtUword num_of_elems)
{
  size_t logsectionsize = GT_BITS_FOR_TYPE(uint<%=bits%>_t);
  gt_assert(GT_BITS_FOR_TYPE(GtUword) > logsectionsize);
  return sizeof (uint<%=bits%>_t) * num_of_elems +
    sizeof (GtUword) * (GT_ELEM2SECTION(maxelement, logsectionsize) + 1);
}

size_t gt_intset_<%=bits%>_size_of_struct(void)
{
  return sizeof (GtIntset<%=bits%>) +
         sizeof (struct GtIntsetClass) +
         sizeof (struct GtIntsetMembers);
}

bool gt_intset_<%=bits%>_file_is_type(GtUword type)
{
  return type == GT_INTSET_<%=bits%>_TYPE;
}

/* map static local methods to interface */
const GtIntsetClass* gt_intset_<%=bits%>_class(void)
{
  static const GtIntsetClass *this_c = NULL;
  if (this_c == NULL) {
    this_c = gt_intset_class_new(sizeof (GtIntset<%=bits%>),
                                 gt_intset_<%=bits%>_add,
                                 gt_intset_<%=bits%>_file_is_type,
                                 gt_intset_<%=bits%>_get,
                                 gt_intset_<%=bits%>_io,
                                 gt_intset_<%=bits%>_get_idx_smallest_geq,
                                 gt_intset_<%=bits%>_is_member,
                                 gt_intset_<%=bits%>_size_of_rep,
                                 gt_intset_<%=bits%>_size,
                                 gt_intset_<%=bits%>_size_of_struct,
                                 gt_intset_<%=bits%>_write,
                                 gt_intset_<%=bits%>_delete);
  }
  return this_c;
}

#define GT_INTSET_TEST_<%=bits%>_BINSEARCH(IDX) \\
gt_ensure(gt_intset_<%=bits%>_get_test(is, IDX) == arr[IDX]); \\
gt_ensure(gt_intset_<%=bits%>_get(is, IDX) == arr[IDX]); \\
gt_ensure( \\
  gt_intset_<%=bits%>_get_idx_smallest_geq_test(is, arr[IDX] + 1) == \\
  IDX + 1); \\
gt_ensure( \\
  gt_intset_<%=bits%>_get_idx_smallest_geq(is, arr[IDX] + 1) == \\
  IDX + 1)

int gt_intset_<%=bits%>_unit_test(GtError *err)
{
  int had_err = 0;
  GtIntset *is;
  GtUword num_of_elems = gt_rand_max(((GtUword) 1) << 10) + 1,
          *arr = gt_malloc(sizeof (*arr) * num_of_elems),
          stepsize = GT_DIV2(num_of_elems <<4 / num_of_elems),
          idx;
  size_t is_size;

  gt_error_check(err);

  arr[0] = gt_rand_max(stepsize) + 1;
  for (idx = (GtUword) 1; idx < num_of_elems; ++idx) {
    arr[idx] = arr[idx - 1] + gt_rand_max(stepsize) + 1;
  }

  is_size = \
gt_intset_<%=bits%>_size_of_rep(arr[num_of_elems - 1], num_of_elems);

  if (!had_err) {
    if (is_size < (size_t) UINT_MAX) {
      is = gt_intset_<%=bits%>_new(arr[num_of_elems - 1], num_of_elems);
      for (idx = 0; idx < num_of_elems; idx++) {
        gt_intset_<%=bits%>_add(is, arr[idx]);
        gt_ensure(idx + 1 == gt_intset_<%=bits%>_size(is));
        if (idx < num_of_elems - 1)
          gt_ensure(gt_intset_<%=bits%>_get_idx_smallest_geq(is,
                                                     \
<% if bits != 8 %> <% end %>arr[idx] + 1) ==
                    num_of_elems);
      }

      gt_ensure(gt_intset_<%=bits%>_elems_is_valid(is));
      gt_ensure(gt_intset_<%=bits%>_secstart_is_valid(is));

      GT_INTSET_TEST_<%=bits%>_BINSEARCH(0);
      for (idx = 1; !had_err && idx < num_of_elems; idx++) {
        GtUword to_find = (arr[idx - 1] == (arr[idx] - 1)) ? idx - 1 : idx;
        gt_ensure(
          gt_intset_<%=bits%>_get_idx_smallest_geq_test(is, arr[idx] - 1) ==
          to_find);
        gt_ensure(
          gt_intset_<%=bits%>_get_idx_smallest_geq(is, arr[idx] - 1) ==
          to_find);
        GT_INTSET_TEST_<%=bits%>_BINSEARCH(idx);
      }
      if (!had_err)
        had_err = gt_intset_unit_test_notinset(is, 0, arr[0] - 1, err);
      if (!had_err)
        had_err = gt_intset_unit_test_check_seqnum(is, 0, arr[0] - 1, 0, err);
      for (idx = (GtUword) 1; !had_err && idx < num_of_elems; idx++) {
        had_err = gt_intset_unit_test_notinset(is, arr[idx - 1] + 1,
                                               arr[idx] - 1, err);
        if (!had_err)
          had_err = gt_intset_unit_test_check_seqnum(is, arr[idx - 1] + 1,
                                                     arr[idx] - 1, idx, err);
      }
      gt_intset_delete(is);
    }
  }
  gt_free(arr);
  return had_err;
}
IMPL
#end=c

#begin=c
HEADER = <<-HEADER
#ifndef INTSET_<%=bits%>_H
#define INTSET_<%=bits%>_H

#include "core/types_api.h"
#include "extended/intset_rep.h"

/* The <GtIntset<%=bits%>> class implements the <GtIntset> interface.
   This class only works if <GtUword> is larger than <%=bits%> bits! */
typedef struct GtIntset<%=bits%> GtIntset<%=bits%>;

/* Map static local methods to interface */
const     GtIntsetClass* gt_intset_<%=bits%>_class(void);

/* Return a new <GtIntset> object, the implementation beeing of type
   <GtIntset<%=bits%>>.
   Fails if <%=bits%> >= bits for (GtUword). */
GtIntset* gt_intset_<%=bits%>_new(GtUword maxelement, GtUword num_of_elems);

/* Returns true, if the <type> read from a file storing a <GtIntset> indicates
   the type of this implementation. */
bool      gt_intset_<%=bits%>_file_is_type(GtUword type);

/* Return a new <GtIntset> object, with data read from <fp> */
GtIntset* gt_intset_<%=bits%>_new_from_file(FILE *fp, GtError *err);

/* Add <elem> to <intset>. <elem> has to be larger than the previous <elem>
   added. */
void      gt_intset_<%=bits%>_add(GtIntset *intset, GtUword elem);

/* Returns the element at index <idx> in the sorted set <intset>. */
GtUword   gt_intset_<%=bits%>_get(GtIntset *intset, GtUword idx);

/* Returns actual number of stored elements */
GtUword   gt_intset_<%=bits%>_size(GtIntset *intset);

/* Returns <true> if <elem> is a member of the set <intset>. */
bool      gt_intset_<%=bits%>_is_member(GtIntset *intset, GtUword elem);

/* Returns the number of the element in <intset> that is the smallest element
   larger than or equal <value> or <num_of_elems> if there is no such <element>.
   */
GtUword   gt_intset_<%=bits%>_get_idx_smallest_geq(GtIntset *intset, \
GtUword value);

/* Write <intset> to file <fp>. Returns <NULL> on error (<intset> will be
   freed). */
GtIntset* gt_intset_<%=bits%>_write(GtIntset *intset, FILE *fp, GtError *err);

/* IO-function to be used if <intset> is part of a larger structure. If <intset>
   is NULL, will attempt to allocate memory and fill a new <GtIntset<%=bits%>>
   object by reading from <fp>. If <intset> is not NULL, will attempt to write
   its content to <fp>.
   Returns <NULL> on error (<intset> will be freed) and sets <err>. */
GtIntset* gt_intset_<%=bits%>_io(GtIntset *intset, FILE *fp, GtError *err);

/* Deletes <intset> and frees all associated space. */
void      gt_intset_<%=bits%>_delete(GtIntset *intset);

/* Returns the size of the representation of an intset with given number of
   elements <num_of_elems> and maximum value <maxelement>, in bytes. This does
   not include the size of the structure.
   Fails if <%=bits%> >= bits for (GtUword). */
size_t    gt_intset_<%=bits%>_size_of_rep(GtUword maxelement, \
GtUword num_of_elems);

/* Returns the size in bytes of the <GtIntset<%=bits%>>-structure. */
size_t    gt_intset_<%=bits%>_size_of_struct(void);

int gt_intset_<%=bits%>_unit_test(GtError *err);
#endif
HEADER
#end=c

name, email = CodeGen.find_user_data()
year = Time.now.year
all_bits = [8, 16, 32]

all_bits.each do |bits|
  fn = "src/extended/intset_#{bits}"
  code = File.open(fn + '.c', 'w')
  code.puts ERB.new(CRIGHT).result(binding)
  code.puts ERB.new(IMPL).result(binding)
  code.close
  header = File.open(fn + '.h', 'w')
  header.puts ERB.new(CRIGHT).result(binding)
  header.puts ERB.new(HEADER).result(binding)
  header.close
end
# if using vim, and if the addon SyntaxRange is available use this command to
# change syntax highlighting to c:
# :call SyntaxRange#Include('^#begin=c', '^#end=c', 'c', 'NonText')
