/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/ensure.h"
#include "libgtcore/log.h"
#include "libgtcore/ma.h"
#include "libgtview/block.h"
#include "libgtview/element.h"

struct Block {
  Dlist *elements;
  Range range;
  Str *caption;
  bool show_caption;
  Strand strand;
  GenomeFeatureType type;
};

/* Compare function used to insert Elements into dlist */
static int elemcmp(const void *a, const void *b)
{
  Element *elem_a = (Element*) a;
  Element *elem_b = (Element*) b;

  Range ra = element_get_range(elem_a);
  Range rb = element_get_range(elem_b);

  return range_compare(ra, rb);
}

Block* block_new(Env *env)
{
  Block *block;
  Range r;
  env_error_check(env);
  block = ma_malloc(sizeof (Block));
  block->elements = dlist_new(elemcmp, env);
  r.start = 0;
  r.end = 0;
  block->range = r;
  block->caption = NULL;
  block->show_caption = true;
  block->strand = STRAND_UNKNOWN;
  return block;
}

Block* block_new_from_node(GenomeNode *node, Env *env)
{
  Block *block;
  assert(node);
  block = block_new(env);
  block_set_range(block, genome_node_get_range(node));
  block_set_strand(block, genome_feature_get_strand((GenomeFeature*) node));
  return block;
}

void block_insert_element(Block *block, GenomeNode *gn, Config *cfg, Env *env)
{
  Dlistelem *elem;
  Range elem_r, gn_r;
  DominateStatus dominates;
  int count = 0;
  Element *element, *e;
  GenomeFeatureType gn_type, e_type;

  assert(block && gn && cfg);

  gn_r = genome_node_get_range(gn);
  gn_type = genome_feature_get_type((GenomeFeature*) gn);

  log_log("inserting %s (%lu-%lu) into block",
          genome_feature_type_get_cstr(gn_type),
          gn_r.start, gn_r.end);

  for (elem = dlist_first(block->elements); elem;
       elem = dlistelem_next(elem)) {
    element = (Element*) dlistelem_get_data(elem);
    elem_r = element_get_range(element);

    if (range_overlap(elem_r, gn_r)) {
      count += 1;
      e_type = element_get_type(element);

      dominates = config_dominates(cfg, e_type, gn_type, env);

      if (gn_r.start == elem_r.start && gn_r.end < elem_r.end) {
        /* Case:    -------------------
                    ---------- */
        switch (dominates) {
          case DOMINATES_FIRST:
            break;
          case DOMINATES_NOT_SPECIFIED:
          case DOMINATES_UNKNOWN_TYPE:
          case DOMINATES_EQUAL:
          case DOMINATES_SECOND:
            elem_r.start = gn_r.end+1;
            if (elem_r.start == elem_r.end+1) {
              dlist_remove(block->elements, elem, env);
              element_delete(element, env);
            }
            else
              element_set_range(element, elem_r);
            e = element_new_empty(env);
            element_set_range(e, gn_r);
            element_set_type(e, gn_type);
            dlist_add(block->elements, e, env);
            elem = dlist_find(block->elements, e);
            break;
        }
      }
      else if (gn_r.start >= elem_r.start && gn_r.end == elem_r.end) {
        /* Case:  --------------
                     -----------  */
        switch (dominates) {
          case DOMINATES_FIRST:
            gn_r.start = elem_r.end;
            break;
          case DOMINATES_NOT_SPECIFIED:
          case DOMINATES_UNKNOWN_TYPE:
          case DOMINATES_EQUAL:
          case DOMINATES_SECOND:
            elem_r.end = gn_r.start-1;
            if (elem_r.start == elem_r.end+1) {
              dlist_remove(block->elements, elem, env);
              element_delete(element, env);
            }
            else
              element_set_range(element, elem_r);
            e = element_new(gn, env);
            element_set_range(e, gn_r);
            element_set_type(e, gn_type);
            dlist_add(block->elements, e, env);
            elem = dlist_find(block->elements, e);
            break;
        }
      }
      else if (elem_r.start <= gn_r.start && elem_r.end < gn_r.end) {
        /* Case: ----------
                 -------------- */
        bool removed = false;
        unsigned long tmp;
        Range gnnew_r;

        switch (dominates) {
          case DOMINATES_FIRST:
            gn_r.start = elem_r.end+1;
            break;
          case DOMINATES_EQUAL:
          case DOMINATES_UNKNOWN_TYPE:
          case DOMINATES_NOT_SPECIFIED:
          case DOMINATES_SECOND:
            tmp = elem_r.end;
            elem_r.end = gn_r.start-1;
            if (elem_r.start == elem_r.end+1) {
              dlist_remove(block->elements, elem, env);
              element_delete(element, env);
              removed = true;
              elem_r.end = tmp;
            }
            else
              element_set_range(element, elem_r);
            gnnew_r = gn_r;
            gnnew_r.end = elem_r.end;
            e = element_new_empty(env);
            element_set_range(e, gnnew_r);
            element_set_type(e, e_type);
            dlist_add(block->elements, e, env);
            gn_r.start = elem_r.end+1;
            if (removed)
              elem = dlist_find(block->elements, e);
            break;
        }
      }
      else if (elem_r.start < gn_r.start && gn_r.end < elem_r.end) {
        /* Case: -------------
                    ------      */
        Range elemnew_r;
        Element *elemnew;

        switch (dominates) {
          case DOMINATES_FIRST:
            break;
          case DOMINATES_EQUAL:
          case DOMINATES_SECOND:
          case DOMINATES_NOT_SPECIFIED:
          case DOMINATES_UNKNOWN_TYPE:
            elemnew_r = elem_r;
            elem_r.end = gn_r.start-1;
            element_set_range(element, elem_r);
            elemnew_r.start = gn_r.end+1;
            elemnew = element_new_empty(env);
            element_set_range(elemnew, elemnew_r);
            element_set_type(elemnew, e_type);
            e = element_new(gn, env);
            dlist_add(block->elements, elemnew, env);
            dlist_add(block->elements, e, env);
            break;
        }
      }
    }
  }
  if (count == 0) {
    e = element_new(gn, env);
    dlist_add(block->elements, e, env);
  }
}

Range block_get_range(const Block *block)
{
   assert(block);
   return block->range;
}

void block_set_range(Block *block, Range r)
{
  assert(block && r.start <= r.end);
  block->range = r;
}

bool block_has_only_one_fullsize_element(Block *block)
{
  bool ret = false;
  assert(block);
  if (dlist_size(block->elements) == 1UL) {
    Range elem_range, block_range;
    assert(dlist_first(block->elements) == dlist_last(block->elements));
    elem_range = element_get_range(
                   dlistelem_get_data(dlist_first(block->elements)));
    block_range = block_get_range(block);
    ret = (range_compare(block_range, elem_range) == 0);
  }
  return ret;
}

void block_set_caption_visibility(Block *block, bool val)
{
  assert(block);
  block->show_caption = val;
}

bool block_caption_is_visible(const Block *block)
{
  assert(block);
  return block->show_caption;
}

void block_set_caption(Block *block, Str *caption)
{
  assert(block && caption);
  block->caption = caption;
}

Str* block_get_caption(const Block *block)
{
  assert(block);
  return block->caption;
}

void block_set_strand(Block *block, Strand strand)
{
  assert(block);
  block->strand = strand;
}

Strand block_get_strand(const Block *block)
{
  assert(block);
  return block->strand;
}

void block_set_type(Block *block, GenomeFeatureType type)
{
  assert(block);
  block->type = type;
}

GenomeFeatureType block_get_type(const Block *block)
{
  assert(block);
  return block->type;
}

Dlist* block_get_elements(const Block *block)
{
  assert(block);
  return block->elements;
}

int block_unit_test(Env* env)
{
  Range r1, r2, r_temp, b_range;
  Dlist* elements;
  int had_err = 0;
  Strand s;
  GenomeNode *gn1, *gn2;
  Element *e1, *e2, *elem;
  Block * b;
  Str *caption1 = str_new_cstr("foo", env);
  Str *caption2 = str_new_cstr("bar", env);
  Config *cfg;

  cfg = config_new(false, env);

  r1.start = 10UL;
  r1.end = 50UL;

  r2.start = 51UL;
  r2.end = 80UL;

  gn1 = genome_feature_new(gft_exon, r1, STRAND_FORWARD, NULL, 0, env);
  gn2 = genome_feature_new(gft_intron, r2, STRAND_FORWARD, NULL, 0, env);

  e1 = element_new(gn1, env);
  e2 = element_new(gn2, env);

  b = block_new(env);

  /* test block_insert_elements */
  ensure(had_err, (0UL == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn1, cfg, env);
  ensure(had_err, (1UL == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn2, cfg, env);
  ensure(had_err, (2UL == dlist_size(block_get_elements(b))));

  /* test block_get_elements */
  elements = block_get_elements(b);
  elem = (Element*) dlistelem_get_data(dlist_first(elements));
  ensure(had_err, elements_are_equal(e1, elem));
  ensure(had_err, !elements_are_equal(e2, (Element*) dlist_first(elements)));
  elem = (Element*) dlistelem_get_data(dlist_last(elements));
  ensure(had_err, !elements_are_equal(e1, elem));
  ensure(had_err, elements_are_equal(e2, elem));

  /* test block_set_range & block_get_range */
  r_temp = range_join(r1, r2);
  block_set_range(b, r_temp);
  b_range = block_get_range(b);
  ensure(had_err, (0 == range_compare(b_range, r_temp)));
  ensure(had_err, (1 == range_compare(r2, r_temp)));

  /* tests block_set_caption & block_get_caption */
  block_set_caption(b, caption1);
  ensure(had_err, (0 == str_cmp(block_get_caption(b), caption1)));
  ensure(had_err, (0 != str_cmp(block_get_caption(b), caption2)));

  /* tests block_set_strand & block_get_range */
  s = block_get_strand(b);
  ensure(had_err, (STRAND_UNKNOWN == s));
  block_set_strand(b, STRAND_FORWARD);
  s = block_get_strand(b);
  ensure(had_err, (STRAND_FORWARD == s));

  config_delete(cfg, env);
  str_delete(caption2);
  element_delete(e1, env);
  element_delete(e2, env);
  block_delete(b, env);
  genome_node_delete(gn1, env);
  genome_node_delete(gn2, env);

  return had_err;
}

void block_delete(Block *block,
                  Env *env)
{
  Dlistelem *delem;
  if (!block) return;
  for (delem = dlist_first(block->elements); delem;
       delem = dlistelem_next(delem)) {
    Element* elem = (Element*) dlistelem_get_data(delem);
    element_delete(elem, env);
  }
  str_delete(block->caption);
  dlist_delete(block->elements, env);
  ma_free(block);
}
