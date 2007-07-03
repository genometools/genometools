/*
   Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
   Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
   See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtview/block.h>
#include <libgtview/element.h>

struct Block
{
  Dlist *elements;
  Range range;
  const char* caption;
  const char* parent_caption;
  bool show_caption;
  Strand strand;
  GenomeFeatureType type;
};

void block_set_caption_visibility(Block *b, bool val)
{
  assert(b);
  b->show_caption = val;
}

bool block_caption_is_visible(Block *b)
{
  assert(b);
  return b->show_caption;
}

/*
Compare function to insert Elements into dlist
*/
int elemcmp(const void *a, const void *b)
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
  env_error_check(env);
  block = env_ma_malloc(env, sizeof (Block));
  block->elements = dlist_new(elemcmp, env);
  Range r;
  r.start = 0;
  r.end = 0;
  block->range = r;
  block->caption = NULL;
  block->parent_caption = NULL;
  block->show_caption = true;
  block->strand = STRAND_UNKNOWN;

  assert(block);
  return block;
}

void block_insert_element(Block *block,
                          GenomeNode *gn,
			  Config *cfg,
			  Env *env)
{
  assert(block && gn && cfg);

  Dlistelem *elem;
  Range elem_r, gn_r;
  int dominates;
  int count = 0;
  Element *element, *e;
  GenomeFeatureType gn_type, e_type;

  gn_r = genome_node_get_range(gn);
  gn_type = genome_feature_get_type((GenomeFeature*) gn);

  for (elem = dlist_first(block->elements); elem != NULL;
      elem = dlistelem_next(elem))
  {
    element = (Element*) dlistelem_get_data(elem);
    elem_r = element_get_range(element);

    if (range_overlap(elem_r, gn_r))
    {
      count += 1;
      e_type = element_get_type(element);

      dominates = config_dominates(cfg, e_type, gn_type, env);
      if (dominates == DOMINATES_EQUAL
         || dominates == DOMINATES_NOT_SPECIFIED
	 || dominates == DOMINATES_UNKNOWN_TYPE)
      {
        dominates = DOMINATES_SECOND;
      }

      /* Fall:    -------------------
                  ---------- */
      if (gn_r.start == elem_r.start && gn_r.end < elem_r.end)
      {
        switch (dominates)
        {
          case DOMINATES_FIRST:
	    break;
	  case DOMINATES_SECOND:
            elem_r.start = gn_r.end+1;
	    element_set_range(element, elem_r);
	    e = element_new(gn, cfg, env);
	    dlist_add(block->elements, e, env);
	    break;
        }
      }

      /* Fall:  --------------
                   -----------  */
      else if (gn_r.start >= elem_r.start && gn_r.end == elem_r.end)
      {
        switch (dominates)
	{
          case DOMINATES_FIRST:
	    gn_r.start = elem_r.end;
	    break;
	  case DOMINATES_SECOND:
	    elem_r.end = gn_r.start-1;
	    if (elem_r.start == elem_r.end+1)
	    {
              dlist_remove(block->elements, elem, env);
	      element_delete(element, env);
	    }
	    else
	    {
              element_set_range(element, elem_r);
	    }
	    e = element_new(gn, cfg, env);
            element_set_range(e, gn_r);
	    dlist_add(block->elements, e, env);
	    elem = dlist_find(block->elements, e);
	    break;
	}
      }

      /* Fall: ----------
               -------------- */
      else if (elem_r.start <= gn_r.start && elem_r.end < gn_r.end)
      {
        bool removed = false;
        switch (dominates)
	{
          case DOMINATES_FIRST:
            gn_r.start = elem_r.end+1;
	    break;
	  case DOMINATES_SECOND:
            elem_r.end = gn_r.start-1;
	    if (elem_r.start == elem_r.end+1)
	    {
              dlist_remove(block->elements, elem, env);
	      element_delete(element, env);
	      removed = true;
	    }
	    else
            {
	      element_set_range(element, elem_r);
	    }
	    Range gnnew_r = gn_r;
	    gnnew_r.end = elem_r.end;
	    e = element_new_empty(cfg, env);
	    element_set_range(e, gnnew_r);
	    element_set_type(e, gn_type);
	    dlist_add(block->elements, e, env);
	    gn_r.start = elem_r.end+1;
	    if (removed)
	    {
              elem = dlist_find(block->elements, e);
	    }
	    break;
	}
      }

      /* Fall: -------------
                  ------      */
      else if (elem_r.start < gn_r.start && gn_r.end < elem_r.end)
      {
        Range elemnew_r;

        switch (dominates)
	{
          case DOMINATES_FIRST:
	    break;
	  case DOMINATES_SECOND:
            elemnew_r = elem_r;
	    elem_r.end = gn_r.start-1;
	    element_set_range(element, elem_r);
	    elemnew_r.start = gn_r.end+1;
            Element *elemnew = element_new_empty(cfg, env);
	    element_set_range(elemnew, elemnew_r);
	    element_set_type(elemnew, e_type);
	    e = element_new(gn, cfg, env);
	    dlist_add(block->elements, elemnew, env);
	    dlist_add(block->elements, e, env);
	    break;
	}
      }

    }

  }
  if (count == 0)
  {
    e = element_new(gn, cfg, env);
    dlist_add(block->elements, e, env);
  }
}

Range block_get_range(Block *block)
{
   assert(block);

   return block->range;
}

void block_set_range(Block *block, Range r)
{
  assert(block && r.start && r.end);

  block->range = r;
}

void block_set_caption(Block *block,
                       const char* caption)
{
  assert(block);
  block->caption = caption;
}

const char* block_get_caption(Block *block)
{
  assert(block);

  return block->caption;
}

void block_set_parent_caption(Block *block,
                       const char* caption)
{
  assert(block);
  block->parent_caption = caption;
}

const char* block_get_parent_caption(Block *block)
{
  assert(block);

  return block->parent_caption;
}

void block_set_strand(Block *block,
                      Strand strand)
{
  assert(block);

  block->strand = strand;
}

Strand block_get_strand(Block *block)
{
  assert(block);

  return block->strand;
}

void block_set_type(Block *block,
                    GenomeFeatureType type)
{
  assert(block);

  block->type = type;
}

GenomeFeatureType block_get_type(Block *block)
{
  assert(block);

  return block->type;
}

Dlist* block_get_elements(Block *block)
{
  return block->elements;
}

int block_unit_test(Env* env)
{
  Range r1, r2, r_temp, b_range;
  Dlist* elements;
  int had_err = 0;
  Strand s;

  Config *cfg;

  cfg = config_new(env, false);

  r1.start = 10;
  r1.end = 50;

  r2.start = 51;
  r2.end = 80;

  GenomeNode* gn1 = genome_feature_new(gft_exon, r1,
                                       STRAND_FORWARD, NULL, 0, env);
  GenomeNode* gn2 = genome_feature_new(gft_intron, r2,
                                       STRAND_FORWARD, NULL, 0, env);

  Element* e1 = element_new(gn1, cfg, env);
  Element* e2 = element_new(gn2, cfg, env);

  Block* b = block_new(env);

  const char* caption1 = "foo";
  const char* caption2 = "bar";

  /* test block_insert_elements */
  ensure(had_err, (0 == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn1, cfg, env);
  ensure(had_err, (1 == dlist_size(block_get_elements(b))));
  block_insert_element(b, gn2, cfg, env);
  ensure(had_err, (2 == dlist_size(block_get_elements(b))));

  /* test block_get_elements */
  elements = block_get_elements(b);
  Element *elem = (Element*) dlistelem_get_data(dlist_first(elements));
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

  /* tests block_set_caption
     & block_get_caption */
  block_set_caption(b, caption1);
  ensure(had_err, (0 == strcmp(block_get_caption(b), caption1)));
  ensure(had_err, (0 != strcmp(block_get_caption(b), caption2)));

  /* tests block_set_strand
     & block_get_range */
  s = block_get_strand(b);
  ensure(had_err, (STRAND_UNKNOWN == s));
  block_set_strand(b, STRAND_FORWARD);
  s = block_get_strand(b);
  ensure(had_err, (STRAND_FORWARD == s));

  config_delete(cfg, env);
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

  for (delem = dlist_first(block->elements); delem != NULL;
      delem = dlistelem_next(delem))
  {
    Element* elem = (Element*) dlistelem_get_data(delem);
    element_delete(elem, env);
  }
  dlist_delete(block->elements, env);
  env_ma_free(block, env);
}

