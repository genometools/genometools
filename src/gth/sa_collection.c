/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "extended/redblack.h"
#include "gth/sa_cmp.h"
#include "gth/sa_collection.h"

/*
  The following structure contains binary search trees of spliced alignments.
  They contain the same elements, but ordered according to different criteria.
  The tree rooted at <rootlist> contains the spliced alignments primarily
  ordered according to their left genomic position (as the linked list in GS2).
  Thereby, the genomic positions always refer to the forward strand (in
  contrary to the rest of the program, where genomic positions always refer to
  the actual strand). The elements of the tree rooted at <rootEST> are
  primarily ordered according to their reference ids. This tree is used to
  ensure that no alignment which ``spans'' another alignment (according to the
  function span_check()) is inserted into the trees.
*/
struct GthSACollection {
  GtRBTnode *rootlist,
            *rootEST;
};

GthSACollection* gth_sa_collection_new(void)
{
  return gt_calloc(1, sizeof (GthSACollection));
}

#define EVALRVAL(R)\
        if ((R) < 0)\
        {\
          return -1;\
        }\
        if ((R) > 0)\
        {\
          return  1;\
        }

typedef enum
{
  NOSPAN = 0,
  DISCARD,
  REPLACE,
  NUMOFSPANACTIONS
} Spancheck;

/*
  The following function checks if one of the two GthSA
  structures spans the other. This is checked on the basis of the left and
  right genomic border. If <saB> is included in <saA>, DISCARD is returned.
  The other way around REPLACE is returned. If none of the both spans the other
  NOSPAN is returned;
*/
static Spancheck span_check(GthSA *saA, GthSA *saB)
{
  GtRange rangeA, rangeB;

  rangeA = gth_sa_range_forward(saA);
  rangeB = gth_sa_range_forward(saB);

  if ((rangeA.start <= rangeB.start) && (rangeA.end >= rangeB.end)) {
    /* discard saB */
    return DISCARD;
  }
  if ((rangeA.start >= rangeB.start) && (rangeA.end <= rangeB.end)) {
    /* replace saA by saB */
    return REPLACE;
  }
  return NOSPAN;
}

static int compare_strands(bool strandsignA, bool strandsignB)
{
  if (strandsignA && !strandsignB)
    return -1;
  if (!strandsignA && strandsignB)
    return 1;
  return 0;
}

/*
  This function is used to check if the tree rooted at <rootEST>
  contains an alignment with the same ids and strand orientations.
  But the tree is ordered according to the function
  compare_duplicate_and_genomic_pos() given below.
*/

static int compare_duplicate(const GtKeytype dataA, const GtKeytype dataB,
                             GT_UNUSED void *cmpinfo)
{
  int rval;
  GthSA *saA = (GthSA*) dataA;
  GthSA *saB = (GthSA*) dataB;

  gt_assert(!cmpinfo);

  rval = gt_str_cmp(gth_sa_gen_id_str(saA), gth_sa_gen_id_str(saB));
  EVALRVAL(rval);

  rval = compare_strands(gth_sa_gen_strand_forward(saA),
                         gth_sa_gen_strand_forward(saB));
  EVALRVAL(rval);

  rval = gt_str_cmp(gth_sa_ref_id_str(saA), gth_sa_ref_id_str(saB));
  EVALRVAL(rval);

  rval = compare_strands(gth_sa_ref_strand_forward(saA),
                         gth_sa_ref_strand_forward(saB));
  EVALRVAL(rval);

  return 0;
}

/*
  This function is used as ordering for the tree rooted at <rootEST>.
*/

static int compare_duplicate_and_genomic_pos(const GtKeytype dataA,
                                             const GtKeytype dataB,
                                             GT_UNUSED void *cmpinfo)
{
  int rval;

  rval = compare_duplicate(dataA, dataB, NULL);
  EVALRVAL(rval);

  return gth_sa_cmp_genomic_forward(dataA, dataB);
}

/*
  The following function is used as ordering for the tree rooted at <rootlist>.
*/

static int compare_sa(const GtKeytype dataA, const GtKeytype dataB,
                      GT_UNUSED void *cmpinfo)
{
  GthSA *saA = (GthSA*) dataA;
  GthSA *saB = (GthSA*) dataB;
  Spancheck action;
  GtRange rangeA, rangeB;
  int rval = 0;

  gt_assert(cmpinfo == NULL);

  rangeA = gth_sa_range_forward(saA);
  rangeB = gth_sa_range_forward(saB);

  /* save the return value to be able to reuse it */
  rval = gt_str_cmp(gth_sa_ref_id_str(saA), gth_sa_ref_id_str(saB));

  if (!(gt_str_cmp(gth_sa_gen_id_str(saA), gth_sa_gen_id_str(saB))) &&
      (rval == 0) &&
      (gth_sa_gen_strand_forward(saA) == gth_sa_gen_strand_forward(saB)) &&
      (gth_sa_ref_strand_forward(saA) == gth_sa_ref_strand_forward(saB))) {
    action = span_check(saA, saB);
    if (action == DISCARD)
      return 0;
    gt_assert(action != REPLACE);
  }

  /* genomic file number comparison */
  if (gth_sa_gen_file_num(saA) < gth_sa_gen_file_num(saB))
    return -1;
  if (gth_sa_gen_file_num(saA) > gth_sa_gen_file_num(saB))
    return 1;

  /* genomic file numbers are equal, compare ranges */
  if ((rangeA.start < rangeB.start) ||
      ((rangeA.start == rangeB.start) && (rangeA.end > rangeB.end))) {
    /* saA '<' saB */
    return -1;
  }
  else if (rangeA.start == rangeB.start && rangeA.end == rangeB.end) {
    /* XXX: check if this stuff is still necessary (MPI has been removed)
       some additional sorting criteria for the case when the two ranges are the
       same, this helps preventing ``random'' order of spliced alignments in an
       MPI run */
    if (gth_sa_score(saA) < gth_sa_score(saB))
      return -1;
    if (gth_sa_score(saA) > gth_sa_score(saB))
      return 1;

    if (gth_sa_cumlen_scored_exons(saA) <
        gth_sa_cumlen_scored_exons(saB)) {
      return -1;
    }
    if (gth_sa_cumlen_scored_exons(saA) >
        gth_sa_cumlen_scored_exons(saB)) {
      return 1;
    }

    if (gth_sa_gen_strand_forward(saA) < gth_sa_gen_strand_forward(saB))
      return -1;
    if (gth_sa_gen_strand_forward(saA) > gth_sa_gen_strand_forward(saB))
      return 1;

    if (gth_sa_ref_strand_forward(saA) < gth_sa_ref_strand_forward(saB))
      return -1;
    if (gth_sa_ref_strand_forward(saA) > gth_sa_ref_strand_forward(saB))
      return 1;

    /* rval contains the value of a string comparison between the two reference
       ids */
    if (rval < 0)
      return -1;
    else if (rval > 0)
      return  1;
    gt_assert(0);
  }

  /* saA '>' saB */
  return  1;
}

static void insert_alignment(GthSACollection *sa_collection, GthSA *saB)
{
  GthSA *saA = NULL;
  bool nodecreated;

  gt_assert(sa_collection && saB);

  /* insert spliced alignment into tree rooted at <rootlist> */
  saA = (GthSA*) gt_rbt_search(saB, &nodecreated, &sa_collection->rootlist,
                               compare_sa, NULL);
  /* insertion into binary tree succeeded */
  gt_assert(saA && nodecreated);

  /* insert spliced alignment into tree rooted at <rootEST> */
  saA = (GthSA*) gt_rbt_search(saB, &nodecreated, &sa_collection->rootEST,
                               compare_duplicate_and_genomic_pos, NULL);
  /* insertion into binary tree succeeded */
  gt_assert(saA && nodecreated);
}

bool gth_sa_collection_insert_sa(GthSACollection *sa_collection, GthSA *saB,
                                 GthSAFilter *sa_filter, GthStat *stat)
{
  GthSA *saA = NULL, *spliced_alignmentptr = NULL, *satodel;
  Spancheck action;
  GtArray *alignmentstodelete;
  bool discard = false,
       replace = false;
  unsigned long i;

  /* filter */
  if (sa_filter) {
    if (gth_sa_filter_filter_sa(sa_filter, saB)) {
      /* discard saB.
         returning false to indicate that no element has been inserted */
      return false;
    }
  }

  saA = (GthSA*) gt_rbt_find(saB, sa_collection->rootEST, compare_duplicate,
                             NULL);
  if (saA == NULL) {
    /* no alignment with the same ids and strand orientations is in the tree,
       insert saB into both trees. */
    insert_alignment(sa_collection, saB);
  }
  else {
    /* one or more alignments with the same ids and strand orientations exist
       in the tree. check if one of them spans the new alignment. */
    alignmentstodelete = gt_array_new(sizeof (GthSA*));

    /* check the actual one */
    action = span_check(saA, saB);
    switch (action) {
      case NOSPAN:
        /* nothing to do */
        break;
      case DISCARD:
        discard = true;
        break;
      case REPLACE:
        replace = true;
        gt_array_add(alignmentstodelete, saA);
        break;
      default: gt_assert(0);
    }

    /* going to the left */
    spliced_alignmentptr = gt_rbt_previouskey(saA, sa_collection->rootEST,
                                              compare_duplicate_and_genomic_pos,
                                              NULL);

    while ((spliced_alignmentptr) &&
           (!compare_duplicate(saB, spliced_alignmentptr, NULL))) {
      /* spliced_alignmentptr has the same ids and strand orientations:
         check if it spans saB */
      action = span_check(spliced_alignmentptr, saB);
      switch (action) {
        case NOSPAN:
          /* nothing to do */
          break;
        case DISCARD:
          discard = true;
          break;
        case REPLACE:
          replace = true;
          gt_array_add(alignmentstodelete, spliced_alignmentptr);
          break;
        default: gt_assert(0);
      }

      /* retrieving next alignment */
      spliced_alignmentptr = gt_rbt_previouskey(spliced_alignmentptr,
                                                sa_collection->rootEST,
                                              compare_duplicate_and_genomic_pos,
                                                NULL);
    }

    /* going to right */
    spliced_alignmentptr = gt_rbt_nextkey(saA, sa_collection->rootEST,
                                          compare_duplicate_and_genomic_pos,
                                          NULL);
    while ((spliced_alignmentptr) &&
           (!compare_duplicate(saB, spliced_alignmentptr, NULL))) {
      /* spliced_alignmentptr has the same ids and strand orientations:
         check if it spans saB */
      action = span_check(spliced_alignmentptr, saB);
      switch (action) {
        case NOSPAN:
          /* nothing to do */
          break;
        case DISCARD:
          discard = true;
          break;
        case REPLACE:
          replace = true;
          gt_array_add(alignmentstodelete, spliced_alignmentptr);
          break;
        default: gt_assert(0);
      }

      /* retrieving next alignment */
      spliced_alignmentptr = gt_rbt_nextkey(spliced_alignmentptr,
                                            sa_collection->rootEST,
                                            compare_duplicate_and_genomic_pos,
                                            NULL);
    }

    /* DISCARD and REPLACE are not true at the same time */
    gt_assert(!(discard && replace));
    if (discard) {
      gt_assert(saA);
      /* saB will be discarded
         returning false to indicate that no element has been inserted */
      gt_array_delete(alignmentstodelete);
      return false;
    }
    if (replace) {
      /* deleting all alignments which need to be replaced by
         saA */
      for (i = 0; i < gt_array_size(alignmentstodelete); i++) {
        satodel = *(GthSA**) gt_array_get(alignmentstodelete, i);
        (void) gt_rbt_delete(satodel, &sa_collection->rootlist, compare_sa,
                             NULL);
        (void) gt_rbt_delete(satodel, &sa_collection->rootEST,
                             compare_duplicate_and_genomic_pos, NULL);

        /* the ith spliced alignment has been removed from the tree now,
           but the memory still needs to be freed */
        gth_sa_delete(satodel);
      }
    }

    /* free */
    gt_array_delete(alignmentstodelete);

    /* insert saB */
    insert_alignment(sa_collection, saB);
  }

  /* returning true to indicate that an element has been inserted */
  if (stat)
    gth_stat_increment_numofSAs(stat);
  return true;
}

static int storealignmentptr(GtKeytype data, GtRbtVisit which,
                              GT_UNUSED unsigned long depth, void *actinfo)
{
  GthSA *sa = (GthSA*) data;
  GtArray *alignments = (GtArray*) actinfo;
  switch (which) {
    case GT_RBT_PREORDER:
    case GT_RBT_ENDORDER:
      break;
    case GT_RBT_POSTORDER:
    case GT_RBT_LEAF:
      gt_array_add(alignments, sa);
      break;
    default: gt_assert(0);
  }
  return 0;
}

/* The following function traverses the tree of alignments in <sa_collection>
   and returns all alignment pointers in the array. */
static GtArray* sa_collection_get_alignments(const GthSACollection
                                             *sa_collection)
{
  GtArray *alignments = gt_array_new(sizeof (GthSA*));
  int had_err;
  /* traverse the tree */
  had_err = gt_rbt_walk(sa_collection->rootlist, storealignmentptr, alignments);
  gt_assert(!had_err); /* storealignmentptr() is sane */
  return alignments;
}

bool gth_sa_collections_are_equal(const GthSACollection *sa_collectionA,
                                  const GthSACollection *sa_collectionB)
{
  GtArray *alignments_from_A, *alignments_from_B;
  unsigned long i;

  /* compute arrays of SAs from the trees */
  alignments_from_A = sa_collection_get_alignments(sa_collectionA);
  alignments_from_B = sa_collection_get_alignments(sa_collectionB);

  /* check if both arrays have the same size */
  if (gt_array_size(alignments_from_A) != gt_array_size(alignments_from_B)) {
    gt_array_delete(alignments_from_A);
    gt_array_delete(alignments_from_B);
    return false;
  }

  /* both arrays have the same size -> compare the individual SAs */
  for (i = 0; i < gt_array_size(alignments_from_A); i++) {
    if (!gth_sas_are_equal(*(GthSA**) gt_array_get(alignments_from_A, i),
                           *(GthSA**) gt_array_get(alignments_from_B, i))) {
      gt_array_delete(alignments_from_A);
      gt_array_delete(alignments_from_B);
      return false;
    }
  }

  /* free space */
  gt_array_delete(alignments_from_A);
  gt_array_delete(alignments_from_B);

  return true;
}

bool gth_sa_collection_contains_sa(const GthSACollection *sa_collection)
{
  gt_assert(sa_collection);
  return sa_collection->rootlist ? true : false;
}

static void sa_free(void *sa, GT_UNUSED void *unused)
{
  gth_sa_delete(sa);
}

void gth_sa_collection_traverse(const GthSACollection *sa_collection,
                                GthSAVisitor *sa_visitor, GthInput *input)
{
  GthSACollectionIterator *iterator;
  unsigned long num_of_sas = 0;
  GthSA *sa;
  gt_assert(sa_collection && sa_visitor && input);
  gth_sa_visitor_preface(sa_visitor);
  iterator = gth_sa_collection_iterator_new(sa_collection);
  while ((sa = gth_sa_collection_iterator_next(iterator))) {
    gth_input_load_genomic_file(input, gth_sa_gen_file_num(sa), false);
    gth_input_load_reference_file(input, gth_sa_ref_file_num(sa), false);
    gth_sa_visitor_visit_sa(sa_visitor, sa);
    num_of_sas++;
  }
  gth_sa_collection_iterator_delete(iterator);
  gth_sa_visitor_proc_number_of_sas(sa_visitor, num_of_sas);
}

void gth_sa_collection_delete(GthSACollection *sa_collection)
{
  if (!sa_collection) return;
  gt_rbt_destroy(true, sa_free, NULL, sa_collection->rootlist);
  /* for the second tree no free function is needed, because all elements have
     been freed already by the previous redblacktreedestroy() */
  gt_rbt_destroy(false, NULL, NULL, sa_collection->rootEST);
  gt_free(sa_collection);
}

struct GthSACollectionIterator {
  GtArray *alignments;
  unsigned long counter;
};

GthSACollectionIterator* gth_sa_collection_iterator_new(const GthSACollection
                                                        *sa_collection)
{
  GthSACollectionIterator *iterator;
  gt_assert(sa_collection);
  iterator = gt_malloc(sizeof *iterator);
  iterator->alignments = sa_collection_get_alignments(sa_collection);
  iterator->counter = 0;
  return iterator;
}

GthSA* gth_sa_collection_iterator_next(GthSACollectionIterator *iterator)
{
  gt_assert(iterator);
  if (iterator->counter < gt_array_size(iterator->alignments))
    return *(GthSA**) gt_array_get(iterator->alignments, iterator->counter++);
  return NULL;
}

void gth_sa_collection_iterator_delete(GthSACollectionIterator *iterator)
{
  if (!iterator) return;
  gt_array_delete(iterator->alignments);
  gt_free(iterator);
}
