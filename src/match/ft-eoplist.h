#ifndef FT_EOPLIST_H
#define FT_EOPLIST_H
#include <stdbool.h>
#include "core/unused_api.h"
#include "core/chardef.h"
#include "match/ft-polish.h"
#include "core/types_api.h"

/* A list of edit operation is representation is represented by the following
   opaque type */

typedef struct GtEoplist GtEoplist;

/* The constructor method */
GtEoplist *gt_eoplist_new(void);

/* the destructor */
void gt_eoplist_delete(GtEoplist *eoplist);

/* reset the list to empty it */
void gt_eoplist_reset(GtEoplist *eoplist);

/* add a match of the given length to the eoplist */
void gt_eoplist_match_add(GtEoplist *eoplist,GtUword length);

/* add a single mismatch to the eoplist */
void gt_eoplist_mismatch_add(GtEoplist *eoplist);

/* add a single deletion to the eoplist */
void gt_eoplist_deletion_add(GtEoplist *eoplist);

/* add a single insertion to the eoplist */
void gt_eoplist_insertion_add(GtEoplist *eoplist);

/* reverse the end of an eoplist beginning with index firstindex. */
void gt_eoplist_reverse_end(GtEoplist *eoplist,GtUword firstindex);

/* obtain length of eoplist */
GtUword gt_eoplist_length(const GtEoplist *eoplist);

/* return number of matches in eoplist */
GtUword gt_eoplist_matches_count(const GtEoplist *eoplist);

/* return number of mismatches in eoplist */
GtUword gt_eoplist_mismatches_count(const GtEoplist *eoplist);

/* return number of deletions in eoplist */
GtUword gt_eoplist_deletions_count(const GtEoplist *eoplist);

/* return number of insertions in eoplist */
GtUword gt_eoplist_insertions_count(const GtEoplist *eoplist);

/* return unit cost of eoplist, which is sum of the number of deletions,
   insertions and mismatches */
GtUword gt_eoplist_unit_cost(const GtEoplist *eoplist);

/* To inspect an edit operation list, one employs the following class */

typedef struct GtEoplistReader GtEoplistReader;

/* verify that the given eoplist represents an alignment of the sequences
   stored by the eoplist->useq and and eoplist->vseq of length
   eoplist->ulen and eoplist->vlen. If eoplist->useq and eoplist->vseq
   have not been set, then checs involing the sequence are not
   perforemed. */

void gt_eoplist_verify(const GtEoplist *eoplist,
                       GtEoplistReader *eoplist_reader,
                       GtUword edist,
                       bool distinguish_mismatch_match);

/* The constructor for the reader with respect to the given eoplist. */
GtEoplistReader *gt_eoplist_reader_new(const GtEoplist *eoplist);

/* An eoplist_reader can be used in two modes when enumerating
   cigar symbols. In the default mode,
   a mismatch is not distinguished from a match. After calling the following
   function, the reader will enumerate cigar symbols in which a mismatch
   has is distinguished from a match. */
void gt_eoplist_reader_distinguish_mismatch_match(
      GtEoplistReader *eoplist_reader);

/* The destructor */
void gt_eoplist_reader_delete(GtEoplistReader *eoplist_reader);

/* There are two ways to access the eoplist:
   1) The first way is to
   enumerate the symbols of the cigar
   string equivalent to the eoplist. We recommend to not store the
   entire cigar string (consisting of all cigar operations), as the eoplist
   is more space efficient. A cigar symbol consists of the type of the
   operation and an iteration value which encodes how many times the operation
   is applied. Hence the following two definitions make sense. */

typedef enum
{
  GtDeletionOp,
  GtInsertionOp,
  GtMismatchOp,
  GtMatchOp,
  GtUndefinedOp
} GtEopType;

/* The following function show an <GtEoptype>-value according to the convention
   of the SAM format. If the flag <distinguish_mismatch_match> is true,
   then for matches and mismatches different character '=' and 'X', respectively
   are returned. Otherwise the same character 'M' is returned.  */
char gt_eoplist_pretty_print(GtEopType eoptype,bool distinguish_mismatch_match);

typedef struct
{
  GtEopType eoptype;
  GtUword iteration;
} GtCigarOp;

/* The following method returns true if and only if there is a next
   cigar operation. If it returns true, then it stores the cigar operation
   in the memory location pointed to by cigar_op. */

bool gt_eoplist_reader_next_cigar(GtCigarOp *cigar_op,
                                  GtEoplistReader *eoplist_reader);

/* The second way of accessing an eoplist is to enumerate the corresponding
   sequence of segments dividing the aligned sequences into non-overlapping
   pairs of segments of u and v. The segments on u have length delta, except
   possibly for the last segment, which only has length delta if
   len_u mod delta == 0. the last segment on u has length len_u mod delta.
   The segments of v have a variable length such that the correspoding
   substring is the substring which aligns with the corresponding segment of
   u. Segments are represented by the following type. */

typedef struct
{
  GtUword aligned_u,
          aligned_v;
} GtEoplistSegment;

/*
  The following function takes an eoplist reader and returns true if and
  only if there are segments. If it returns true, the segment length are
  stored in the structure pointed to by <segment>. */

bool gt_eoplist_reader_next_segment(GtEoplistSegment *segment,
                                    GtEoplistReader *eoplist_reader,
                                    GtUword delta);

/* The characters in the following string correspond to the 0-based integer
   representing the edit operations */

/* for example, the following code snipped outputs the cigarstring
   in a single line to stdout

  GtEoplistReader *eoplist_reader = gt_eoplist_reader_new(eoplist);
  GtCigarOp co;

  while (gt_eoplist_reader_next_cigar(&co,eoplist_reader))
  {
    printf("" GT_WU "%c",co.iteration,
                   gt_eoplist_pretty_print(co.eoptype,
                                           distinguish_mismatch_match));
  }
  printf("\n");
*/

void gt_eoplist_reader_reset(GtEoplistReader *eoplist_reader,
                             const GtEoplist *eoplist);

void gt_eoplist_reader_reset_width(GtEoplistReader *eoplist_reader,
                                   unsigned int width);

void gt_eoplist_set_sequences(GtEoplist *eoplist,
                              const GtUchar *useq,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vlen);

void gt_eoplist_format_generic(FILE *fp,
                               const GtEoplist *eoplist,
                               GtEoplistReader *eoplist_reader,
                               bool distinguish_mismatch_match,
                               const GtUchar *characters,
                               GtUchar wildcardshow);

void gt_eoplist_format_exact(FILE *fp,
                             const GtEoplist *eoplist,
                             GtEoplistReader *eoplist_reader,
                             const GtUchar *characters);

GtUword gt_eoplist_num_segments(const GtEoplist *eoplist,GtUword delta);

double gt_eoplist_segments_entropy(const GtEoplist *eoplist,GtUword delta);

GtEoplist *gt_eoplist_new_from_cigar(const char *cigarstring,GtUword length);

char *gt_eoplist2cigar_string(const GtEoplist *eoplist,
                              bool distinguish_mismatch_match);

void gt_eoplist_seed_display_set(GtEoplist *eoplist);

void gt_eoplist_set_seedoffset(GtEoplist *eoplist,
                               GtUword useedoffset,
                               GtUword seedlen);

void gt_eoplist_show_plain(const GtEoplist *eoplist,FILE *fp);

#ifndef OUTSIDE_OF_GT
void gt_eoplist_polished_ends(GtEoplist *eoplist,
                              const Polishing_info *pol_info,
                              bool withpolcheck);
#endif
#endif
