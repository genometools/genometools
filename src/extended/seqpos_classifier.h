/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQPOS_CLASSIFIER_H
#define SEQPOS_CLASSIFIER_H

typedef struct GtSeqposClassifier GtSeqposClassifier;

/* Given an annotation (gff3) of a sequence s,
 * by iterating over all positions of the s allows to
 * answer the question: does the sequence position i
 * belong to a given type of feature (e.g. CDS)?
 *
 * The 1-based coordinates in the gff3 are converted
 * to 0-based sequence positions.
 *
 * The requirements/limitations are:
 * - the features of the given type must be available
 *   sorted by starting position in the gff3 file
 * - the _is_inside_feature(i) method is called for each
 *   position i starting from 0, sequentially
 * - a single type of feature is considered, specified
 *   when calling the _new() method
 * - currently a single sequence is supported (to remove this
 *   limitation, a mapping of sequence regions to sequence positions
 *   would be needed)
 *
 */

/*
 * filename: name of the (sorted) gff3 annotation file
 * feature_type: type of feature to use, e.g. gt_ft_CDS
 *
 */

GtSeqposClassifier *gt_seqpos_classifier_new(const char *filename,
    const char *feature_type);

/*
 * i: the position to test (the function must be called sequentially
 *    for each value of i from 0 to the last position to test)
 *
 * inside: this is set to true if the position is inside a feature
 *         of the feature type specified in the _new function (false otherwise)
 *
 * end_of_annotation: this is set to true if the annotation has been
 *                    completely parsed (false otherwise)
 *                    (in this case, the value of inside is unspecified)
 *
 * return value: 0 on success, another value if an error occurred
 *               in which case <err> is also set
 *
 */
int gt_seqpos_classifier_position_is_inside_feature(
    GtSeqposClassifier *seqpos_classifier, unsigned long i, bool *inside,
    bool *end_of_annotation, GtError *err);

/* returns the number of feature nodes of the specifed feature type found */
unsigned long gt_seqpos_classifier_nof_features_found(
    GtSeqposClassifier *seqpos_classifier);

void gt_seqpos_classifier_delete(GtSeqposClassifier *seqpos_classifier);

#endif
