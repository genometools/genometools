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

/* The <GtSeqposClassifier> class.
   Given an annotation (GFF3) of a sequence s, it answers the question whether
   a sequence position i belong to a given type of feature (e.g. CDS) by
   iterating over all positions of s. The 1-based coordinates in the gff3 are
   converted to 0-based sequence positions in the process.
   The features of the given type must be sorted by starting position in the
   GFF3 file, such that when <gt_seqpos_classifier_is_inside_feature()> method
   is called for each, i must increase.
   Currently only a single sequence is supported. To remove this limitation,
   a mapping of sequence regions to sequence positions would be needed. */
typedef struct GtSeqposClassifier GtSeqposClassifier;

/* Creates a new <GtSeqposClassifier> using the annotation in <filename>,
   recognizing the feature type <feature_type>. If <filename> equals <NULL>,
   the data is read from <stdin>. */
GtSeqposClassifier* gt_seqpos_classifier_new(const char *filename,
                                             const char *feature_type);

/* Tests whether <i> is inside a feature in an annotation as defined by
   <seqpos_classifier>. If so, <inside> is set to true. If the end of the
   annotation has been reached, <end_of_annotation> is set to true. In this
   case, the value of <inside> is undefined. Note that <i> must be increasing
   for sequential calls of this function.
   Returns 0 on success, otherwise a different value is returned and <err> is
   set accordingly. */
int                 gt_seqpos_classifier_position_is_inside_feature(
                                          GtSeqposClassifier *seqpos_classifier,
                                          unsigned long i, bool *inside,
                                          bool *end_of_annotation,
                                          GtError *err);

/* Returns the number of features of the feature type modeled in
   <seqpos_classifier> found up to this point, referring to calls of
   <gt_seqpos_classifier_position_is_inside_feature()>. */
unsigned long       gt_seqpos_classifier_nof_features_found(
                                         GtSeqposClassifier *seqpos_classifier);

void                gt_seqpos_classifier_delete(
                                         GtSeqposClassifier *seqpos_classifier);

#endif
