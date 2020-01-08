/*
  Copyright (c) 2005-2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef CONSENSUS_SA_API_H
#define CONSENSUS_SA_API_H

#include "core/range_api.h"
#include "core/strand_api.h"

/* ConsensusSplicedAlignment module */

/* Function to obtain the genomic range from spliced alignment <sa>. */
typedef GtRange  (*GtGetGenomicRangeFunc)(const void *sa);
/* Function to obtain the strand from spliced alignment <sa>. */
typedef GtStrand (*GtGetStrandFunc)(const void *sa);
/* Function to obtain exon ranges from spliced alignment <sa>. */
typedef void   (*GtGetExonsFunc)(GtArray *exon_ranges, const void *sa);
/* Function to process the generated from the consensus spliced alignment
   process. */
typedef void   (*GtProcessSpliceFormFunc)(GtArray *spliced_alignments_in_form,
                                        const void *set_of_sas,
                                        GtUword number_of_sas,
                                        size_t size_of_sa,
                                        void *userdata);

/*
  Construct consensus spliced alignments according to:

  B.J. Haas, A.L. Delcher, S.M. Mount, J.R. Wortman, R.K. Smith Jr,
  L.I. Hannick, R. Maiti, C.M. Ronning, D.B. Rusch, C.D. Town, S.L. Salzberg,
  and O. White. Improving the Arabidopsis genome annotation using maximal
  transcript alignment assemblies. Nucleic Acids Res., 31(19):5654-5666, 2003.

  following the description on page 972 and 973 of the paper:

  G. Gremme, V. Brendel, M.E. Sparks, and S. Kurtz. Engineering a Software Tool
  for Gene Structure Prediction in Higher Organisms. Information and Software
  Technology, 47(15):965-978, 2005.
*/
void gt_consensus_sa(const void *set_of_sas, GtUword number_of_sas,
                     size_t size_of_sa, GtGetGenomicRangeFunc, GtGetStrandFunc,
                     GtGetExonsFunc, GtProcessSpliceFormFunc, void *userdata);

#endif
