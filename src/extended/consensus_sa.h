/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef CONSENSUS_SA_H
#define CONSENSUS_SA_H

#include "core/range.h"
#include "core/strand_api.h"

/*
  This module implements the method to construct consensus spliced alignments

  B.J. Haas, A.L. Delcher, S.M. Mount, J.R. Wortman, R.K. Smith Jr,
  L.I. Hannick, R. Maiti, C.M. Ronning, D.B. Rusch, C.D. Town, S.L. Salzberg,
  and O. White. Improving the Arabidopsis genome annotation using maximal
  transcript alignment assemblies. Nucleic Acids Res., 31(19):5654-5666, 2003.

  following the description on page 972 and 973 of the paper

  G. Gremme, V. Brendel, M.E. Sparks, and S. Kurtz. Engineering a Software Tool
  for Gene Structure Prediction in Higher Organisms. Information and Software
  Technology, 47(15):965-978, 2005.
*/

typedef GtRange  (*GetGenomicRangeFunc)(const void *sa);
typedef GtStrand (*GetStrandFunc)(const void *sa);
typedef void   (*GetExonsFunc)(GtArray *exon_ranges, const void *sa);
typedef void   (*ProcessSpliceFormFunc)(GtArray *spliced_alignments_in_form,
                                        const void *set_of_sas,
                                        unsigned long number_of_sas,
                                        size_t size_of_sa,
                                        void *userdata);

void gt_consensus_sa(const void *set_of_sas, unsigned long number_of_sas,
                     size_t size_of_sa, GetGenomicRangeFunc, GetStrandFunc,
                     GetExonsFunc, ProcessSpliceFormFunc, void *userdata);

#endif
