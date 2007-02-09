/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONSENSUS_SA_H
#define CONSENSUS_SA_H

#include "array.h"
#include "log.h"
#include "range.h"
#include "strand.h"

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

typedef Range (*Get_genomic_range_func)(const void *sa);
typedef Strand (*Get_strand_func)(const void *sa);
typedef void (*Get_exons_func)(Array *exon_ranges, const void *sa);
typedef void (*Process_splice_form_func)(Array *spliced_alignments_in_form,
                                         const void *set_of_sas,
                                         unsigned long number_of_sas,
                                         size_t size_of_sa,
                                         void *userdata);

void consensus_sa(const void *set_of_sas,
                  unsigned long number_of_sas,
                  size_t size_of_sa,
                  Get_genomic_range_func get_genomic_range,
                  Get_strand_func get_strand,
                  Get_exons_func get_exons,
                  Process_splice_form_func process_splice_form,
                  void *userdata,
                  Log *l);

#endif
