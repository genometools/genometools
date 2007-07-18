/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONSENSUS_SA_H
#define CONSENSUS_SA_H

#include "libgtcore/range.h"
#include "libgtcore/strand.h"

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

typedef Range  (*GetGenomicRangeFunc)(const void *sa);
typedef Strand (*GetStrandFunc)(const void *sa);
typedef void   (*GetExonsFunc)(Array *exon_ranges, const void *sa, Env*);
typedef void   (*ProcessSpliceFormFunc)(Array *spliced_alignments_in_form,
                                        const void *set_of_sas,
                                        unsigned long number_of_sas,
                                        size_t size_of_sa,
                                        void *userdata, Env*);

void consensus_sa(const void *set_of_sas, unsigned long number_of_sas,
                  size_t size_of_sa, GetGenomicRangeFunc get_genomic_range,
                  GetStrandFunc get_strand, GetExonsFunc get_exons,
                  ProcessSpliceFormFunc process_splice_form, void *userdata,
                  Env *env);

#endif
