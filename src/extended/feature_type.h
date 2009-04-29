/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef FEATURE_TYPE_H
#define FEATURE_TYPE_H

#include <stdbool.h>

/* Some predefined feature type strings. */
#define gt_ft_CDS                      "CDS"
#define gt_ft_EST_match                "EST_match"
#define gt_ft_LTR_retrotransposon      "LTR_retrotransposon"
#define gt_ft_SNP                      "SNP"
#define gt_ft_TF_binding_site          "TF_binding_site"
#define gt_ft_cDNA_match               "cDNA_match"
#define gt_ft_exon                     "exon"
#define gt_ft_five_prime_UTR           "five_prime_UTR"
#define gt_ft_five_prime_splice_site   "five_prime_splice_site"
#define gt_ft_gene                     "gene"
#define gt_ft_intron                   "intron"
#define gt_ft_inverted_repeat          "inverted_repeat"
#define gt_ft_long_terminal_repeat     "long_terminal_repeat"
#define gt_ft_mRNA                     "mRNA"
#define gt_ft_protein_match            "protein_match"
#define gt_ft_repeat_region            "repeat_region"
#define gt_ft_target_site_duplication  "target_site_duplication"
#define gt_ft_three_prime_UTR          "three_prime_UTR"
#define gt_ft_three_prime_splice_site  "three_prime_splice_site"
#define gt_ft_transcript               "transcript"

#endif
