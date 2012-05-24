/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef GFF3_DEFINES_H
#define GFF3_DEFINES_H

#define GT_GFF_META_CHARS          "##"

/* version */
#define GT_GFF_VERSION             3
#define GT_GFF_VERSION_STRING      "3"
#define GT_GVF_VERSION             "1.06"

/* pragmas */
#define GT_GFF_VERSION_DIRECTIVE   "gff-version"
#define GT_GFF_VERSION_PREFIX      GT_GFF_META_CHARS GT_GFF_VERSION_DIRECTIVE
#define GT_GFF_FASTA_DIRECTIVE     "##FASTA"
#define GT_GFF_SEQUENCE_REGION     "##sequence-region"
#define GT_GFF_SPECIES             "##species"
#define GT_GFF_FEATURE_ONTOLOGY    "##feature-ontology"
#define GT_GFF_ATTRIBUTE_ONTOLOGY  "##attribute-ontology"
#define GT_GFF_SOURCE_ONTOLOGY     "##source-ontology"
#define GT_GFF_NCBI_TAXONOMY_URI   "##NCBI_Taxonomy_URI"
#define GT_GFF_GENOME_BUILD        "##genome-build"
#define GT_GFF_TERMINATOR          "###"
#define GT_GVF_VERSION_DIRECTIVE   "gvf-version"
#define GT_GVF_VERSION_PREFIX      GT_GFF_META_CHARS GT_GVF_VERSION_DIRECTIVE
#define GT_GVF_REFERENCE_FASTA     "##reference-fasta"
#define GT_GVF_FEATURE_GFF3        "##feature-gff3"
#define GT_GVF_FILE_VERSION        "##file-version"
#define GT_GVF_FILE_DATE           "##file-date"
#define GT_GVF_INDIVIDUAL_ID       "##individual-id"
#define GT_GVF_POPULATION          "##population"
#define GT_GVF_SEX                 "##sex"
#define GT_GVF_TECHNOLOGY_PLATFORM "##technology-platform"
#define GT_GVF_TECHNOLOGY_PLATFORM_CLASS \
                                   "##technology-platform-class"
#define GT_GVF_TECHNOLOGY_PLATFORM_NAME \
                                   "##technology-platform-name"
#define GT_GVF_TECHNOLOGY_PLATFORM_VERSION \
                                   "##technology-platform-version"
#define GT_GVF_TECHNOLOGY_PLATFORM_MACHINE_ID \
                                   "##technology-platform-machine-id"
#define GT_GVF_TECHNOLOGY_PLATFORM_READ_LENGTH \
                                   "##technology-platform-read-length"
#define GT_GVF_TECHNOLOGY_PLATFORM_READ_TYPE \
                                   "##technology-platform-read-type"
#define GT_GVF_TECHNOLOGY_PLATFORM_READ_PAIR_SPAN \
                                   "##technology-platform-read-pair-span"
#define GT_GVF_TECHNOLOGY_PLATFORM_AVERAGE_COVERAGE \
                                   "##technology-platform-average-coverage"
#define GT_GVF_SEQUENCING_SCOPE    "##sequencing-scope"
#define GT_GVF_CAPTURE_METHOD      "##capture-method"
#define GT_GVF_CAPTURE_REGIONS     "##capture-regions"
#define GT_GVF_SEQUENCE_ALIGNMENT  "##sequence-alignment"
#define GT_GVF_VARIANT_CALLING     "##variant-calling"
#define GT_GVF_SAMPLE_DESCRIPTION  "##sample-description"
#define GT_GVF_GENOMIC_SOURCE      "##genomic-source"
#define GT_GVF_MULTI_INDIVIDUAL    "##multi-individual"
#define GT_GVF_DATA_SOURCE         "##data-source"
#define GT_GVF_SCORE_METHOD        "##score-method"
#define GT_GVF_SOURCE_METHOD       "##source-method"
#define GT_GVF_ATTRIBUTE_METHOD    "##attribute-method"
#define GT_GVF_PHENOTYPE_DESCRIPTION \
                                   "##phenotype-description"
#define GT_GVF_PHASED_GENOTYPES    "##phased-genotypes"

/* predefined attributes */
#define GT_GFF_ID                 "ID"
#define GT_GFF_NAME               "Name"
#define GT_GFF_ALIAS              "Alias"
#define GT_GFF_PARENT             "Parent"
#define GT_GFF_TARGET             "Target"
#define GT_GFF_GAP                "Gap"
#define GT_GFF_DERIVES_FROM       "Derives_from"
#define GT_GFF_NOTE               "Note"
#define GT_GFF_DBXREF             "Dbxref"
#define GT_GFF_ONTOLOGY_TERM      "Ontology_term"
#define GT_GFF_IS_CIRCULAR        "Is_circular"
#define GT_GVF_GENOTYPE           "Genotype"
#define GT_GVF_REFERENCE_SEQ      "Reference_seq"
#define GT_GVF_VARIANT_SEQ        "Variant_seq"
#define GT_GVF_VARIANT_FREQ       "Variant_freq"
#define GT_GVF_VARIANT_EFFECT     "Variant_effect"
#define GT_GVF_VARIANT_READS      "Variant_reads"
#define GT_GVF_TOTAL_READS        "Total_reads"
#define GT_GVF_PHASED             "Phased"
#define GT_GVF_START_RANGE        "Start_range"
#define GT_GVF_END_RANGE          "End_range"
#define GT_GVF_INDIVIDUAL         "Individual"
#define GT_GVF_REFERENCE_CODON    "Reference_codon"
#define GT_GVF_VARIANT_CODON      "Variant_codon"
#define GT_GVF_REFERENCE_AA       "Reference_aa"
#define GT_GVF_VARIANT_AA         "Variant_aa"
#define GT_GVF_BREAKPOINT_DETAIL  "Breakpoint_detail"
#define GT_GVF_SEQUENCE_CONTEXT   "Sequence_context"
#define GT_GVF_ZYGOSITY           "Zygosity"

#endif
