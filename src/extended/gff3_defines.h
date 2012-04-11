/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GFF3_DEFINES_H
#define GFF3_DEFINES_H

/* version */
#define GT_GFF_VERSION             3

/* pragmas */
#define GT_GFF_VERSION_PREFIX      "##gff-version"
#define GT_GFF_FASTA_DIRECTIVE     "##FASTA"
#define GT_GFF_SEQUENCE_REGION     "##sequence-region"
#define GT_GFF_SPECIES             "##species"
#define GT_GFF_FEATURE_ONTOLOGY    "##feature-ontology"
#define GT_GFF_ATTRIBUTE_ONTOLOGY  "##attribute-ontology"
#define GT_GFF_SOURCE_ONTOLOGY     "##source-ontology"
#define GT_GFF_NCBI_TAXONOMY_URI   "##NCBI_Taxonomy_URI"
#define GT_GFF_GENOME_BUILD        "##genome-build"
#define GT_GFF_TERMINATOR          "###"

/* predfined attributes */
#define GT_GFF_ID              "ID"
#define GT_GFF_NAME            "Name"
#define GT_GFF_PARENT          "Parent"
#define GT_GFF_TARGET          "Target"

#endif
