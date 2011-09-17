/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHOUTPUT_H
#define GTHOUTPUT_H

/* This file bundles all output related structures and definitions. */

#include "core/error.h"
#include "core/file.h"
#include "gth/gthstrandchar.h"

/* The initial XML indent level */
#define INITIAL_XML_INDENTLEVEL         0

/* This character is shown at the beginnig of comment lines. */
#define COMMENTCHAR '$'

/* This offset is added when printing out sequence positions.
   Usually it equals 1, because in the output we are counting from 1 on.
   Within the program we are always counting from 0 on. */
#define OUTPUTOFFSET 1

#define SHOWGENPOS(FORWARD, TOTALLENGTH, OFFSET, P)\
        (FORWARD)\
        ?  ((P) - (OFFSET) + OUTPUTOFFSET)\
        :  ((TOTALLENGTH) - 1 - ((P) - (OFFSET)) + OUTPUTOFFSET)

#define DELIMITERLINELENGTH     80

/* This width is mainly used for the output of scores */
#define SCORE_WIDTHTYPE         5.3

/* The number of bases in an alignment line */
#define ALIGNMENTLINEWIDTH              60

/* Two times the number of bases in an alignment line */
#define DOUBLEALIGNMENTLINEWIDTH        (2 * ALIGNMENTLINEWIDTH)

/* The source tag used for GFF3 output */
#define GTHSOURCETAG            "gth"

typedef void (*GthShowVerbose)(const char*);
typedef void (*GthShowVerboseVM)(char*);

/* This structure bundles all output related variables. */
typedef struct {
  char *outputfile;                 /* name of output file */
  bool verboseseqs,                 /* show additional sequence information */
                                    /* (for debugging purposes) */
       skipalignmentout,            /* skip the output of spliced alignments */
       showseqnums,                 /* show sequence numbers in output */
       pglgentemplate,              /* show genomic template in PGL lines */
       xmlout,                      /* show output in XML format */
       gff3out,                     /* show output in GFF3 format */
       gff3descranges,              /* use description ranges for GFF3 output */
       gs2out,                      /* output in deprecated GeneSeqer2 format */
       md5ids,                      /* show MD5 fingerprints as sequence IDs */
       comments,                    /* output (additional) comments */
       showeops,                    /* show edit operations after (protein) DP
                                     */
       sortags,                     /* sort AGSs */
       start_codon,                 /* ORF must begin with a start codon */
       final_stop_codon;            /* final ORF must end with a stop codon */
  double sortagswf;                 /* weight factor for the sorting of AGSs */
  unsigned int maxagsnum;           /* the maximum number of AGSs per PGL */
  unsigned long minORFlength,       /* minimum ORF length shown in assembly */
                showintronmaxlen;   /* up to  length an intron is shown
                                       completly, otherwise a part in the middle
                                       is not shown.
                                       if 0 all introns are shown completely */
  int widthforgenpos;               /* width which is used to format genomic
                                       sequence positions (GS2=ifwdth) */
  GthShowVerbose showverbose;       /* function to show status info */
  GthShowVerboseVM showverboseVM;   /* function to show vmatch status info */
  GtFile *outfp;                    /* output file pointer */
} GthOutput;

GthOutput* gthoutput_new(void);
void       gthoutput_delete(GthOutput*);

#endif
