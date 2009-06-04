/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <string.h>
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/str.h"
#include "core/option.h"
#include "core/undef.h"
#include "core/versionfunc.h"
#include "core/symboldef.h"
#include "match/alphadef.h"
#include "ltrharvest-opt.h"
#include "repeattypes.h"

/*
 The following function shows all options that are set by default or from
 the user on stdout.
*/
void showuserdefinedoptionsandvalues(const LTRharvestoptions *lo)
{
  printf("# user defined options and values:\n");
  if (lo->verbosemode)
  {
    printf("#   verbosemode: On\n");
  }
  else
  {
    printf("#   verbosemode: Off\n");
  }
  printf("#   indexname: %s\n", gt_str_get(lo->str_indexname));
  if (lo->fastaoutput)
  {
    printf("#   outputfile: %s\n", gt_str_get(lo->str_fastaoutputfilename));
  }
  if (lo->fastaoutputinnerregion)
  {
    printf("#   outputfile inner region: %s\n",
        gt_str_get(lo->str_fastaoutputfilenameinnerregion));
  }
  if (lo->gff3output)
  {
    printf("#   outputfile gff3 format: %s\n",
           gt_str_get(lo->str_gff3filename));
  }
  printf("#   xdropbelowscore: %d\n", lo->xdropbelowscore);
  printf("#   similaritythreshold: %.2f\n", lo->similaritythreshold);
  printf("#   minseedlength: %lu\n", lo->minseedlength);
  printf("#   matchscore: %d\n", lo->arbitscores.mat);
  printf("#   mismatchscore: %d\n", lo->arbitscores.mis);
  printf("#   insertionscore: %d\n", lo->arbitscores.ins);
  printf("#   deletionscore: %d\n", lo->arbitscores.del);
  printf("#   minLTRlength: %lu\n",  lo->repeatinfo.lmin);
  printf("#   maxLTRlength: %lu\n",  lo->repeatinfo.lmax);
  printf("#   minLTRdistance: %lu\n",  lo->repeatinfo.dmin);
  printf("#   maxLTRdistance: %lu\n",  lo->repeatinfo.dmax);
  if (lo->nooverlapallowed)
  {
    printf("#   overlaps: no\n");
  }
  else
  {
    if (lo->bestofoverlap)
    {
      printf("#   overlaps: best\n");
    }
    else
    {
      printf("#   overlaps: all\n");
    }
  }
  printf("#   minTSDlength: %u\n",  lo->minlengthTSD);
  printf("#   maxTSDlength: %u\n",  lo->maxlengthTSD);
  printf("#   palindromic motif: %s\n", gt_str_get(lo->motif.str_motif));
  printf("#   motifmismatchesallowed: %u\n", lo->motif.allowedmismatches);
  printf("#   vicinity: " FormatSeqpos " nt\n",
          PRINTSeqposcast(lo->vicinityforcorrectboundaries));
  if (lo->repeatinfo.ltrsearchseqrange.start != 0 ||
      lo->repeatinfo.ltrsearchseqrange.end != 0)
  {
    printf("# ltrsearchseqrange=(%lu,%lu)\n",
          PRINTSeqposcast(lo->repeatinfo.ltrsearchseqrange.start),
          PRINTSeqposcast(lo->repeatinfo.ltrsearchseqrange.end));
  }
}

/*
 This function prints the arguments from argv on standard output.
 */
void printargsline(const char **argv, int argc)
{
  int i;

  printf("# args=");
  for (i=1; i<argc; i++)
  {
    printf("%s",argv[i]);
    if (i == (argc-1))
    {
      printf("\n");
    } else
    {
      printf(" ");
    }
  }
}

/* test the motif and encode the characters by using alpha */
int testmotifandencodemotif (Motif *motif, const Encodedsequence *encseq,
                             GtError *err)
{
  const GtUchar *symbolmap;
  GtUchar c_tab[UCHAR_MAX+1];
  unsigned int i;

  symbolmap = getencseqAlphabetsymbolmap(encseq);
  if ( symbolmap[(unsigned int)motif->firstleft] == (GtUchar) UNDEFCHAR)
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondleft] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->firstright] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstright);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondright] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondright);
    return -1;
  }

  for (i=0; i<=(unsigned int) UCHAR_MAX; i++)
  {
    c_tab[i] = (GtUchar) UNDEFCHAR;
  }
  /* define complementary symbols */
  c_tab[symbolmap['a']] = symbolmap['t'];
  c_tab[symbolmap['c']] = symbolmap['g'];
  c_tab[symbolmap['g']] = symbolmap['c'];
  c_tab[symbolmap['t']] = symbolmap['a'];

  /* if motif is not palindromic */
  if ( (c_tab[symbolmap[(unsigned int)motif->firstleft]] !=
       c_tab[c_tab[symbolmap[(unsigned int)motif->secondright]]])
           ||
      (c_tab[symbolmap[(unsigned int)motif->secondleft]] !=
       c_tab[c_tab[symbolmap[(unsigned int)motif->firstright]]]) )
  {
    gt_error_set(err, "Illegal motif, motif not palindromic");
    return -1;
  }

  /* encode the symbols */
  motif->firstleft = symbolmap[(unsigned int)motif->firstleft];
  motif->secondleft = symbolmap[(unsigned int)motif->secondleft];
  motif->firstright = symbolmap[(unsigned int)motif->firstright];
  motif->secondright = symbolmap[(unsigned int)motif->secondright];

  return 0;
}

static OPrval parse_options(int *parsed_args,
                            LTRharvestoptions *lo,
                            int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionindex,
         *optionltrsearchseqrange,
         *optionseed,
         *optionminlenltr,
         *optionmaxlenltr,
         *optionmindistltr,
         *optionmaxdistltr,
         *optionmintsd,
         *optionmaxtsd,
         *optionsimilar,
         *optionmotif,
         *optionmotifmis,
         *optionvic,
         *optionoverlaps,
         *optionxdrop,
         *optionmat,
         *optionmis,
         *optionins,
         *optiondel,
         *optionv,
         *optionlongoutput,
         *optionout,
         *optionoutinner,
         *optiongff3;
  OPrval oprval;
  GtRange default_ltrsearchseqrange = {0,0};
  unsigned int vicinityforcorrectboundaries;

  static const char *overlaps[] = {
    "best", /* the default */
    "no",
    "all",
    NULL
  };

  gt_error_check(err);
  op = gt_option_parser_new("[option ...] -index filenameindex",
                         "Predict LTR retrotransposons.");

  /* -index */
  lo->str_indexname = gt_str_new();
  optionindex = gt_option_new_string("index",
                             "specify the name of the enhanced suffix "
                             "array index (mandatory)",
                             lo->str_indexname, NULL);
  gt_option_is_mandatory(optionindex);
  gt_option_parser_add_option(op, optionindex);

  /* -range */
  optionltrsearchseqrange
    = gt_option_new_range("range",
                          "specify sequence range in which LTRs are searched",
                          &lo->repeatinfo.ltrsearchseqrange,
                          &default_ltrsearchseqrange);
  gt_option_parser_add_option(op, optionltrsearchseqrange);

  /* -seed */
  optionseed = gt_option_new_ulong_min("seed",
                               "specify minimum seed length for"
                               " exact repeats",
                               &lo->minseedlength,
                               30UL,
                               1UL);
  gt_option_parser_add_option(op, optionseed);

  /* -minlenltr */
  optionminlenltr = gt_option_new_ulong_min_max("minlenltr",
                               "specify minimum length for each LTR",
                               &lo->repeatinfo.lmin,
                               100UL,
                               1UL,
                               UNDEF_ULONG);
  gt_option_parser_add_option(op, optionminlenltr);

  /* -maxlenltr */
  optionmaxlenltr = gt_option_new_ulong_min_max("maxlenltr",
                               "specify maximum length for each LTR",
                               &lo->repeatinfo.lmax,
                               1000UL,
                               1UL,
                               UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxlenltr);

  /* -mindistltr */
  optionmindistltr = gt_option_new_ulong_min_max("mindistltr",
                               "specify minimum distance of "
                               "LTR startpositions",
                               &lo->repeatinfo.dmin,
                               1000UL,
                               1UL,
                               UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmindistltr);

  /* -maxdistltr */
  optionmaxdistltr = gt_option_new_ulong_min_max("maxdistltr",
                               "specify maximum distance of "
                               "LTR startpositions",
                               &lo->repeatinfo.dmax,
                               15000UL,
                               1UL,
                               UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxdistltr);

  /* -similar */
  optionsimilar = gt_option_new_double_min_max("similar",
                               "specify similaritythreshold in "
                               "range [1..100%]",
                               &lo->similaritythreshold,
                               (double) 85.0,
                               (double) 0.0,
                               100.0);
  gt_option_parser_add_option(op, optionsimilar);

  /* -mintsd */
  optionmintsd = gt_option_new_uint_min_max("mintsd",
                              "specify minimum length for each TSD",
                               &lo->minlengthTSD,
                               4U,
                               0,
                               UNDEF_UINT);
  gt_option_parser_add_option(op, optionmintsd);

  /* -maxtsd */
  optionmaxtsd = gt_option_new_uint_min_max("maxtsd",
                              "specify maximum length for each TSD",
                               &lo->maxlengthTSD,
                               20U,
                               0,
                               UNDEF_UINT);
  gt_option_parser_add_option(op, optionmaxtsd);

  /* -motif */
  /* characters will be tranformed later
     into characters from virtualtree alphabet */
  lo->motif.firstleft   = (GtUchar) 't';
  lo->motif.secondleft  = (GtUchar) 'g';
  lo->motif.firstright  = (GtUchar) 'c';
  lo->motif.secondright = (GtUchar) 'a';
  lo->motif.str_motif = gt_str_new();
  optionmotif = gt_option_new_string("motif",
                             "specify 2 nucleotides startmotif + "
                             "2 nucleotides endmotif: ****",
                             lo->motif.str_motif, NULL);
  gt_option_parser_add_option(op, optionmotif);

  /* -motifmis */
  optionmotifmis = gt_option_new_uint_min_max("motifmis",
                             "specify maximum number of "
                             "mismatches in motif [0,3]",
                             &lo->motif.allowedmismatches,
                             4U,
                             0,
                             3U);
  gt_option_parser_add_option(op, optionmotifmis);

  /* -vic */
  optionvic = gt_option_new_uint_min_max("vic",
                        "specify the number of nucleotides (to the left and "
                        "to the right) that will be searched "
                        "for TSDs and/or motifs around 5' and 3' boundary "
                        "of predicted LTR retrotransposons",
                        &vicinityforcorrectboundaries,
                        60U,
                        1U,
                        500U);
  gt_option_parser_add_option(op, optionvic);

  /* -overlaps */
  lo->str_overlaps = gt_str_new();
  optionoverlaps = gt_option_new_choice("overlaps",
               "specify no|best|all",
               lo->str_overlaps,
               overlaps[0],
               overlaps);
  gt_option_parser_add_option(op, optionoverlaps);

  /* -xdrop */
  optionxdrop = gt_option_new_int_min("xdrop",
                        "specify xdropbelowscore for extension-alignment",
                        &lo->xdropbelowscore,
                        (int)5,
                        (int)0);
  gt_option_parser_add_option(op, optionxdrop);

  /* -mat */
  lo->arbitscores.gcd  = (int) 1;      /* set only for initialization,
                                        do not change! */
  optionmat = gt_option_new_int_min("mat",
                        "specify matchscore for extension-alignment",
                        &lo->arbitscores.mat,
                        (int)2,
                        (int)1);
  gt_option_parser_add_option(op, optionmat);

  /* -mis */
  optionmis = gt_option_new_int_max("mis",
                        "specify mismatchscore for extension-alignment",
                        &lo->arbitscores.mis,
                        (int)-2,
                        (int)-1);
  gt_option_parser_add_option(op, optionmis);

  /* -ins */
  optionins = gt_option_new_int_max("ins",
                        "specify insertionscore for extension-alignment",
                        &lo->arbitscores.ins,
                        (int)-3,
                        (int)-1);
  gt_option_parser_add_option(op, optionins);

  /* -del */
  optiondel = gt_option_new_int_max("del",
                        "specify deletionscore for extension-alignment",
                        &lo->arbitscores.del,
                        (int)-3,
                        (int)-1);
  gt_option_parser_add_option(op, optiondel);

  /* -v */
  optionv = gt_option_new_bool("v",
                           "verbose mode",
                           &lo->verbosemode,
                           false);
  gt_option_parser_add_option(op, optionv);

  /* -longoutput */
  optionlongoutput = gt_option_new_bool("longoutput",
                           "additional motif/TSD output",
                           &lo->longoutput,
                           false);
  gt_option_parser_add_option(op, optionlongoutput);

  /* -out */
  lo->fastaoutput = false;      /* by default no FASTA output */
  lo->str_fastaoutputfilename = gt_str_new();
  optionout = gt_option_new_string("out",
                             "specify FASTA outputfilename",
                             lo->str_fastaoutputfilename, NULL);
  gt_option_parser_add_option(op, optionout);

  /* -outinner */
  lo->fastaoutputinnerregion = false;
  lo->str_fastaoutputfilenameinnerregion = gt_str_new();
  optionoutinner = gt_option_new_string("outinner",
                             "specify FASTA outputfilename for inner regions",
                             lo->str_fastaoutputfilenameinnerregion, NULL);
  gt_option_parser_add_option(op, optionoutinner);

  /* -gff3 */
  lo->gff3output = false;       /* by default no gff3 output */
  lo->str_gff3filename = gt_str_new();
  optiongff3 = gt_option_new_string("gff3",
                             "specify GFF3 outputfilename",
                             lo->str_gff3filename, NULL);
  gt_option_parser_add_option(op, optiongff3);

  /* implications */
  gt_option_imply(optionmaxtsd, optionmintsd);
  gt_option_imply(optionmotifmis, optionmotif);

  gt_option_imply_either_2(optionlongoutput, optionmintsd, optionmotif);

  gt_option_parser_refer_to_manual(op);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  lo->vicinityforcorrectboundaries = (Seqpos) vicinityforcorrectboundaries;
  if (oprval == OPTIONPARSER_OK)
  {
    if (lo->repeatinfo.lmin > lo->repeatinfo.lmax)
    {
      gt_error_set(err,"argument of -minlenltr is greater than argument of"
          " -maxlenltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if (lo->repeatinfo.dmin > lo->repeatinfo.dmax)
    {
      gt_error_set(err,
          "argument of -mindistltr is greater than argument of -maxdistltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if (lo->repeatinfo.lmax > lo->repeatinfo.dmin)
    {
      gt_error_set(err,"argument of -maxlenltr is greater than argument of"
                    " -mindistltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if (lo->minlengthTSD > lo->maxlengthTSD)
    {
      gt_error_set(err,
          "argument of -mintsd is greater than argument of -maxtsd");
      oprval = OPTIONPARSER_ERROR;
    }

    /* If option motif is set,
       store characters, transform them later */
    if (gt_option_is_set(optionmotif))
    {
      if (gt_str_length(lo->motif.str_motif) != 4UL)
      {
        gt_error_set(err,
            "argument of -motif has not exactly 4 characters");
        oprval = OPTIONPARSER_ERROR;
      }
      lo->motif.firstleft = (GtUchar)  gt_str_get(lo->motif.str_motif)[0];
      lo->motif.secondleft = (GtUchar)  gt_str_get(lo->motif.str_motif)[1];
      lo->motif.firstright = (GtUchar)  gt_str_get(lo->motif.str_motif)[2];
      lo->motif.secondright = (GtUchar)  gt_str_get(lo->motif.str_motif)[3];
      /* default if motif specified */
      if (!gt_option_is_set(optionmotifmis))
      {
        lo->motif.allowedmismatches = 0;
      }
    }

    /* If option overlaps is set */
    if (gt_option_is_set(optionoverlaps))
    {
      if (strcmp(gt_str_get(lo->str_overlaps), "no") == 0)
      {
        lo->bestofoverlap = false;
        lo->nooverlapallowed = true;
      }
      else if (strcmp(gt_str_get(lo->str_overlaps), "best") == 0 )
      {
        lo->bestofoverlap = true;
        lo->nooverlapallowed = false;
      }
      else if (strcmp(gt_str_get(lo->str_overlaps), "all") == 0 )
      {
        lo->bestofoverlap = false;
        lo->nooverlapallowed = false;
      }
      else
      {
        gt_assert(0); /* cannot happen */
      }
    }
    else
    {
      /* default is "best" */
      lo->bestofoverlap = true;     /* take best prediction
                                       if overlap occurs, default */
      lo->nooverlapallowed = false; /* overlapping predictions (not)allowed*/
    }

    /* if FASTA output is set */
    if (gt_option_is_set(optionout))
    {
      lo->fastaoutput = true;
    }

    /* if FASTA output inner region is set */
    if (gt_option_is_set(optionoutinner))
    {
      lo->fastaoutputinnerregion = true;
    }

    /* if GFF3 output is set */
    if (gt_option_is_set(optiongff3))
    {
      lo->gff3output = true;
    }
  }

  gt_option_parser_delete(op);
  return oprval;
}

void wrapltrharvestoptions(LTRharvestoptions *lo)
{
  /* no checking if error occurs, since errors have been output before */
  gt_str_delete(lo->str_indexname);
  gt_str_delete(lo->str_fastaoutputfilename);
  gt_str_delete(lo->str_fastaoutputfilenameinnerregion);
  gt_str_delete(lo->str_gff3filename);
  gt_str_delete(lo->str_overlaps);
  gt_str_delete(lo->motif.str_motif);
}

int ltrharvestoptions(LTRharvestoptions *lo, int argc, const char **argv,
                      GtError *err)
{
  int parsed_args;
  OPrval rval;

  gt_error_check(err);

  /** init LTRharvestoptions lo **/
  rval = parse_options(&parsed_args, lo, argc, argv, err);
  if (rval == OPTIONPARSER_OK)
  {
    if (parsed_args != argc)
    {
      gt_error_set(err, "Listing of options and arguments is not correct.");
      rval = OPTIONPARSER_ERROR;
    }
  }
  return (rval == OPTIONPARSER_OK) ? 0: -1;
}
