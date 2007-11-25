/*
  Copyright (c) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
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
#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "libgtcore/env.h"
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/undef.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/symboldef.h"
#include "libgtmatch/alphadef.h"
#include "ltrharvest-opt.h"
#include "repeattypes.h"

/*
 The following function shows all options that are set by default or from
 the user on stdout.
*/
void showuserdefinedoptionsandvalues(LTRharvestoptions *lo)
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
  printf("#   indexname: %s\n", str_get(lo->str_indexname));
  if (lo->fastaoutput)
  {
    printf("#   outputfile: %s\n", str_get(lo->str_fastaoutputfilename));
  }
  if (lo->fastaoutputinnerregion)
  {
    printf("#   outputfile inner region: %s\n",
        str_get(lo->str_fastaoutputfilenameinnerregion));
  }
  if (lo->gff3output)
  {
    printf("#   outputfile gff3 format: %s\n", str_get(lo->str_gff3filename));
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
  printf("#   palindromic motif: %s\n", str_get(lo->motif.str_motif));
  printf("#   motifmismatchesallowed: %u\n", lo->motif.allowedmismatches);
  printf("#   vicinity: " FormatSeqpos " nt\n",
          PRINTSeqposcast(lo->vicinityforcorrectboundaries));
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
int testmotifandencodemotif (Motif *motif, const Alphabet *alpha, Env *env)
{
  const Uchar *symbolmap;
  Uchar c_tab[UCHAR_MAX+1];
  unsigned int i;

  symbolmap = getsymbolmapAlphabet(alpha);
  if ( symbolmap[(unsigned int)motif->firstleft] == (Uchar) UNDEFCHAR)
  {
    env_error_set(env,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondleft] == (Uchar) UNDEFCHAR )
  {
    env_error_set(env,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->firstright] == (Uchar) UNDEFCHAR )
  {
    env_error_set(env,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstright);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondright] == (Uchar) UNDEFCHAR )
  {
    env_error_set(env,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondright);
    return -1;
  }

  for (i=0; i<=(unsigned int) UCHAR_MAX; i++)
  {
    c_tab[i] = (Uchar) UNDEFCHAR;
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
    env_error_set(env, "Illegal motif, motif not palindromic");
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
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *optionindex,
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
  unsigned int vicinityforcorrectboundaries;

  static const char *overlaps[] = {
    "best", /* the default */
    "no",
    "all",
    NULL
  };

  env_error_check(env);
  op = option_parser_new("[option ...] -index filenameindex",
                         "Predict LTR retrotransposons.", env);
  option_parser_set_mailaddress(op,"<dellinghaus@zbh.uni-hamburg.de>");

  /* -index */
  lo->str_indexname = str_new(env);
  optionindex = option_new_string("index",
                             "specify the name of the enhanced suffix "
                             "array index (mandatory)",
                             lo->str_indexname, NULL, env);
  option_is_mandatory(optionindex);
  option_parser_add_option(op, optionindex, env);

  /* -seed */
  optionseed = option_new_ulong_min("seed",
                               "specify minimum seed length for"
                               " exact repeats",
                               &lo->minseedlength,
                               (unsigned long) 30,
                               (unsigned long) 1,
                               env);
  option_parser_add_option(op, optionseed, env);

  /* -minlenltr */
  optionminlenltr = option_new_ulong_min_max("minlenltr",
                               "specify minimum length for each LTR",
                               &lo->repeatinfo.lmin,
                               (unsigned long) 100,
                               (unsigned long) 1,
                               UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionminlenltr, env);

  /* -maxlenltr */
  optionmaxlenltr = option_new_ulong_min_max("maxlenltr",
                               "specify maximum length for each LTR",
                               &lo->repeatinfo.lmax,
                               (unsigned long) 1000,
                               (unsigned long) 1,
                               UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionmaxlenltr, env);

  /* -mindistltr */
  optionmindistltr = option_new_ulong_min_max("mindistltr",
                               "specify minimum distance of "
                               "LTR startpositions",
                               &lo->repeatinfo.dmin,
                               (unsigned long) 1000,
                               (unsigned long) 1,
                               UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionmindistltr, env);

  /* -maxdistltr */
  optionmaxdistltr = option_new_ulong_min_max("maxdistltr",
                               "specify maximum distance of "
                               "LTR startpositions",
                               &lo->repeatinfo.dmax,
                               (unsigned long) 15000,
                               (unsigned long) 1,
                               UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionmaxdistltr, env);

  /* -similar */
  optionsimilar = option_new_double_min_max("similar",
                               "specify similaritythreshold in "
                               "range [1..100%]",
                               &lo->similaritythreshold,
                               (double) 85.0,
                               (double) 0.0,
                               100.0,
                               env);
  option_parser_add_option(op, optionsimilar, env);

  /* -mintsd */
  optionmintsd = option_new_uint_min_max("mintsd",
                              "specify minimum length for each TSD",
                               &lo->minlengthTSD,
                               0,
                               0,
                               UNDEF_UINT,
                               env);
  option_parser_add_option(op, optionmintsd, env);

  /* -maxtsd */
  optionmaxtsd = option_new_uint_min_max("maxtsd",
                              "specify maximum length for each TSD",
                               &lo->maxlengthTSD,
                               20U,
                               0,
                               (unsigned int) UNDEF_INT,
                               env);
  option_parser_add_option(op, optionmaxtsd, env);

  /* -motif */
  /* characters will be tranformed later
     into characters from virtualtree alphabet */
  lo->motif.firstleft   = (Uchar) 't';
  lo->motif.secondleft  = (Uchar) 'g';
  lo->motif.firstright  = (Uchar) 'c';
  lo->motif.secondright = (Uchar) 'a';
  lo->motif.str_motif = str_new(env);
  optionmotif = option_new_string("motif",
                             "specify 2 nucleotides startmotif + "
                             "2 nucleotides endmotif: ****",
                             lo->motif.str_motif, NULL, env);
  option_parser_add_option(op, optionmotif, env);

  /* -motifmis */
  optionmotifmis = option_new_uint_min_max("motifmis",
                             "specify maximum number of "
                             "mismatches in motif [0,3]",
                             &lo->motif.allowedmismatches,
                             (unsigned int)4,
                             (unsigned int)0,
                             (unsigned int)3,
                             env);
  option_parser_add_option(op, optionmotifmis, env);

  /* -vic */
  optionvic = option_new_uint_min_max("vic",
                        "specify the number of nucleotides (to the left and "
                        "to the right) that will be searched "
                        "for TSDs and/or motifs around 5' and 3' boundary "
                        "of predicted LTR retrotransposons",
                        &vicinityforcorrectboundaries,
                        60U,
                        1U,
                        500U,
                        env);
  option_parser_add_option(op, optionvic, env);

  /* -overlaps */
  lo->str_overlaps = str_new(env);
  optionoverlaps = option_new_choice("overlaps",
               "specify no|best|all",
               lo->str_overlaps,
               overlaps[0],
               overlaps,
               env);
  option_parser_add_option(op, optionoverlaps, env);

  /* -xdrop */
  optionxdrop = option_new_int_min("xdrop",
                        "specify xdropbelowscore for extension-alignment",
                        &lo->xdropbelowscore,
                        (int)5,
                        (int)0,
                        env);
  option_parser_add_option(op, optionxdrop, env);

  /* -mat */
  lo->arbitscores.gcd  = (int) 1;      /* set only for initialization,
                                        do not change! */
  optionmat = option_new_int_min("mat",
                        "specify matchscore for extension-alignment",
                        &lo->arbitscores.mat,
                        (int)2,
                        (int)1,
                        env);
  option_parser_add_option(op, optionmat, env);

  /* -mis */
  optionmis = option_new_int_max("mis",
                        "specify mismatchscore for extension-alignment",
                        &lo->arbitscores.mis,
                        (int)-2,
                        (int)-1,
                        env);
  option_parser_add_option(op, optionmis, env);

  /* -ins */
  optionins = option_new_int_max("ins",
                        "specify insertionscore for extension-alignment",
                        &lo->arbitscores.ins,
                        (int)-3,
                        (int)-1,
                        env);
  option_parser_add_option(op, optionins, env);

  /* -del */
  optiondel = option_new_int_max("del",
                        "specify deletionscore for extension-alignment",
                        &lo->arbitscores.del,
                        (int)-3,
                        (int)-1,
                        env);
  option_parser_add_option(op, optiondel, env);

  /* -v */
  optionv = option_new_bool("v",
                           "verbose mode",
                           &lo->verbosemode,
                           false,
                           env);
  option_parser_add_option(op, optionv, env);

  /* -longoutput */
  optionlongoutput = option_new_bool("longoutput",
                           "additional motif/TSD output",
                           &lo->longoutput,
                           false,
                           env);
  option_parser_add_option(op, optionlongoutput, env);

  /* -out */
  lo->fastaoutput = false;      /* by default no FASTA output */
  lo->str_fastaoutputfilename = str_new(env);
  optionout = option_new_string("out",
                             "specify FASTA outputfilename",
                             lo->str_fastaoutputfilename, NULL, env);
  option_parser_add_option(op, optionout, env);

  /* -outinner */
  lo->fastaoutputinnerregion = false;
  lo->str_fastaoutputfilenameinnerregion = str_new(env);
  optionoutinner = option_new_string("outinner",
                             "specify FASTA outputfilename for inner regions",
                             lo->str_fastaoutputfilenameinnerregion, NULL, env);
  option_parser_add_option(op, optionoutinner, env);

  /* -gff3 */
  lo->gff3output = false;       /* by default no gff3 output */
  lo->str_gff3filename = str_new(env);
  optiongff3 = option_new_string("gff3",
                             "specify GFF3 outputfilename",
                             lo->str_gff3filename, NULL, env);
  option_parser_add_option(op, optiongff3, env);

  /* implications */
  option_imply(optionmaxtsd, optionmintsd, env);
  option_imply(optionmotifmis, optionmotif, env);

  option_imply_either_2(optionlongoutput, optionmintsd, optionmotif, env);

  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  lo->vicinityforcorrectboundaries = (Seqpos) vicinityforcorrectboundaries;
  if (oprval == OPTIONPARSER_OK)
  {
    if (lo->repeatinfo.lmin > lo->repeatinfo.lmax)
    {
      env_error_set(env,"argument of -minlenltr is greater than argument of"
          " -maxlenltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if (lo->repeatinfo.dmin > lo->repeatinfo.dmax)
    {
      env_error_set(env,
          "argument of -mindistltr is greater than argument of -maxdistltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if (lo->minlengthTSD > lo->maxlengthTSD)
    {
      env_error_set(env,
          "argument of -mintsd is greater than argument of -maxtsd");
      oprval = OPTIONPARSER_ERROR;
    }

    /* If option motif is set,
       store characters, transform them later */
    if (option_is_set(optionmotif))
    {
      if (str_length(lo->motif.str_motif) != (unsigned long) 4)
      {
        env_error_set(env,
            "argument of -motif has not exactly 4 characters");
        oprval = OPTIONPARSER_ERROR;
      }
      lo->motif.firstleft = (Uchar)  str_get(lo->motif.str_motif)[0];
      lo->motif.secondleft = (Uchar)  str_get(lo->motif.str_motif)[1];
      lo->motif.firstright = (Uchar)  str_get(lo->motif.str_motif)[2];
      lo->motif.secondright = (Uchar)  str_get(lo->motif.str_motif)[3];
      /* default if motif specified */
      if (!option_is_set(optionmotifmis))
      {
        lo->motif.allowedmismatches = (unsigned int)0;
      }
    }

    /* If option overlaps is set */
    if (option_is_set(optionoverlaps))
    {
      if ( strcmp(str_get(lo->str_overlaps), "no") == 0)
      {
        lo->bestofoverlap = false;
        lo->nooverlapallowed = true;
      }
      else {
        if ( strcmp(str_get(lo->str_overlaps), "best") == 0 )
        {
          lo->bestofoverlap = true;
          lo->nooverlapallowed = false;
        }
        else {
          if ( strcmp(str_get(lo->str_overlaps), "all") == 0 )
          {
            lo->bestofoverlap = false;
            lo->nooverlapallowed = false;
          }
          else
          {
            env_error_set(env,
               "argument of -overlaps is not valid. Choose best|all|no");
            oprval = OPTIONPARSER_ERROR;
          }
        }
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
    if (option_is_set(optionout))
    {
      lo->fastaoutput = true;
    }

    /* if FASTA output inner region is set */
    if (option_is_set(optionoutinner))
    {
      lo->fastaoutputinnerregion = true;
    }

    /* if GFF3 output is set */
    if (option_is_set(optiongff3))
    {
      lo->gff3output = true;
    }
  }
  /* Error from optionparser occurred */
  else
  {
    /*oprval = OPTIONPARSER_ERROR;*/
  }

  option_parser_delete(op, env);
  return oprval;
}

void wrapltrharvestoptions(LTRharvestoptions *lo,Env *env)
{
  /* no checking if error occurs, since errors have been output before */
  str_delete(lo->str_indexname,env);
  str_delete(lo->str_fastaoutputfilename,env);
  str_delete(lo->str_fastaoutputfilenameinnerregion,env);
  str_delete(lo->str_gff3filename,env);
  str_delete(lo->str_overlaps,env);
  str_delete(lo->motif.str_motif,env);
}

int ltrharvestoptions(LTRharvestoptions *lo, int argc, const char **argv,
                         Env *env)
{
  int parsed_args;
  OPrval rval;

  env_error_check(env);

  /** init LTRharvestoptions lo **/
  rval = parse_options(&parsed_args, lo, argc, argv, env);
  if (rval == OPTIONPARSER_OK)
  {
    if (parsed_args != argc)
    {
      env_error_set(env, "Listing of options and arguments is not correct.");
      rval = OPTIONPARSER_ERROR;
    }
  }
  return (rval == OPTIONPARSER_OK) ? 0: - 1;
}
