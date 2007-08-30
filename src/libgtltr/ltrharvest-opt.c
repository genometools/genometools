/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
#include <stdio.h>

#include "libgtcore/env.h"
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/undef.h"

#include "libgtmatch/symboldef.h"
#include "libgtmatch/arraydef.h"
#include "libgtmatch/chardef.h"
#include "libgtmatch/alphadef.h"
#include "libgtmatch/alphabet.pr"

#include "ltrharvest-opt.h"
#include "repeattypes.h"

/*
 The following function shows all options that are set by default or from
  the user on stdout.
*/

void showuserdefinedoptionsandvalues(LTRharvestoptions *lo
/*
    bool verbosemode,
    bool bestofoverlap,
    bool nooverlapallowed,
    char *indexname,
    Arbitraryscores *arbitscores,
    unsigned int minseedlength,
    Xdropscore xdropbelowscore,
    double similaritythreshold,
    RepeatInfo *repeatinfo,
    unsigned int minlengthTSD,
    unsigned int maxlengthTSD,
    Motif *motif,
    unsigned int vicinityforcorrectboundaries,
    Alphabet *alpha,
    bool fastaoutput,
    char *fastaoutputfilename,
    bool fastaoutputinnerregion,
    char *fastaoutputfilenameinnerregion,
    bool gff3output,
    char *gff3filename
*/
)
{
  printf("# user defined options and values:\n");
  if(lo->verbosemode)
  {
    printf("#   verbosemode: On\n");
  }
  else
  {
    printf("#   verbosemode: Off\n");
  }
  printf("#   indexname: %s\n", str_get(lo->str_indexname));
  if(lo->fastaoutput)
  {
    printf("#   outputfile: %s\n", str_get(lo->str_fastaoutputfilename));
  }
  if(lo->fastaoutputinnerregion)
  {
    printf("#   outputfile inner region: %s\n",
	str_get(lo->str_fastaoutputfilenameinnerregion));
  }
  if(lo->gff3output)
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
  if(lo->nooverlapallowed)
  {
    printf("#   overlaps: no\n");
  }
  else
  {
    if(lo->bestofoverlap)
    {
      printf("#   overlaps: best\n");
    }
    else
    {
      printf("#   overlaps: all\n");
    }
  }
  printf("#   minTSDlength: %lu\n",  lo->minlengthTSD);
  printf("#   maxTSDlength: %lu\n",  lo->maxlengthTSD);
  printf("#   leftboundarymotifLTR: %c%c\n", lo->motif.firstleft,
                                             lo->motif.secondleft);
  printf("#   rightboundarymotifLTR: %c%c\n", lo->motif.firstright,
                                              lo->motif.secondright);
  printf("#   motifmismatchesallowed: %u\n", lo->motif.allowedmismatches);
  printf("#   vicinity: %u nt\n", lo->vicinityforcorrectboundaries);
}

/*
   This function prints the arguments from argv on standard output.
 */
void printargsline(const char **argv, int argc)
{
  int i;

  printf("# args=");
  for(i=1; i<argc; i++)
  {
    printf("%s",argv[i]);
    if(i == (argc-1))
    {
      printf("\n");
    } else
    {
      printf(" ");
    }
  }
}

int testmotif(Motif *motif, Alphabet *alpha, Env *env)
{
  const Uchar *symbolmap;

  symbolmap = getsymbolmapAlphabet(alpha);
  if( UNDEFCHAR == symbolmap[(unsigned int)motif->firstleft])
  {
    env_error_set(env,"Illegal nucleotide character %c " 
                      "as argument to option -motif", motif->firstleft);
    //ERROR1("Illegal nucleotide character %c as argument to option -motif",
    //	   motif->firstleft);
    return -1;
  }
  if( UNDEFCHAR == symbolmap[(unsigned int)motif->secondleft] )
  {
    env_error_set(env,"Illegal nucleotide character %c " 
                      "as argument to option -motif", motif->secondleft);
    //ERROR1("Illegal nucleotide character %c as argument to option -motif",
    //	motif->secondleft);
    return -1;
  }
  if( UNDEFCHAR == symbolmap[(unsigned int)motif->firstright] )
  {
    env_error_set(env,"Illegal nucleotide character %c " 
                      "as argument to option -motif", motif->firstright);
    //ERROR1("Illegal nucleotide character %c as argument to option -motif",
    //	motif->firstright);
    return -1;
  }
  if( UNDEFCHAR == symbolmap[(unsigned int)motif->secondright] )
  {
    env_error_set(env,"Illegal nucleotide character %c " 
                      "as argument to option -motif", motif->secondright);
    //ERROR1("Illegal nucleotide character %c as argument to option -motif",
    //	motif->secondright);
    return -1;
  }

  return 0;
}

static void versionfunc(const char *progname)
{
  printf("%s version 1.0\n",progname);
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

  static const char *overlaps[] = {
    "best", /* the default */ 
    "no",
    "all"
  };

  env_error_check(env);
  op = option_parser_new("[option ...] -index filenameindex",
                         "Predict LTR retrotransposons.", env);
  option_parser_set_mailaddress(op,"<ellinghaus@stud.zbh.uni-hamburg.de>");

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
  optionmintsd = option_new_ulong_min_max("mintsd",
                              "specify minimum length for each TSD", 
                               &lo->minlengthTSD,
                               (unsigned long) 0,
                               (unsigned long) 0,
			       UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionmintsd, env);

  /* -maxtsd */
  optionmaxtsd = option_new_ulong_min_max("maxtsd",
                              "specify maximum length for each TSD", 
                               &lo->maxlengthTSD,
                               (unsigned long) 20,
                               (unsigned long) 0,
			       UNDEF_ULONG,
                               env);
  option_parser_add_option(op, optionmaxtsd, env);

  /* -motif */
  /* characters will be tranformed later 
     into characters from virtualtree alphabet */ 
  lo->motif.firstleft   = (unsigned char)'t'; 
  lo->motif.secondleft  = (unsigned char)'g';  
  lo->motif.firstright  = (unsigned char)'c';
  lo->motif.secondright = (unsigned char)'a';
  lo->motif.str_motif = str_new(env);
  optionmotif = option_new_string("motif",
                             "specify 2 nucleotides startmotif + "
                             "2 nucleotides endmotif : ****",
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
                        &lo->vicinityforcorrectboundaries,
		        (unsigned int)60,
		        (unsigned int)1,
		        (unsigned int)500,
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
  //option_imply(optionminlenltr, optionmaxlenltr, env);
  //option_imply(optionmaxlenltr, optionminlenltr, env);
  //option_imply(optionmindistltr, optionmaxdistltr, env);
  //option_imply(optionmaxdistltr, optionmindistltr, env);
  //option_imply(optionmintsd, optionmaxtsd, env);
  //option_imply(optionmaxtsd, optionmintsd, env);
  option_imply(optionmotifmis, optionmotif, env);
  
  option_imply_either_2(optionlongoutput, optionmintsd, optionmotif, env);

  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  //printf("oprval = %d\n", oprval);
  if (oprval == OPTIONPARSER_OK)
  {
    if (lo->repeatinfo.lmin > lo->repeatinfo.lmax)
    {
      env_error_set(env,"argument of -minlenltr is greater than argument of"
	  " -maxlenltr");
      oprval = OPTIONPARSER_ERROR;
    } 
    if(lo->repeatinfo.dmin > lo->repeatinfo.dmax) 
    {
      env_error_set(env,
	  "argument of -mindistltr is greater than argument of -maxdistltr");
      oprval = OPTIONPARSER_ERROR;
    }
    if(lo->minlengthTSD > lo->maxlengthTSD)
    {
      env_error_set(env,
	  "argument of -mintsd is greater than argument of -maxtsd");
      oprval = OPTIONPARSER_ERROR;
    }

    /* If option motif is set,
       store characters, transform them later,
       siehe ehemalige main.c Z.301 */
    if (option_is_set(optionmotif))
    {
      if (str_length(lo->motif.str_motif) > (unsigned long) 4)
      {
	env_error_set(env,
	    "argument of -motif has more than 4 characters");
	oprval = OPTIONPARSER_ERROR;
      }
      lo->motif.firstleft = 
	(unsigned char) str_get(lo->motif.str_motif)[0];
      lo->motif.secondleft = 
	(unsigned char) str_get(lo->motif.str_motif)[1];
      lo->motif.firstright = 
	(unsigned char) str_get(lo->motif.str_motif)[2];
      lo->motif.secondright = 
	(unsigned char) str_get(lo->motif.str_motif)[3];
      /* default if motif specified */
      if (!option_is_set(optionmotifmis))
      {
	lo->motif.allowedmismatches = (unsigned int)0;
      }
    }

    /* If option overlaps is set */
    if (option_is_set(optionoverlaps))
    {
      if( strcmp(str_get(lo->str_overlaps), "no") == 0)
      {
	lo->bestofoverlap = false; 
	lo->nooverlapallowed = true;
      }
      else { 
	if( strcmp(str_get(lo->str_overlaps), "best") == 0 )
	{
	  lo->bestofoverlap = true; 
	  lo->nooverlapallowed = false;
	}
	else { 
	  if( strcmp(str_get(lo->str_overlaps), "all") == 0 )
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
    oprval = OPTIONPARSER_ERROR;
  }

  option_parser_delete(op, env);
  return oprval;
}

/*
 The following function shows all options that are set by default or from 
 the user on stdout.
 */
/*
static void showuserdefinedoptionsandvalues(LTRharvestoptions *lo)
{
  printf("# user defined options and values:\n");
  if(lo->verbosemode)
  {
    printf("#   verbosemode: On\n");
  }
  else
  {
    printf("#   verbosemode: Off\n");
  }
  printf("#   indexname: %s\n", str_get(lo->str_indexname));
  if(lo->fastaoutput)
  {
    printf("#   outputfile: %s\n", str_get(lo->str_fastaoutputfilename));
  }
  if(lo->fastaoutputinnerregion)
  {
    printf("#   outputfile inner region: %s\n", 
           str_get(lo->str_fastaoutputfilenameinnerregion));
  }
  if(lo->gff3output)
  {
    printf("#   outputfile gff3 format: %s\n", str_get(lo->str_gff3filename));
  }
  printf("#   xdropbelowscore: %u\n", lo->xdropbelowscore);
  printf("#   similaritythreshold: %.2f\n", 
      lo->similaritythreshold);
  printf("#   minseedlength: %lu\n", lo->minseedlength);
  printf("#   matchscore: %d\n",  lo->arbitscores.mat);
  printf("#   mismatchscore: %d\n", lo->arbitscores.mis);
  printf("#   insertionscore: %d\n", lo->arbitscores.ins);
  printf("#   deletionscore: %d\n", lo->arbitscores.del);
  printf("#   minLTRlength: %lu\n", lo->repeatinfo.lmin);
  printf("#   maxLTRlength: %lu\n", lo->repeatinfo.lmax);
  printf("#   minLTRdistance: %lu\n", lo->repeatinfo.dmin);
  printf("#   maxLTRdistance: %lu\n", lo->repeatinfo.dmax);
  if(lo->nooverlapallowed)
  {
    printf("#   overlaps: no\n");
  }
  else
  {
    if(lo->bestofoverlap)
    {
      printf("#   overlaps: best\n");
    }
    else
    {
      printf("#   overlaps: all\n");
    }
  } 
  printf("#   minTSDlength: %lu\n", lo->minlengthTSD);
  printf("#   maxTSDlength: %lu\n", lo->maxlengthTSD);
  //printf("#   leftboundarymotifLTR: %c%c\n", 
  //               alpha->characters[motif->firstleft],
  //		 alpha->characters[motif->secondleft]);
  //printf("#   rightboundarymotifLTR: %c%c\n", 
  //               alpha->characters[motif->firstright],
  //		 alpha->characters[motif->secondright]);
  printf("#   motifmismatchesallowed: %u\n", lo->motif.allowedmismatches);
  printf("#   vicinity: %u nt\n",  lo->vicinityforcorrectboundaries);
}
*/

int wrapltrharvestoptions(LTRharvestoptions *lo,Env *env)
{
  /* no checking if error occurs, since errors have been output before */
  lo->env = NULL;
  str_delete(lo->str_indexname,env);
  str_delete(lo->str_fastaoutputfilename,env);
  str_delete(lo->str_fastaoutputfilenameinnerregion,env);
  str_delete(lo->str_gff3filename,env);
  str_delete(lo->str_overlaps,env);
  str_delete(lo->motif.str_motif,env);
  return 0;
}

int ltrharvestoptions(LTRharvestoptions *lo, int argc, const char **argv,
                        Env *env)
{
  int parsed_args, retval = 0;
  OPrval rval;

  env_error_check(env);
 
  /** init LTRharvestoptions lo **/
  //Virtualtree virtualtree; // ersetzen
  //lo->demand = TISTAB | SUFTAB | LCPTAB | BWTTAB | DESTAB; //ersetzen
  
  rval = parse_options(&parsed_args, lo, argc, argv, env);
  if (rval == OPTIONPARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == OPTIONPARSER_REQUESTS_EXIT)
    {
      retval = 2;
    } else
    {
      if (rval == OPTIONPARSER_OK)
      {
	assert(parsed_args == argc);
	/* transform motif from argv */
	/*if( transformmotif(lo) != 0)
	{
	  STANDARDMESSAGE;
	}
	printargsline(argv,argc); 
        */
	/* if verbosemode, show chosen options and values */
	//showuserdefinedoptionsandvalues(lo); 
      }
    }
  }
  return retval;
}
