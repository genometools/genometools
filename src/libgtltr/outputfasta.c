/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/alphadef.h"
//#include "libgtmatch/alphabet.pr"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-seqread.pr"
#include "libgtmatch/pos2seqnum.pr"
#include "libgtmatch/echoseq.pr"

//#include "fhandledef.h"
//#include "inputsymbol.h"

//#include "multiseq.pr"
//#include "multiseq-adv.pr"
//#include "filehandle.pr"

#include "ltrharvest-opt.h"
#include "repeattypes.h"

/*
 The datetype Fastaoutinfo aggregates the values needed to show a fasta file.
 */

#define FASTASEPARATOR '>'

typedef struct
{
  Sequentialsuffixarrayreader *ssar;
  const Encodedsequence *encseq; // encoded sequence
  const Alphabet *alpha;     // the alphabet
  const Uchar *characters;   // for visible characters
  const char *destab;        // pointer on descriptions
  unsigned long *descendtab; // positions of desc-separators
  Seqpos totallength;  // totallength of encseq
  unsigned long numofdbsequences; // num of sequences in suffix array
  unsigned int linewidth; // the line width to show the alignment
  Seqpos *markpos; // positions of SEPARATOR symbols
  bool showseqnum;     // with or without the sequence number
  FILE *formatout;     // file pointer to show the alignment
} Fastaoutinfo;


/*
 The following function processes one predicted LTR element, i.e.
 a FASTA entry is created from the predicted LTR element.
 */
static int showpredictionfastasequence(Fastaoutinfo *info, Seqpos startpos, 
                    Seqpos len, /*@unused@*/Str *str_indexfilename, Env *env)
{
  Seqpos i,
         k,
	 offset;
  unsigned long seqnum = 
                  getrecordnumSeqpos(info->markpos, info->numofdbsequences, 
		                     info->totallength, startpos, env);
  unsigned long desclen;
  const char *desptr;

  (void) putc(FASTASEPARATOR,info->formatout);
  if(info->showseqnum)
  {
    fprintf(info->formatout,"%lu ", seqnum);
  }
  if(seqnum == 0)
  {
    offset = 0;
  }
  else
  {
    offset = info->markpos[seqnum-1]+1;
  }
  fprintf(info->formatout,"(seq-nr) [" FormatSeqpos ","
                                               FormatSeqpos "] ",
		       // increase by one for output 
                       PRINTSeqposcast(startpos - offset + 1),   
		       // increase by one for output 
                       PRINTSeqposcast(startpos - offset + len) );
  //fprintf(fastaoutinfo->formatout,"from %s ", indexfilename);

  // if there are sequence descriptions
  //check if descentab in demand !!!!
  desptr = retriesequencedescription(&desclen,
                                     info->destab,
                                     info->descendtab,
                                     seqnum);
  for(i=0; i < desclen; i++) 
  {
    fprintf(info->formatout, "%c", desptr[i]);
  }

  /*if(fastaoutinfo->multiseq->descspace.spaceUchar != NULL)
  {
    desclength = DESCRIPTIONLENGTH(fastaoutinfo->multiseq,seqnum);
    if(WRITETOFILEHANDLE(DESCRIPTIONPTR(fastaoutinfo->multiseq,seqnum),
	  (Uint) sizeof(Uchar),
	  desclength,
	  fastaoutinfo->formatout) != 0)
    {
      return (Sint) -1;
    }
  }*/

  (void) putc('\n',info->formatout);
  for(k=0, i = startpos; i < startpos + len; i++, k++)
  {
    if(k >= info->linewidth)
    {
      fprintf(info->formatout, "\n");
      k=0;
    }
    fprintf(info->formatout, "%c", 
      info->characters[getencodedchar(info->encseq, i, Forwardmode)] );
  }
  (void) putc('\n',info->formatout);

  return 0;
}

/*
 The following function processes all predicted LTR elements with
 the function from the apply pointer.
 */
static int overallpredictionsequences(const LTRharvestoptions *lo,
    Sequentialsuffixarrayreader *ssar,
    bool innerregion,
    void *applyinfo,
    int(*apply)(Fastaoutinfo *,Seqpos, Seqpos, Str*, Env *env),
    Env *env)
{
  unsigned int i;
  Seqpos start,
         end;
  LTRboundaries *boundaries;

  for(i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
  {
    boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
    if( !boundaries->skipped )
    {
      if(innerregion)
      {
        start = boundaries->leftLTR_3 + 1;
        end = boundaries->rightLTR_5 - 1;
        //DEBUG3(2,"overallpredictionsequences: pred=%lu, start=%lu, end=%lu\n",
	//  (Showuint) i,
	//  (Showuint) (start-seq),(Showuint) (end-seq));
      }
      else
      {
        start = boundaries->leftLTR_5;
        end = boundaries->rightLTR_3;
        //DEBUG3(2,"overallpredictionsequences: pred=%lu, start=%lu, end=%lu\n",
	//  (Showuint) i,
	//  (Showuint) (start-seq),(Showuint) (end-seq));
      }
      if(apply(applyinfo, 
	       innerregion ? boundaries->leftLTR_3 + 1: boundaries->leftLTR_5,
	       end - start + 1, lo->str_indexname, env) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

static int openoutfile(FILE **fp, char *filename, Env *env)
{
  *fp = fopen(filename,"w");
  if(*fp == NULL)
  {
    env_error_set(env, "cannot open file \"%s\"", filename);
    return -1;
  }
  return 0;
}

/*
 The following function prints sequences in multiple FASTA format to the 
 specified output.
 */
int showpredictionsmultiplefasta(const LTRharvestoptions *lo,
		       Seqpos *markpos,
		       bool innerregion,
		       unsigned int linewidth,
                       Sequentialsuffixarrayreader *ssar,
		       bool showseqnum,
		       Env *env)
{
  Fastaoutinfo fastaoutinfo;
  FILE *formatout = NULL;
  unsigned long *descendtab = NULL, 
                destablength; 
  const char *destab = NULL;

  if(openoutfile(&formatout, 
     innerregion ? str_get(lo->str_fastaoutputfilenameinnerregion) 
                           : str_get(lo->str_fastaoutputfilename), 
     env) != 0 )
  {
    return -1;
  }

  fastaoutinfo.ssar = ssar;
  fastaoutinfo.encseq = encseqSequentialsuffixarrayreader(ssar);
  fastaoutinfo.alpha = alphabetSequentialsuffixarrayreader(ssar);
  fastaoutinfo.characters = getcharactersAlphabet(fastaoutinfo.alpha);
  fastaoutinfo.totallength = getencseqtotallength(
                               encseqSequentialsuffixarrayreader(ssar)); 
  fastaoutinfo.numofdbsequences = 
                      numofdbsequencesSequentialsuffixarrayreader(ssar);
  fastaoutinfo.linewidth = linewidth;
  fastaoutinfo.showseqnum = showseqnum;
  fastaoutinfo.formatout = formatout;
  fastaoutinfo.markpos = markpos;

  destablength = destablengthSequentialsuffixarrayreader(ssar); //unbedingt einkommentieren
  destab = destabSequentialsuffixarrayreader(ssar); //unbedingt einkommentieren 
  descendtab = calcdescendpositions(destab,
                                    destablength,
                                    fastaoutinfo.numofdbsequences,
                                    env);
  fastaoutinfo.destab = destab;
  fastaoutinfo.descendtab = descendtab;

  if(overallpredictionsequences(lo,
				ssar,
                                innerregion, 
                                &fastaoutinfo,
				showpredictionfastasequence,
				env) != 0)
  {
    return -1;
  }

  FREESPACE(descendtab);

  if(fclose(formatout) != 0)
  {
    env_error_set(env, "cannot close file \"%s\"", 
                       str_get(lo->str_fastaoutputfilename));
    return -1;
  }
  
  return 0;
}
