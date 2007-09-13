/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/esa-seqread.h"

#include "ltrharvest-opt.h"
#include "repeattypes.h"

int printgff3format(
  LTRharvestoptions *lo,
  Sequentialsuffixarrayreader *ssar,
  const Seqpos *markpos,
  Env *env)
{
  LTRboundaries *boundaries;
  Seqpos contiglen,
         offset;
  unsigned long h, 
                i,
		contignumber,
       idcounterRepregion = 0,
       idcounterRetrotrans = 0,
       idcounterLTR = 0,
       idcounterTSD = 0,
       idcounterMotif = 0;
       //desclength;
 
 unsigned long numofdbsequences = 
                numofdbsequencesSequentialsuffixarrayreader(ssar);
  Seqpos totallength = getencseqtotallength(
                               encseqSequentialsuffixarrayreader(ssar)); 

  FILE *fp = NULL;
  fp = fopen(str_get(lo->str_gff3filename),"w");
  if(fp == NULL)
  {
    env_error_set(env, "cannot open file \"%s\"", str_get(lo->str_gff3filename));
    return -1;
  }

  if(lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
  {
    /* no LTR-pairs predicted */
    return 0;
  }
  else
  {
    fprintf(fp, "##gff-version 3\n");
    /* print output sorted by contignumber*/
    for(h = 0; h < numofdbsequences; h++)
    {
      /* contig is first sequence, and only one sequence in multiseq */
      if( h == 0 && numofdbsequences == 1)
      {
	contiglen = totallength;
      }
      else
      {
	/* first sequence and more than one sequence in suffixarray */
	if( h == 0)
	{
	  contiglen = markpos[h];
	}
	else
	{
	  /* last sequence in suffixarray */
	  if(h == numofdbsequences - 1)  
	  {
	    contiglen = totallength - 1 - markpos[h-1];
	  }
	  else
	  {
	    contiglen = markpos[h] - markpos[h-1] - 1;
	  }
	}
      }
      fprintf(fp, "##sequence-region seq%lu 1 " FormatSeqpos "\n", 
                  h, PRINTSeqposcast(contiglen));
      /*fprintf(fp, "# ");
      // if there are sequence descriptions  //HIER DESCRIPTIONS ERGAENZEN!!!
      if(virtualtree->multiseq.descspace.spaceUchar != NULL)
      {
	desclength = DESCRIPTIONLENGTH(&virtualtree->multiseq,h);
	if(WRITETOFILEHANDLE(DESCRIPTIONPTR(&virtualtree->multiseq,h),
	      (Uint) sizeof(Uchar),
	      desclength,
	      fp) != 0)
	{
	  return (Sint) -1;
	}
      }*/

      for(i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
      { 
	boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
	contignumber = boundaries->contignumber;
	if( (!boundaries->skipped) && contignumber == h)
	{
	  if( contignumber == 0)
	  {
	    offset = 0;
	  }
	  else
	  {
	    offset = markpos[contignumber-1]+1;
	  }

          /* repeat-region */
	  fprintf(fp, "seq%lu\tLTRharvest\trepeat_region\t" FormatSeqpos
	              "\t" FormatSeqpos "\t" ".\t?\t.\tID=RepeatReg%lu\n",
	      contignumber,
	      /* increase boundary position by one for output */
	      PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1 
		          - boundaries->lenleftTSD),// increase by 1 
	      PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1
		          + boundaries->lenrightTSD),// increase by 1 
	      idcounterRepregion++ );

          /* LTR retrotransposon */
	  fprintf(fp, "seq%lu\tLTRharvest\tLTR_retrotransposon\t" 
	              FormatSeqpos "\t" FormatSeqpos "\t"
	              ".\t?\t.\tID=LTRret%lu;Parent=RepeatReg%lu\n",
	      contignumber,
	      /* increase boundary position by one for output */
	      PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),// increase by 1 
	      PRINTSeqposcast(boundaries->rightLTR_3 -offset  + 1),// increase by 1 
	      idcounterRetrotrans++,
	      idcounterRepregion-1 );

          /* LTRs */
	  fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t" 
	              FormatSeqpos "\t" FormatSeqpos "\t"
	              ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
	      contignumber,
	      /* increase boundary position by one for output */
	      PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),// increase by 1 
	      PRINTSeqposcast(boundaries->leftLTR_3 -offset + 1),// increase by 1 
	      idcounterLTR++,
	      idcounterRetrotrans-1 );
	  fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t" 
	              FormatSeqpos "\t" FormatSeqpos "\t"
	              ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
	      contignumber,
	      /* increase boundary position by one for output */
	      PRINTSeqposcast(boundaries->rightLTR_5 -offset + 1),// increase by 1 
	      PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1),// increase by 1 
	      idcounterLTR++,
	      idcounterRetrotrans-1 );
	      
	  
	  if(lo->minlengthTSD > 1)
	  {
            /* TSDs */
	    fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1 
		          - boundaries->lenleftTSD),// increase by 1 
		PRINTSeqposcast(boundaries->leftLTR_5 -offset),// increase by 1 
		idcounterTSD++,
	        idcounterRepregion-1 );
		
	    fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->rightLTR_3 -offset + 2),// increase by 1 
		PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1
		          + boundaries->lenrightTSD),// increase by 1 
		idcounterTSD++,
	        idcounterRepregion-1 );
	  }
	  
	  if(lo->motif.allowedmismatches < 4)
	  {
	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),// increase by 1 
		PRINTSeqposcast(boundaries->leftLTR_5 -offset + 2),// increase by 1 
		idcounterMotif++,
		idcounterRepregion-1 );
	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->leftLTR_3 -offset),// increase by 1 
		PRINTSeqposcast(boundaries->leftLTR_3 -offset + 1),// increase by 1 
		idcounterMotif++,
		idcounterRepregion-1 );
	    
	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->rightLTR_5 -offset + 1),// increase by 1 
		PRINTSeqposcast(boundaries->rightLTR_5 -offset + 2),// increase by 1 
		idcounterMotif++,
		idcounterRepregion-1 );
	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
	                FormatSeqpos "\t" FormatSeqpos "\t"
		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
		contignumber,
		/* increase boundary position by one for output */
		PRINTSeqposcast(boundaries->rightLTR_3 -offset),// increase by 1 
		PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1),// increase by 1 
		idcounterMotif++,
		idcounterRepregion-1 );
	  }
	}
      }
    }
  }
  if(fclose(fp) != 0)
  {
    env_error_set(env, "cannot close file \"%s\"", str_get(lo->str_gff3filename));
    return -1;
  }
  return 0;
}
