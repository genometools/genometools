/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/readmode-def.h"
#include "libgtmatch/pos2seqnum.pr"
#include "libgtmatch/symboldef.h"
#include "libgtmatch/alphabet.pr"

#include "ltrharvest-opt.h"

/*
   The following function prints the predicted LTR retrotransposon
   results to stdout.
 */
/* CAUTION: For output the positions will be decreased by one,
 *          since normally in annotation first base in sequence
 *          is position 1 instead of 0.
 */
int showinfoiffoundfullLTRs(LTRharvestoptions *lo,
    Suffixarray *suffixarray,
    Env *env)
/*
    ArrayLTRboundaries *arrayLTRboundaries,
    Virtualtree *virtualtree,
    BOOL longoutput,
    unsigned int minlengthTSD,
    Motif *motif
*/
{
  LTRboundaries *boundaries;
  unsigned long h, 
                i;
  unsigned int contignumber;
  Seqpos offset;
  Seqpos *markpos = NULL;
  unsigned long numofdbsequences = suffixarray->numofdbsequences;
  const Uchar *characters;
 
  /* in order to get to visible dna characters */
  characters = getcharactersAlphabet(suffixarray->alpha);

  /* count unique LTR predictions */
  //unsigned int counter = 0;
  /*for(i = 0; i < arrayLTRboundaries->nextfreeLTRboundaries; i++)
    {
    boundaries = &(arrayLTRboundaries->spaceLTRboundaries[i]);
    if( !boundaries->skipped )
    {
    counter++;
    }
    }*/

  if( numofdbsequences > (unsigned long)1 )
  {
    markpos = calculatemarkpositions(suffixarray->encseq, 
                                     numofdbsequences, 
				     env);
    if(markpos == NULL)
    {
      return -1;
    }
  }

  if(lo->longoutput)
  {
    if(lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
    {
      printf("No full LTR-pair predicted.\n");
    }
    else
    {
      printf("# predictions are reported in the following way\n");
      printf("# s(ret) e(ret) l(ret) ");
      printf("s(lLTR) e(lLTR) l(lLTR)");
      if(lo->minlengthTSD > (unsigned long) 1)
      {
	printf(" TSD l(TSD)");
      }
      if(lo->motif.allowedmismatches < (unsigned int)4)
      {
	printf(" m(lLTR)");
      }
      printf(" s(rLTR) e(rLTR) l(rLTR)");
      if(lo->minlengthTSD > (unsigned long) 1)
      {
	printf(" TSD l(TSD)");
      }
      if(lo->motif.allowedmismatches < (unsigned int)4)
      {
	printf(" m(rLTR)");
      }
      printf(" sim(LTRs)");
      printf(" seq-nr");
      printf("\n# where:\n");
      printf("# s = starting position\n");
      printf("# e = ending position\n");
      printf("# l = length\n");
      if(lo->motif.allowedmismatches < (unsigned int)4)
      {
	printf("# m = motif\n");
      }
      printf("# ret = LTR-retrotransposon\n");
      printf("# lLTR = left LTR\n");
      printf("# rLTR = right LTR\n");
      if(lo->minlengthTSD > (unsigned long) 1)
      {
	printf("# TSD = target site duplication\n");
      }
      printf("# sim = similarity\n");
      printf("# seq-nr = sequence number\n");

      /* print output sorted by contignumber*/
      for(h = 0; h < numofdbsequences; h++)
      {
	for(i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
	{
	  boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
	  contignumber = boundaries->contignumber;
	  if( (!boundaries->skipped) && (unsigned long)contignumber == h)
	  {
	    if( contignumber == 0)
	    {
	      offset = (Seqpos)0;
	    }
	    else
	    {
	      //offset =
	      //(virtualtree->multiseq.markpos.spaceunsigned int[contignumber-1])+1;
	      offset = markpos[contignumber-1]+(Seqpos)1;
	    }
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->leftLTR_5 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->rightLTR_3 -offset  + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast((boundaries->rightLTR_3 - boundaries->leftLTR_5 + 1)));
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->leftLTR_5 -offset  + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->leftLTR_3 -offset  + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast((boundaries->leftLTR_3 - boundaries->leftLTR_5 + (Seqpos)1)));
	    if(lo->minlengthTSD > (unsigned long) 1)
	    {
	      Seqpos j;

	      for(j = 0; j < boundaries->lenleftTSD; j++)
	      {
		printf("%c", 
		 characters[getencodedchar(suffixarray->encseq,
		    boundaries->leftLTR_5 - boundaries->lenleftTSD + j,
		    Forwardmode)]);
	      }
	      printf("  " FormatSeqpos "  ", 
	             PRINTSeqposcast(boundaries->lenleftTSD));
	    }
	    if(lo->motif.allowedmismatches < (unsigned int)4)
	    {
	      printf("%c%c..%c%c  ",
		  characters[getencodedchar(suffixarray->encseq,
                                 boundaries->leftLTR_5,
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
                                 boundaries->leftLTR_5+(Seqpos)1,
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
                                 boundaries->leftLTR_3-(Seqpos)1,
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
                                 boundaries->leftLTR_3,
				 Forwardmode)] );
		  /* 
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->leftLTR_5)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->leftLTR_5+1)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->leftLTR_3-1)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->leftLTR_3)] );
	          */
	    }
	    // increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->rightLTR_5 -offset + (Seqpos)1));
	    // increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast(boundaries->rightLTR_3 -offset + (Seqpos)1));
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast((boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1)));
	    if(lo->minlengthTSD > (unsigned long) 1)
	    {
	      Seqpos j;

	      for(j = 0; j < boundaries->lenrightTSD; j++)
	      {
		printf("%c", characters[getencodedchar(suffixarray->encseq, 
		    boundaries->rightLTR_3 + j + 1,
		    Forwardmode)]);
	      }
	      printf("  " FormatSeqpos "  ", 
	             PRINTSeqposcast(boundaries->lenrightTSD));
	    }
	    if(lo->motif.allowedmismatches < (unsigned int)4)
	    {
	      printf("%c%c..%c%c",
		  characters[getencodedchar(suffixarray->encseq,
                                 boundaries->rightLTR_5,
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
			         boundaries->rightLTR_5+(Seqpos)1,
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
			         boundaries->rightLTR_3-(Seqpos)1, 	  
				 Forwardmode)],
		  characters[getencodedchar(suffixarray->encseq,
			         boundaries->rightLTR_3,
				 Forwardmode)] );
		  /*
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->rightLTR_5)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->rightLTR_5+1)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->rightLTR_3-1)],
		  virtualtree->alpha.characters[
		  *(sequence + boundaries->rightLTR_3)] );
	          */
	    }
	    /* print similarity */
	    printf("  %.2f", boundaries->similarity);
	    /* print sequence number */
	    printf("  %u\n", contignumber);
	  }
	}
      }
    }
  }
  /* short output */
  else
  {
    if(lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
    {
      /* nothing */
    }
    else
    {
      /* print short output of full length LTR-retrotransposon(s) */
      printf("# predictions are reported in the following way\n");
      printf("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR)"
	  " s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr \n");
      printf("# where:\n");
      printf("# s = starting position\n");
      printf("# e = ending position\n");
      printf("# l = length\n");
      printf("# ret = LTR-retrotransposon\n");
      printf("# lLTR = left LTR\n");
      printf("# rLTR = right LTR\n");
      printf("# sim = similarity\n");
      printf("# seq-nr = sequence number\n");

      /* print output sorted by contignumber*/
      for(h = 0; h < numofdbsequences; h++)
      {
	for(i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
	{
	  boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
	  contignumber = boundaries->contignumber;
	  if( (!boundaries->skipped) && contignumber == h)
	  {
	    if( contignumber == 0)
	    {
	      offset = (Seqpos)0;
	    }
	    else
	    {
	      //offset = virtualtree->multiseq.markpos.spaceunsigned int[contignumber-1]+1;
	      offset = markpos[contignumber-1]+(Seqpos)1;
	    }

	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->leftLTR_5 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->rightLTR_3 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( (boundaries->rightLTR_3 - boundaries->leftLTR_5 + (Seqpos)1)));
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->leftLTR_5 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->leftLTR_3 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( (boundaries->leftLTR_3 - boundaries->leftLTR_5 + (Seqpos)1)));
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->rightLTR_5 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( boundaries->rightLTR_3 -offset + (Seqpos)1));// increase by 1
	    printf(FormatSeqpos "  ",
		PRINTSeqposcast( (boundaries->rightLTR_3 - boundaries->rightLTR_5 + (Seqpos)1)));
	    /* print similarity */
	    printf("%.2f  ", boundaries->similarity);
	    /* print sequence number */
	    printf("%u\n", contignumber);
	  }
	}
      }
    }
  }

  FREESPACE(markpos);

  return 0;
}

//int printgff3format( // Funktion noch unvollstaendig
//  LTRharvestoptions *lo,
//  Suffixarray *suffixarray,
//  const Seqpos *markpos)
//{
//  LTRboundaries *boundaries;
//  unsigned long h, i;
//  Seqpos contiglen;
//  unsigned int contignumber,
//       idcounterRepregion = 0,
//       idcounterRetrotrans = 0,
//       idcounterLTR = 0,
//       idcounterTSD = 0,
//       idcounterMotif = 0,
//       offset,
//       desclength;
//  unsigned long numofdbsequences = suffixarray->numofdbsequences;
//  //Uint numofsequences = virtualtree->multiseq.numofsequences;
//  Seqpos totallength = getencseqtotallength(suffixarray->encseq);
//
//  FILE *fp = NULL;
//  fp = fopen(str_get(lo->str_gff3filename),"w");
//  if(fp == NULL)
//  {
//    //ERROR1("cannot open file \"%s\"", gff3filename);
//    return -1;
//  }
//
//  if(lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
//  {
//    /* no LTR-pairs predicted */
//    return 0;
//  }
//  else
//  {
//    fprintf(fp, "##gff-version 3\n");
//    /* print output sorted by contignumber*/
//    for(h = 0; h < numofdbsequences; h++)
//    {
//      /* contig is first sequence, and only one sequence in multiseq */
//      if( h == 0 && numofdbsequences == 1)
//      {
//	contiglen = totallength;
//      }
//      else
//      {
//	/* first sequence and more than one sequence in multiseq */
//	if( h == 0)
//	{
//	  contiglen = markpos[h];
//	}
//	else
//	{
//	  /* last sequence in multiseq */
//	  if(h == numofdbsequences - 1)  
//	  {
//	    contiglen = totallength - 1 - markpos[h-1];
//	  }
//	  else
//	  {
//	    contiglen = markpos[h] - markpos[h-1] - 1;
//	  }
//	}
//      }
//      fprintf(fp, "##sequence-region seq%lu 1 " FormatSeqpos "\n", 
//                  h, PRINTSeqposcast(contiglen));
//      fprintf(fp, "# ");
//      // if there are sequence descriptions    // HIER WEITER MACHEN
//      if(virtualtree->multiseq.descspace.spaceUchar != NULL)
//      {
//	desclength = DESCRIPTIONLENGTH(&virtualtree->multiseq,h);
//	if(WRITETOFILEHANDLE(DESCRIPTIONPTR(&virtualtree->multiseq,h),
//	      (Uint) sizeof(Uchar),
//	      desclength,
//	      fp) != 0)
//	{
//	  return (Sint) -1;
//	}
//      }
//
//      for(i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
//      { 
//	boundaries = &(lo->arrayLTRboundaries->spaceLTRboundaries[i]);
//	contignumber = boundaries->contignumber;
//	if( (!boundaries->skipped) && contignumber == h)
//	{
//	  if( contignumber == 0)
//	  {
//	    offset = 0;
//	  }
//	  else
//	  {
//	    offset = 
//	      (virtualtree->multiseq.markpos.spaceUint[contignumber-1])+1;
//	  }
//
//          /* repeat-region */
//	  fprintf(fp, "seq%lu\tLTRharvest\trepeat_region\t%lu\t%lu\t"
//	      ".\t?\t.\tID=RepeatReg%lu\n",
//	      (Showuint)contignumber,
//	      /* increase boundary position by one for output */
//	      (Showuint)boundaries->leftLTR_5 -offset + 1 
//		          - (Uint)boundaries->lenleftTSD,// increase by 1 
//	      (Showuint)boundaries->rightLTR_3 -offset + 1
//		          + (Uint)boundaries->lenrightTSD,// increase by 1 
//	      (Showuint)idcounterRepregion++ );
//
//          /* LTR retrotransposon */
//	  fprintf(fp, "seq%lu\tLTRharvest\tLTR_retrotransposon\t%lu\t%lu\t"
//	      ".\t?\t.\tID=LTRret%lu;Parent=RepeatReg%lu\n",
//	      (Showuint)contignumber,
//	      /* increase boundary position by one for output */
//	      (Showuint)boundaries->leftLTR_5 -offset + 1,// increase by 1 
//	      (Showuint)boundaries->rightLTR_3 -offset  + 1,// increase by 1 
//	      (Showuint)idcounterRetrotrans++,
//	      (Showuint)idcounterRepregion-1 );
//
//          /* LTRs */
//	  fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t%lu\t%lu\t"
//	      ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
//	      (Showuint)contignumber,
//	      /* increase boundary position by one for output */
//	      (Showuint)boundaries->leftLTR_5 -offset + 1,// increase by 1 
//	      (Showuint)boundaries->leftLTR_3 -offset + 1,// increase by 1 
//	      (Showuint)idcounterLTR++,
//	      (Showuint)idcounterRetrotrans-1 );
//	  fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t%lu\t%lu\t"
//	      ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
//	      (Showuint)contignumber,
//	      /* increase boundary position by one for output */
//	      (Showuint)boundaries->rightLTR_5 -offset + 1,// increase by 1 
//	      (Showuint)boundaries->rightLTR_3 -offset + 1,// increase by 1 
//	      (Showuint)idcounterLTR++,
//	      (Showuint)idcounterRetrotrans-1 );
//	      
//	  
//	  if(minlengthTSD > (Uint) 1)
//	  {
//            /* TSDs */
//	    fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->leftLTR_5 -offset + 1 
//		          - (Uint)boundaries->lenleftTSD,// increase by 1 
//		(Showuint)boundaries->leftLTR_5 -offset,// increase by 1 
//		(Showuint)idcounterTSD++,
//	        (Showuint)idcounterRepregion-1 );
//		
//	    fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->rightLTR_3 -offset + 2,// increase by 1 
//		(Showuint)boundaries->rightLTR_3 -offset + 1
//		          + (Uint)boundaries->lenrightTSD,// increase by 1 
//		(Showuint)idcounterTSD++,
//	        (Showuint)idcounterRepregion-1 );
//	  }
//	  
//	  if(motif->allowedmismatches < (Ushort)4)
//	  {
//	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->leftLTR_5 -offset + 1,// increase by 1 
//		(Showuint)boundaries->leftLTR_5 -offset + 2,// increase by 1 
//		(Showuint)idcounterMotif++,
//		(Showuint)idcounterRepregion-1 );
//	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->leftLTR_3 -offset,// increase by 1 
//		(Showuint)boundaries->leftLTR_3 -offset + 1,// increase by 1 
//		(Showuint)idcounterMotif++,
//		(Showuint)idcounterRepregion-1 );
//	    
//	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->rightLTR_5 -offset + 1,// increase by 1 
//		(Showuint)boundaries->rightLTR_5 -offset + 2,// increase by 1 
//		(Showuint)idcounterMotif++,
//		(Showuint)idcounterRepregion-1 );
//	    fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
//	                "%lu\t%lu\t"
//		        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
//		(Showuint)contignumber,
//		/* increase boundary position by one for output */
//		(Showuint)boundaries->rightLTR_3 -offset,// increase by 1 
//		(Showuint)boundaries->rightLTR_3 -offset + 1,// increase by 1 
//		(Showuint)idcounterMotif++,
//		(Showuint)idcounterRepregion-1 );
//	  }
//	}
//      }
//    }
//  }
//  if(fclose(fp) != 0)
//  {
//    ERROR0("cannot close file");
//    return (Sint) -1;
//  }
//  return (Sint) 0;
//}
