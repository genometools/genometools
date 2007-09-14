/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>

#include "libgtcore/env.h"

#include "libgtmatch/encseq-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/arraydef.h"
#include "libgtmatch/symboldef.h"
#include "libgtmatch/spacedef.h"
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-seqread.pr"
#include "libgtmatch/esa-mmsearch-def.h"
#include "libgtmatch/sfx-map.pr"
#include "libgtmatch/esa-mmsearch.pr"

#include "ltrharvest-opt.h"
#include "minmax.h"
#include "repeattypes.h"
#include "repeats.h"

//#include "libgtmatch/encseq-def.h"
//#include "libgtmatch/pos2seqnum.pr"

//#include "searchforLTRs.h"
//#include "repeats.h"
//#include "myxdrop.h"


/* Prof. Kurtz files */
/*
#include "types.h"
#include "minmax.h"
#include "debugdef.h"
#include "errordef.h"
#include "arraydef.h"
#include "virtualdef.h"
#include "vmatfind-def.h"
#include "xdropdef.h"
#include "multidef.h"

#include "alphabet.pr"
#include "multiseq-adv.pr"
#include "multiseq.pr"
#include "readvirt.pr"
#include "detpfxlen.pr"

// my files
#include "searchforLTRs.h"
#include "repeattypes.h"
#include "repeats.pr"
*/

/*
 The following function searches for TSDs and/or a specified palindromic 
 motif at the 5'-border of left LTR and 3'-border of right LTR. Thereby, 
 all maximal repeats from the vicinity are processed one after another
 to find the TSD with the minimum deviation with regard to the boundary
 position from the x-drop alignment. If also a motif is searched,
 a simple motif check at the boundaries of the TSDs is performed.
 */

static void searchforbestTSDandormotifatborders(
    SubRepeatInfo *info,
    LTRharvestoptions *lo,
    Sequentialsuffixarrayreader *ssar,
    Seqpos *markpos,
    LTRboundaries *boundaries,
    unsigned int *motifmismatchesleftLTR,
    unsigned int *motifmismatchesrightLTR,
    Env *env
    )
{
  unsigned long i;
  Seqpos offset,
         motifpos1,
         motifpos2,
	 dbseqpos1,
	 dbseqpos2;
  Seqpos back, forward;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR;
  Seqpos oldleftLTR_5  = boundaries->leftLTR_5,
       oldrightLTR_3 = boundaries->rightLTR_3,
       difffromoldboundary1 = 0,
       difffromoldboundary2 = 0;
  unsigned int hitcounter = 0;
  const Encodedsequence *encseq = encseqSequentialsuffixarrayreader(ssar);

  env_error_check(env);

  if( boundaries->contignumber == 0)
  {
    offset = 0;
  }
  else
  {
    offset = markpos[boundaries->contignumber-1]+1;
  }

  if( info->repeats.nextfreeRepeat > 0 )
  {
    boundaries->tsd = true;
  }
  boundaries->motif_near_tsd = false;
  
  /*printf("  old boundary pos leftLTR_5  = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->leftLTR_5 - offset));
  printf("  old boundary pos rightLTR_3 = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->rightLTR_3 - offset));
  */
 
  for(i = 0; i < info->repeats.nextfreeRepeat; i++)
  {
    // dbseqpos1 is the first position after the left repeat
    motifpos1 = info->repeats.spaceRepeat[i].pos1 + 
                           info->repeats.spaceRepeat[i].len;
    // dbseqpos2 is two positions before the right repeat
    motifpos2 = info->repeats.spaceRepeat[i].pos1
		  + info->repeats.spaceRepeat[i].offset - 2;
    dbseqpos1 = motifpos1;
    dbseqpos2 = motifpos2;

    for(back = 0; 
        back < info->repeats.spaceRepeat[i].len - info->lmin + 1; 
	back++)
    {
      for(forward = 0; 
	  forward < info->repeats.spaceRepeat[i].len - 
	            info->lmin + 1 - back;
	  forward++)
      {
	tmp_motifmismatchesleftLTR = tmp_motifmismatchesrightLTR = 0;
	//if( *(dbseqpos1 - back) != motif->firstleft) 
	if( getencodedchar(encseq, dbseqpos1 - back, Forwardmode) 
	    != lo->motif.firstleft) 
	{
	  tmp_motifmismatchesleftLTR++;
	}
	if( getencodedchar(encseq, dbseqpos1 + 1 - back, Forwardmode) 
	    != lo->motif.secondleft)
	{
	  tmp_motifmismatchesleftLTR++;
	}
        if( getencodedchar(encseq,   dbseqpos2 + forward, Forwardmode) 
	    != lo->motif.firstright)
	{
	  tmp_motifmismatchesrightLTR++;
	}
        if( getencodedchar(encseq, dbseqpos2 + 1 + forward, Forwardmode) 
	    != lo->motif.secondright)
	{
	  tmp_motifmismatchesrightLTR++;
	}
	
	if(tmp_motifmismatchesleftLTR <= lo->motif.allowedmismatches
	   &&
	   tmp_motifmismatchesrightLTR <= lo->motif.allowedmismatches
	  )
	{
	   Seqpos tsd_len;
	   tsd_len = info->repeats.spaceRepeat[i].len - back - forward;
	 
	   // TSD length too big 
	   if(tsd_len > (Seqpos)info->lmax)
	   {
             // nothing
	   }
	   // first hit
	   else if( !boundaries->motif_near_tsd )
	   {
	     unsigned int max, min;
	   
             // save number of mismatches
             *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
             *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

	     // adjust boundaries
	     boundaries->motif_near_tsd = true;
	     boundaries->leftLTR_5  = motifpos1 - back;
	     boundaries->rightLTR_3 = motifpos2 + 1 + forward;

             // store TSD length
             boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;
	     /*printf("first hit\n");
	     printf("boundaries->leftLTR_5 = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->leftLTR_5 - offset) );
	     printf("boundaries->rightLTR_3 = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->rightLTR_3 - offset) );
             printf("length left TSD: " FormatSeqpos "\n", 
	             PRINTSeqposcast(boundaries->lenleftTSD));
             printf("length right TSD: " FormatSeqpos "\n",
	             PRINTSeqposcast(boundaries->lenrightTSD));*/
	     
	     max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
	     min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
             difffromoldboundary1 = max - min; 

             max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
	     min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
             difffromoldboundary2 = max - min; 

	     hitcounter++;
	   }
	   else
	   {
	     unsigned int max, 
		  min,
		  difffromnewboundary1,
		  difffromnewboundary2;
	      
	     // test if hit is nearer to old boundaries than previous hit
	     max = MAX(oldleftLTR_5, (motifpos1 - back));
	     min = MIN(oldleftLTR_5, (motifpos1 - back));
             difffromnewboundary1 = max - min;
	     max = MAX(oldrightLTR_3, (motifpos2 + 1 + forward));
             min = MIN(oldrightLTR_3, (motifpos2 + 1 + forward));
	     difffromnewboundary2 = max - min;

	     if( (difffromnewboundary1 + difffromnewboundary2) < 
		 (difffromoldboundary1 + difffromoldboundary2) )
	     {
                 
	         // save number of mismatches
                 *motifmismatchesleftLTR  = tmp_motifmismatchesleftLTR;
                 *motifmismatchesrightLTR = tmp_motifmismatchesrightLTR;

	         // adjust boundaries
                 boundaries->leftLTR_5  = motifpos1 - back;
	         boundaries->rightLTR_3 = motifpos2 + 1 + forward;
	     
		 // store TSD length
                 boundaries->lenleftTSD = boundaries->lenrightTSD = tsd_len;
		 /*printf("next better hit\n");
		 printf("boundaries->leftLTR_5 = " FormatSeqpos "\n",
		     PRINTSeqposcast(boundaries->leftLTR_5 - offset) );
		 printf("boundaries->rightLTR_3 = " FormatSeqpos "\n",
		     PRINTSeqposcast(boundaries->rightLTR_3 - offset) );
		 printf("length left TSD: " FormatSeqpos "\n", 
		     PRINTSeqposcast(boundaries->lenleftTSD));
		 printf("length right TSD: " FormatSeqpos "\n",
		     PRINTSeqposcast(boundaries->lenrightTSD));
                 //DEBUG1(1, "boundaries->leftLTR_5 = %lu\n",
		 //       (Showuint) boundaries->leftLTR_5 - multiseqoffset);
	         //DEBUG1(1, "boundaries->rightLTR_3 = %lu\n",
		 //       (Showuint) boundaries->rightLTR_3 - multiseqoffset);
                 //DEBUG1(1, "length left TSD: %lu\n", (Showuint)
		 //                                    boundaries->lenleftTSD);
                 //DEBUG1(1, "length right TSD: %lu\n",(Showuint)
		 //                                    boundaries->lenrightTSD);
                 printf("difffromnewboundary1 + difffromnewboundary2 ="
		           FormatSeqpos  " + " FormatSeqpos "<\n"
		 	   "difffromoldboundary1 + difffromoldboundary2 ="
		 	   FormatSeqpos  " + " FormatSeqpos "\n",
	                    PRINTSeqposcast(difffromnewboundary1),
	                    PRINTSeqposcast(difffromnewboundary2),
	                    PRINTSeqposcast(difffromoldboundary1),
	                    PRINTSeqposcast(difffromoldboundary2) );*/
                 
		 difffromoldboundary1 = difffromnewboundary1;
		 difffromoldboundary2 = difffromnewboundary2;
		 hitcounter++;
	     }
	   }
	}
      }
    }
  }

  /*
  if( boundaries->motif_near_tsd )
  {
    printf("found %hu time(s) TSD and motif:\n",
	              hitcounter);
    printf("  new boundary pos leftLTR_5  = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->leftLTR_5 - offset));
    printf("  new boundary pos rightLTR_3 = " FormatSeqpos "\n",
	         PRINTSeqposcast(boundaries->rightLTR_3 - offset));
  }
  else
  {
    printf("no TSD and motif found.\n");   
  }
  */
}

/*
 The following function searches only for a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static void searchformotifonlyborders(LTRharvestoptions *lo,
    LTRboundaries *boundaries,
    Sequentialsuffixarrayreader *ssar,
    Seqpos *markpos,
    Seqpos startleftLTR,
    Seqpos endleftLTR,
    Seqpos startrightLTR,
    Seqpos endrightLTR,
    unsigned int *motifmismatchesleftLTR,
    unsigned int *motifmismatchesrightLTR,
    Env *env
    )
{
  Seqpos offset = 0,
         i;
  bool motif1 = false, 
       motif2 = false;
  unsigned int tmp_motifmismatchesleftLTR,
         tmp_motifmismatchesrightLTR,
	 motifmismatches_frombestmatch = 0; 
  const Encodedsequence *encseq = encseqSequentialsuffixarrayreader(ssar);
  Seqpos oldleftLTR_5  = boundaries->leftLTR_5,
         oldrightLTR_3 = boundaries->rightLTR_3,
         difffromoldboundary = 0;

  env_error_check(env);

  if( boundaries->contignumber == 0)
  {
    offset = 0;
  }
  else
  {
    offset = markpos[boundaries->contignumber-1]+(Seqpos)1;
  }
 
  ///// search for left motif around leftLTR_5 /////
  //DEBUG2(1, "search for left motif between pos = %lu up to "
  //	"pos end = %lu\n",
  //	(Showuint) startleftLTR - multiseqoffset,
  //	(Showuint) endleftLTR - multiseqoffset);
  
  for(i = startleftLTR; i < endleftLTR; i++)
  {
    tmp_motifmismatchesleftLTR = (unsigned int)0;
    if(getencodedchar(encseq, i, Forwardmode) != lo->motif.firstleft)
    {
      tmp_motifmismatchesleftLTR++; 
    }  
    if(getencodedchar(encseq, i+1, Forwardmode) != lo->motif.secondleft)
    {
      tmp_motifmismatchesleftLTR++; 
    }  
     
    if(tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR) 
                                <= lo->motif.allowedmismatches)
    {
       // first hit //
       if( !motif1 )
       {
         unsigned int max, min;

	 motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
	 boundaries->leftLTR_5 = i;
	 motif1 = true;
	 //DEBUG0(1, "found left motif:\n");
	 //DEBUG1(1, "  new boundary pos leftLTR_5 = %lu\n",
	 //          (Showuint) boundaries->leftLTR_5 - multiseqoffset);
	 max = MAX(oldleftLTR_5, boundaries->leftLTR_5);
	 min = MIN(oldleftLTR_5, boundaries->leftLTR_5);
	 difffromoldboundary = max - min; 
       }
       // next hit //
       else
       {
         Seqpos max, min, difffromnewboundary;
        
	 // test if hit is nearer to old boundaries than previous hit //
	 max = MAX(oldleftLTR_5, i);
	 min = MIN(oldleftLTR_5, i);
	 difffromnewboundary = max - min;

	 if( difffromnewboundary < difffromoldboundary )
	 {
	   motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
	   boundaries->leftLTR_5 = i;
           //DEBUG0(1, "next better hit\n");
	   //DEBUG1(1, "  new boundary pos leftLTR_5 = %lu\n",
	   //          (Showuint) boundaries->leftLTR_5 - multiseqoffset);
           difffromoldboundary = difffromnewboundary;
	 }
       }
    }
  }
#ifdef DEBUG
  if(i == endleftLTR && (!motif1) )
  {
    printf("no left motif found.\n");
  }
#endif
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  ///// search for right motif around rightLTR_3 /////
  //DEBUG2(1, "search right motif between pos = %lu"
  //          " up to pos end = %lu\n",
  //          (Showuint) startrightLTR - multiseqoffset,
  //	    (Showuint) endrightLTR - multiseqoffset);

  for(i = startrightLTR + (Seqpos)1; i <= endrightLTR; i++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if(getencodedchar(encseq, i, Forwardmode) != lo->motif.firstright)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if(getencodedchar(encseq, i-1, Forwardmode) != lo->motif.secondright)
    {
      tmp_motifmismatchesrightLTR++;
    }

    if(tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif.allowedmismatches)
    {
       // first hit //
       if( !motif2 )
       {
	 Seqpos max, min;

	 motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
	 boundaries->rightLTR_3 = i;
	 motif2 = true;
	 //DEBUG0(1, "found right motif:\n");
	 //DEBUG1(1, "  new boundary pos rightLTR_3 = %lu\n",
	 //          (Showuint) boundaries->rightLTR_3 - multiseqoffset);
	 max = MAX(oldrightLTR_3, boundaries->rightLTR_3);
	 min = MIN(oldrightLTR_3, boundaries->rightLTR_3);
	 difffromoldboundary = max - min; 
       }
       // next hit //
       else
       {
         Seqpos max, min, difffromnewboundary;

	 // test if hit is nearer to old boundaries than previous hit //
	 max = MAX(oldrightLTR_3, i);
	 min = MIN(oldrightLTR_3, i);
	 difffromnewboundary = max - min;
         if( difffromnewboundary < difffromoldboundary )
	 {
	   motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
	   boundaries->rightLTR_3 = i;
           //DEBUG0(1, "next better hit\n");
	   //DEBUG1(1, "  new boundary pos rightLTR_3 = %lu\n",
	   //          (Showuint) boundaries->rightLTR_3 - multiseqoffset);
           difffromoldboundary = difffromnewboundary;
	 }
       }
    }
  }
#ifdef DEBUG
  if(i > endrightLTR && (!motif2))    
  {
    printf("no right motif found.\n");
  }
#endif
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;
  
  if(motif1 && motif2) 
  {
    boundaries->motif_near_tsd = true; 
  }
  else
  {
    boundaries->motif_near_tsd = false;
  }

}

/*
 The following function searches for a specified palindromic motif at the 
 3'-border of left LTR and the 5'-border of right LTR.
 */

static void searchformotifonlyinside(LTRharvestoptions *lo,
    LTRboundaries *boundaries,
    Sequentialsuffixarrayreader *ssar,
    Seqpos *markpos,
    unsigned int *motifmismatchesleftLTR,
    unsigned int *motifmismatchesrightLTR,
    Env *env)
{
  bool motif1 = false, 
       motif2 = false;
  Seqpos startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         oldleftLTR_3  = boundaries->leftLTR_3,
         oldrightLTR_5 = boundaries->rightLTR_5,
         difffromoldboundary = 0;
  Seqpos offset = 0;
  unsigned int tmp_motifmismatchesleftLTR,
               tmp_motifmismatchesrightLTR,
               motifmismatches_frombestmatch = 0;
  
  Seqpos i;
  const Encodedsequence *encseq = encseqSequentialsuffixarrayreader(ssar);
  //Uchar *charptr;
  //Uchar *dbseq = virtualtree->multiseq.sequence;
 
  env_error_check(env);

  if( boundaries->contignumber == 0)
  {
    offset = 0;
  }
  else
  {
    offset = markpos[boundaries->contignumber-1]+(Seqpos)1;
  }

  //DEBUG2(1, "try vicinity left len = %lu, try vicinity right len = %lu"
  //    " (might be pruned if too large)\n",
  //    (Showuint) vicinity_left, (Showuint) vicinity_right);

  /** vicinity of 3'-border of left LTR **/
  /* do not align over 5'border of left LTR,
     in case of need decrease alignment length */
  if( (startleftLTR = boundaries->leftLTR_3 - 
         lo->vicinityforcorrectboundaries) <
      boundaries->leftLTR_5 + 2)
  {
    startleftLTR = boundaries->leftLTR_5 + 2;
  }
  /* do not align over 5'-border of right LTR */
  if( (endleftLTR = boundaries->leftLTR_3 + 
       lo->vicinityforcorrectboundaries) >
      boundaries->rightLTR_5 - 1)
  {
    endleftLTR = boundaries->rightLTR_5 - 1; 
  }
  /** vicinity of 5'-border of right LTR **/
  /* do not align over 3'-border of left LTR */
  if( (startrightLTR = boundaries->rightLTR_5 - 
         lo->vicinityforcorrectboundaries)
       < boundaries->leftLTR_3 + 1)
  {
    startrightLTR = boundaries->leftLTR_3 + 1;
  }
  /* do not align over 3'border of right LTR */
  if( (endrightLTR = boundaries->rightLTR_5 + 
       lo->vicinityforcorrectboundaries) >
      boundaries->rightLTR_3 - 2)
  {
    endrightLTR = boundaries->rightLTR_3 - 2; 
  }
  
  /**** search for right motif around leftLTR_3 ****/
  /*DEBUG2(1, "search for right motif between pos = %lu up to "
	"pos end = %lu\n",
	(Showuint) startleftLTR - multiseqoffset, 
	(Showuint) endleftLTR - multiseqoffset);*/
 
  for(i = startleftLTR + (Seqpos)1; i <= endleftLTR; i++)
  {
    tmp_motifmismatchesleftLTR = (unsigned int)0;
    if(getencodedchar(encseq, i, Forwardmode) != lo->motif.secondright)
    {
      tmp_motifmismatchesleftLTR++; 
    }
    if(getencodedchar(encseq, i-1, Forwardmode) != lo->motif.firstright)
    {
      tmp_motifmismatchesleftLTR++;
    } 
    
    if(tmp_motifmismatchesleftLTR + (*motifmismatchesleftLTR) 
                                <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if( !motif1 )
       {
         Seqpos max, min;

	 motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
	 boundaries->leftLTR_3 = i;
	 motif1 = true;
	 //DEBUG0(1, "found right motif:\n");
	 //DEBUG1(1, "  new boundary pos leftLTR_3 = %lu\n",
	 //          (Showuint) boundaries->leftLTR_3 - multiseqoffset);
	 max = MAX(oldleftLTR_3, boundaries->leftLTR_3);
	 min = MIN(oldleftLTR_3, boundaries->leftLTR_3);
	 difffromoldboundary = max - min; 
       }
       /* next hit */
       else
       {
         Seqpos max, min, difffromnewboundary;
        
	 /* test if hit is nearer to old boundaries than previous hit */
	 max = MAX(oldleftLTR_3, i);
	 min = MIN(oldleftLTR_3, i);
	 difffromnewboundary = max - min;

	 if( difffromnewboundary < difffromoldboundary )
	 {
	   motifmismatches_frombestmatch = tmp_motifmismatchesleftLTR;
	   boundaries->leftLTR_3 = i;
           //DEBUG0(1, "next better hit\n");
	   //DEBUG1(1, "  new boundary pos leftLTR_3 = %lu\n",
	   //          (Showuint) boundaries->leftLTR_3 - multiseqoffset);
           difffromoldboundary = difffromnewboundary;
	 }
       }
    }
  }
#ifdef DEBUG
  if(i > endleftLTR && (!motif1))
  {
    //DEBUG0(1, "no right motif found.\n");
  }
#endif
  *motifmismatchesleftLTR += motifmismatches_frombestmatch;
  motifmismatches_frombestmatch = 0;

  /**** search for left motif around rightLTR_5 ****/
  //DEBUG2(1, "search left motif between pos = %lu"
  //          " up to pos end = %lu\n",
  //          (Showuint) startrightLTR - multiseqoffset,
  //	    (Showuint) endrightLTR - multiseqoffset);

  for(i = startrightLTR ; i < endrightLTR; i++)
  {
    tmp_motifmismatchesrightLTR = 0;
    if(getencodedchar(encseq, i, Forwardmode) != lo->motif.firstleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if(getencodedchar(encseq, i+1, Forwardmode) != lo->motif.secondleft)
    {
      tmp_motifmismatchesrightLTR++;
    }
    if(tmp_motifmismatchesrightLTR + (*motifmismatchesrightLTR)
                                   <= lo->motif.allowedmismatches)
    {
       /* first hit */
       if( !motif2 )
       {
         unsigned int max, min;

	 motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
	 boundaries->rightLTR_5 = i;
	 motif2 = true;
	 //DEBUG0(1, "found left motif:\n");
	 //DEBUG1(1, "  new boundary pos rightLTR_5 = %lu\n",
	 //          (Showuint) boundaries->rightLTR_5 - multiseqoffset);
	 max = MAX(oldrightLTR_5, boundaries->rightLTR_5);
	 min = MIN(oldrightLTR_5, boundaries->rightLTR_5);
	 difffromoldboundary = max - min; 
       }
       /* next hit */
       else
       {
         unsigned int max, min, difffromnewboundary;
        
	 /* test if hit is nearer to old boundaries than previous hit */
	 max = MAX(oldrightLTR_5, i);
	 min = MIN(oldrightLTR_5, i);
	 difffromnewboundary = max - min;

	 if( difffromnewboundary < difffromoldboundary )
	 {
	   motifmismatches_frombestmatch = tmp_motifmismatchesrightLTR;
	   boundaries->rightLTR_5 = i;
           //DEBUG0(1, "next better hit\n");
	   //DEBUG1(1, "  new boundary pos rightLTR_5 = %lu\n",
	   //          (Showuint) boundaries->rightLTR_5 - multiseqoffset);
           difffromoldboundary = difffromnewboundary;
	 }
       }
    }
  }
#ifdef DEBUG
  if(i == endrightLTR && (!motif2))    
  {
    printf("no left motif found.\n");
    //DEBUG0(1, "no left motif found.\n");
  }
#endif
  *motifmismatchesrightLTR += motifmismatches_frombestmatch;
  
  if(motif1 && motif2) 
  {
    boundaries->motif_far_tsd = true; 
  }
  else
  {
    boundaries->motif_far_tsd = false;
  }
}

/*
 The following function searches for TSDs and/or a specified palindromic motif
 at the 5'-border of left LTR and 3'-border of right LTR.
 */

static int searchforTSDandorMotifoutside(
  LTRharvestoptions *lo,
  LTRboundaries *boundaries,
  Sequentialsuffixarrayreader *ssar,
  Seqpos *markpos,
  unsigned int *motifmismatchesleftLTR,
  unsigned int *motifmismatchesrightLTR,
  Env *env)
{
  Seqpos startleftLTR,
         endleftLTR,
         startrightLTR,
         endrightLTR,
         leftlen,
         rightlen;
  unsigned long contignumber = boundaries->contignumber;
  Seqpos offset;
  unsigned long numofdbsequences 
                  = numofdbsequencesSequentialsuffixarrayreader(ssar);
  Seqpos totallength = 
            getencseqtotallength(encseqSequentialsuffixarrayreader(ssar));
  SubRepeatInfo subrepeatinfo;
  const Encodedsequence *encseq = encseqSequentialsuffixarrayreader(ssar);
  
  ///*test
  const Uchar *characters;
  // in order to get to visible dna characters
  characters = getcharactersAlphabet(
                 alphabetSequentialsuffixarrayreader(ssar));
  //test*/

  //Multiseq *multiseqptr = &(virtualtree->multiseq);
  //Uchar *sequence  = multiseqptr->sequence;
  //Virtualtree subvirtualtree;
  //unsigned int demand = LCPTAB | TISTAB | SUFTAB | BWTTAB;

  //DEBUG2(1, "try vicinity left len = %lu, try vicinity right len = %lu"
  //    " (might be pruned if too large)\n",
  //    (Showuint) vicinity, (Showuint) vicinity);

  env_error_check(env);

  if( contignumber == 0)
  {
    offset = 0;
  }
  else
  {
    offset = markpos[contignumber-1]+1;
  }

  // check border cases //

  // vicinity of 5'-border of left LTR ///
  if(contignumber == 0)
  {
    // do not align over left sequence boundary,
    //   in case of need decrease alignment length //
    if( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries)
    {
      startleftLTR = 0;
    }
    else
    {
      startleftLTR = 
        boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries;
    }
  }
  else
  {
    // do not align over left separator symbol 
    //   at markpos.spaceunsigned int[contignumber-1],
    //   in case of need decrease alignment length //
    //DEBUG1(1, "boundaries->leftLTR_5 abs. position = %lu\n", 
    //       (Showuint) boundaries->leftLTR_5);
    //DEBUG1(1, "vicinity = %lu\n", (Showuint) vicinity);
    //DEBUG1(1, "multiseqptr->markpos.spaceunsigned int[contignumber-1]+1 = %ld\n", 
    //      (Showuint) multiseqptr->markpos.spaceunsigned int[contignumber-1]+1);
   
    if( boundaries->leftLTR_5 < lo->vicinityforcorrectboundaries )
    {
      startleftLTR = markpos[contignumber-1]+1;
    }
    else
    {
      if( ((startleftLTR = 
	      boundaries->leftLTR_5 - lo->vicinityforcorrectboundaries) <
	        markpos[contignumber-1]+1)
	    &&
	  (boundaries->leftLTR_5 >= markpos[contignumber-1]+1)
	)
      {
	startleftLTR = markpos[contignumber-1]+1;
      }
    }
  }
  // do not align over 3'-border of left LTR //
  if( (endleftLTR = boundaries->leftLTR_5 + lo->vicinityforcorrectboundaries)
       > boundaries->leftLTR_3 - 2 // -2 because of possible motif 
    )
  {
    endleftLTR = boundaries->leftLTR_3 - 2; 
  }
  leftlen = endleftLTR - startleftLTR + 1;

  /// vicinity of 3'-border of right LTR ///
  // do not align over 5'border of right LTR //
  if( (startrightLTR = 
         boundaries->rightLTR_3 - lo->vicinityforcorrectboundaries) <
      boundaries->rightLTR_5 + 2  // +2 because of possible motif //
    )
  {
    startrightLTR = boundaries->rightLTR_5 + 2;
  }
  if(contignumber == numofdbsequences - 1)
  {
    // do not align over right sequence boundary,
    //   in case of need decrease alignment length
    if( (endrightLTR = 
           boundaries->rightLTR_3 + lo->vicinityforcorrectboundaries) 
	> totallength - 1)
    {
      endrightLTR = totallength - 1;
    }
  }
  else
  {
    // do not align over right separator symbol 
    //   at markpos.spaceunsigned int[contignumber],
    //   in case of need decrease alignment length
    if( ((endrightLTR = 
            boundaries->rightLTR_3 + lo->vicinityforcorrectboundaries) >
	  markpos[contignumber]-1) 
        &&
        (boundaries->rightLTR_3 < markpos[contignumber])
      )
    {
      endrightLTR = markpos[contignumber]-1; 
    }
  }
  rightlen = endrightLTR - startrightLTR + 1; 


  // now, search for correct boundaries //

  ///// search for TSDs and/or motif ///// 
  if(lo->minlengthTSD > (unsigned long) 1) 
  {
    //DEBUG2(1, "searching for TSD and motif between pos = %lu up to "
    //	"pos end = %lu\n",
    //	(Showuint) startleftLTR - multiseqoffset, 
    //	(Showuint) endleftLTR - multiseqoffset);
    //DEBUG1(1, "use total vicinity len = %lu\n", 
    //	(Showuint) leftlen); 
    //DEBUG2(1, "searching for TSD and motif between pos = %lu"
    //	" up to pos end = %lu\n",
    //	(Showuint) startrightLTR - multiseqoffset, 
    //	(Showuint) endrightLTR - multiseqoffset);
    //DEBUG1(1, "use total vicinity len = %lu\n", 
    //	(Showuint) rightlen); 
    
    //makeemptyvirtualtree(&subvirtualtree);

    Uchar *dbseq, *query;
    Seqpos i;
    unsigned long k = 0;
    ALLOCASSIGNSPACE(dbseq,NULL,Uchar,leftlen);
    ALLOCASSIGNSPACE(query,NULL,Uchar,rightlen);

    for(i = startleftLTR; i <= endleftLTR; i++, k++) 
    {
      dbseq[k] = getencodedchar(encseq, i, Forwardmode);
    }
   
    for(k=0, i = startrightLTR; i <= endrightLTR; i++, k++)
    {
      query[k] = getencodedchar(encseq, i, Forwardmode); 
    }
 
    /*test
    printf("dbseq from =  " FormatSeqpos " to " FormatSeqpos "\n",
           PRINTSeqposcast(startleftLTR),
           PRINTSeqposcast(endleftLTR)
	   );
    for(i = 0; i<leftlen; i++)
    {
      printf("%c", characters[dbseq[i]]);
    }
    printf("\n");
    printf("query = " FormatSeqpos " to " FormatSeqpos "\n",
           PRINTSeqposcast(startrightLTR),
           PRINTSeqposcast(endrightLTR)
	   );
    for(i = 0; i<rightlen; i++)
    {
      printf("%c", characters[query[i]]);
    }
    printf("\n");
    test*/

    INITARRAY(&subrepeatinfo.repeats, Repeat);
    subrepeatinfo.lmin = lo->minlengthTSD;
    subrepeatinfo.lmax = lo->maxlengthTSD;
    assert(startleftLTR < startrightLTR);
    subrepeatinfo.offset1 = startleftLTR;
    subrepeatinfo.offset2 = startrightLTR;
    subrepeatinfo.envptr = env;

    if (sarrquerysubstringmatch(dbseq,
	  leftlen,
	  query,
	  (unsigned long)rightlen,
	  (unsigned int)lo->minlengthTSD,
	  alphabetSequentialsuffixarrayreader(ssar),
	  (void*)subsimpleexactselfmatchstore,//(void*)subshowrepeats,
	  &subrepeatinfo,
	  env) != 0)
    {
       return -1;
    }

    // 1. prepare multiseq //
    //initmultiseq(&subvirtualtree.multiseq);

    //addmultiseq(&subvirtualtree.multiseq, sequence + startleftLTR, leftlen);
    //addmultiseq(&subvirtualtree.multiseq, sequence + startrightLTR, rightlen);
    //
    //subvirtualtree.multiseq.numofqueryfiles = 0;
    //subvirtualtree.multiseq.numofquerysequences = 0;
    //subvirtualtree.multiseq.totalquerylength = 0;
    //
    //// 2. copy alphabet of vitualtree //
    //copyAlphabet(&subvirtualtree.alpha, &virtualtree->alpha);
    //
    //// 3. copy specialsymbols //
    //subvirtualtree.specialsymbols = virtualtree->specialsymbols; 

    // 4. determine prefixlength
    //   ACHTUNG: bei prefixlength != 0 bricht das Programm ab //
    //subvirtualtree.prefixlength = 0;

    ////
    //// funktioniert leider nicht
    ////subvirtualtree.prefixlength = 
    ////recommendedprefixlength(
    ////subvirtualtree.alpha.mapsize-1,
    ////subvirtualtree.multiseq.totallength,
    ////SIZEOFBCKENTRY
    ////);
    // //
    //
    //// fill subvirtualtree-tabs //
    //if( completevirtualtree(&subvirtualtree,
    //	  demand,
    // 	  NULL //showonstdout
    //	  ) != 0 )
    //{
    //  ERROR0("completevirtualtree failed.\n");
    //  return (Sint) -1;
    //}

    //DEBUG0(1, "showvirtualtreestatus (for TSD and Motif search):\n");
//#ifdef DEBUG
    //if (showvirtualtreestatus (&subvirtualtree, "(null)", showonstdout) != 0)
    //{
    //  ERROR0("showvirtualtreestatus failed.\n");
    //  return (Sint) -1;
    //}
//#endif

    //INITARRAY (&subrepeatinfo.repeats, Repeat);
    //subrepeatinfo.lmin = minlengthTSD;
    //subrepeatinfo.lmax = maxlengthTSD;
    //subrepeatinfo.separatorpos = leftlen;
    //subrepeatinfo.offset1 = startleftLTR;
    //subrepeatinfo.offset2 = startrightLTR;
    //DEBUG0(1, "the following are repeat candidates:\n");
    //
    //if( vmatmaxoutdynamic (&subvirtualtree, unsigned intConst (1),
    //      minlengthTSD,
    //      NULL, &subrepeatinfo,     // NULL,
    //      (void*)subsimpleexactselfmatchstore) != 0)
    //{
    //  ERROR0("vmatmaxoutdynamic failed.\n");
    //  return (Sint) -1;
    //} 
    //
    //freemultiseq(&subvirtualtree.multiseq);
    //if( freevirtualtree(&subvirtualtree) != 0 )
    //{
    //  ERROR0("freevirtualtree failed.\n");
    //  return (Sint) -1;
    //}

    FREESPACE(dbseq);
    FREESPACE(query);

    searchforbestTSDandormotifatborders(&subrepeatinfo, 
	lo,
	ssar,
	markpos,
	boundaries,
	motifmismatchesleftLTR,
	motifmismatchesrightLTR,
	env
	);

    FREEARRAY (&subrepeatinfo.repeats, Repeat);

  } else // no search for TSDs, search for motif only /////
  {
    
    /*DEBUG2(1, "searching for and motif only between pos = %lu up to "
	"pos end = %lu\n",
	(Showuint) startleftLTR - multiseqoffset, 
	(Showuint) endleftLTR - multiseqoffset);
    DEBUG1(1, "use total vicinity len = %lu\n", 
	(Showuint) leftlen); 

    DEBUG2(1, "searching for motif only between pos = %lu"
	" up to pos end = %lu\n",
	(Showuint) startrightLTR - multiseqoffset, 
	(Showuint) endrightLTR - multiseqoffset);
    DEBUG1(1, "use total vicinity len = %lu\n", 
	(Showuint) rightlen);*/
    
    searchformotifonlyborders(lo,
        boundaries, 
	ssar,
	markpos,
	startleftLTR,
	endleftLTR,
	startrightLTR,
	endrightLTR,
	motifmismatchesleftLTR,
	motifmismatchesrightLTR,
	env);
  }

  return 0;
}

/*
 The following function searches for TSD and/or a specified palindromic motif
 at the borders of left LTR and the right LTR, respectively.
 */
int findcorrectboundaries(
    LTRharvestoptions *lo,
    LTRboundaries *boundaries,
    Sequentialsuffixarrayreader *ssar,
    Seqpos *markpos,
    Env *env)
{
  unsigned int motifmismatchesleftLTR = 0,
               motifmismatchesrightLTR = 0;
 
  env_error_check(env);

  //DEBUG0(1, "searching for correct boundaries in vicinity...\n");

  /****** first: 5'-border of left LTR and 3'-border of right LTR *****/

  if( searchforTSDandorMotifoutside(lo,
                                    boundaries, 
				    ssar,
				    markpos,
				    &motifmismatchesleftLTR,
	                            &motifmismatchesrightLTR,
				    env) != 0 )
  {
    return -1;
  }

  /****** second: 3'-border of left LTR and 5'-border of right LTR *****/

  if( lo->motif.allowedmismatches < (unsigned int)4 )
  {
    //DEBUG0(1, "second: searching for motif only around 3'border of left"
    //	      " LTR and 5'-border of right LTR...\n");
    searchformotifonlyinside(lo,
        boundaries, 
	ssar,
	markpos,
	&motifmismatchesleftLTR,
	&motifmismatchesrightLTR,
	env);
  }

  return 0;
}
