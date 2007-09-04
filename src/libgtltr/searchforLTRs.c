/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtmatch/arraydef.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/pos2seqnum.pr"

#include "searchforLTRs.h"
#include "repeattypes.h"
#include "repeats.h"
#include "myxdrop.h"
#include "minmax.h"

/* 
 The following function checks, if the remaining candidate pairs still
 fulfill the length and distance constraints of LTRs.
 */

/*
static int checklengthanddistanceconstraints(
    LTRboundaries *boundaries,
    RepeatInfo *repeatinfo
    )
{
  unsigned int ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1,
       vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1,
       dist_between_LTRs = boundaries->rightLTR_5 - boundaries->leftLTR_5;

  if(ulen > repeatinfo->lmax || vlen > repeatinfo->lmax ||
     ulen < repeatinfo->lmin || vlen < repeatinfo->lmin ||
     dist_between_LTRs > repeatinfo->dmax || 
     dist_between_LTRs < repeatinfo->dmin ||
     boundaries->leftLTR_3 >= boundaries->rightLTR_5)
  {
    DEBUG0(1, "boundaries exceed length contraints or overlap.\n");
    boundaries->lengthdistconstraint = False;
    boundaries->similarity = 0.0;
    return (int) -1;
  }
  else
  {
    boundaries->lengthdistconstraint = True;
  }
  
  return 0;
}
*/

/*
 The following function determines the similarity of the two predicted LTRs.
 */

/*
static int determinesimilarityfromcandidate(
    LTRboundaries *boundaries,
    Virtualtree *virtualtree)
{
  unsigned int ulen = boundaries->leftLTR_3  - boundaries->leftLTR_5  + 1,
       vlen = boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1;
  Uchar *dbseq = virtualtree->multiseq.sequence;

  int edist;
  Greedyalignreservoir greedyalignreservoir;

#ifdef DEBUG
  unsigned int saflag;
#endif
  DEBUG0(1, "determine similarity threshold...\n");
  DEBUG1(1, "ulen = %lu\n", (Showuint) ulen);
  DEBUG1(1, "vlen = %lu\n", (Showuint) vlen);

  initgreedyalignreservoir(&greedyalignreservoir);
  edist = greedyedistalign(&greedyalignreservoir,
      False,
      0,
      dbseq + boundaries->leftLTR_5,
      (int) ulen,
      dbseq + boundaries->rightLTR_5,
      (int) vlen);
  if( edist < 0 )
  {
    return (int) -1;
  }

#ifdef DEBUG
  MAKESAFLAG(saflag,True,False,False);
  DEBUG1(2, "edist = %ld\n", (Showsint) edist);
  showalignment(saflag,
      stdout,
      UintConst(60),
      greedyalignreservoir.alignment.spaceEditoperation,
      greedyalignreservoir.alignment.nextfreeEditoperation,
      (unsigned int) edist,
      dbseq + boundaries->leftLTR_5,
      dbseq + boundaries->leftLTR_5,
      ulen,
      dbseq + boundaries->rightLTR_5, 
      dbseq + boundaries->rightLTR_5, 
      vlen,
      0,
      0);
#endif

  // determine similarity
  boundaries->similarity = 100.0 *
    (1 - (((double)(edist))/(MAX(ulen,vlen))));

  wraptgreedyalignreservoir(&greedyalignreservoir);

  return 0;
}
*/

void adjustboundariesfromXdropextension(
    Myxdropbest xdropbest_left,
    Myxdropbest xdropbest_right,
    Seqpos seed1_startpos,
    Seqpos seed2_startpos,
    Seqpos seed1_endpos,
    Seqpos seed2_endpos,
    LTRboundaries *boundaries,
    Seqpos offset
    )
{
  // left alignment
  boundaries->leftLTR_5  = seed1_startpos - xdropbest_left.ivalue;
  boundaries->rightLTR_5 = seed2_startpos - xdropbest_left.jvalue;
 
  // right alignment
  boundaries->leftLTR_3  = seed1_endpos + xdropbest_right.ivalue;
  boundaries->rightLTR_3 = seed2_endpos + xdropbest_right.jvalue;

  /*
  printf("Adjusted boundaries from Xdrop-alignment-extension:\n");
  printf("boundaries->leftLTR_5  = " FormatSeqpos "\n", 
            PRINTSeqposcast( boundaries->leftLTR_5 - offset));
  printf("boundaries->leftLTR_3  = " FormatSeqpos "\n", 
            PRINTSeqposcast( boundaries->leftLTR_3 - offset));
  printf("len leftLTR = " FormatSeqpos "\n", 
            PRINTSeqposcast( (boundaries->leftLTR_3 - boundaries->leftLTR_5 + 1)));
  printf("boundaries->rightLTR_5 = " FormatSeqpos "\n", 
            PRINTSeqposcast( boundaries->rightLTR_5 - offset));
  printf("boundaries->rightLTR_3 = " FormatSeqpos "\n", 
            PRINTSeqposcast( boundaries->rightLTR_3 - offset));
  printf("len rightLTR = " FormatSeqpos "\n", 
            PRINTSeqposcast( (boundaries->rightLTR_3 - boundaries->rightLTR_5 + 1)));
 */
}

/*
void shownewboundariesifadjusted(LTRboundaries *boundaries,
                                 unsigned int multiseqoffset)
{
  if(boundaries->tsd || 
    (boundaries->motif_near_tsd && boundaries->motif_far_tsd))
  {
    DEBUG0(1, "adjusted boundaries from TSD and Motif search:\n");
    DEBUG1(1, "boundaries->leftLTR_5  = %lu\n", 
	      (Showuint) boundaries->leftLTR_5 - multiseqoffset);
    DEBUG1(1, "boundaries->leftLTR_3  = %lu\n", 
	      (Showuint) boundaries->leftLTR_3 - multiseqoffset);
    DEBUG1(1, "boundaries->rightLTR_5 = %lu\n", 
	      (Showuint) boundaries->rightLTR_5 - multiseqoffset);
    DEBUG1(1, "boundaries->rightLTR_3 = %lu\n", 
	      (Showuint) boundaries->rightLTR_3 - multiseqoffset);
  }
  else
  {
    DEBUG0(1, "boundaries not adjusted.\n");
  }
}
*/

#define INITBOUNDARIES(B)\
	(B)->contignumber = repeatptr->contignumber;\
        (B)->leftLTR_5 = (Seqpos)0;\
        (B)->leftLTR_3 = (Seqpos)0;\
        (B)->rightLTR_5 = (Seqpos)0;\
        (B)->rightLTR_3 = (Seqpos)0;\
        (B)->lenleftTSD = (Seqpos)0;\
        (B)->lenrightTSD = (Seqpos)0;\
        (B)->tsd = false;\
        (B)->motif_near_tsd = false;\
        (B)->motif_far_tsd = false;\
        (B)->skipped = false;\
        (B)->similarity = 0.0;

/*
 The following function applies the filter algorithms one after another
 to all candidate pairs.
*/

int searchforLTRs (
  Suffixarray *suffixarray,
  LTRharvestoptions *lo,
  Env *env
  )
{
  
  unsigned int repeatcounter;
  ArrayMyfrontvalue fronts;
  Myxdropbest xdropbest_left;
  Myxdropbest xdropbest_right;
  Seqpos alilen = 0,
         totallength; 
  Repeat *repeatptr;
  LTRboundaries *boundaries;
  Seqpos offset = 0;
  unsigned long numofdbsequences = suffixarray->numofdbsequences;
  //Uchar *seq = virtualtree->multiseq.sequence;
  Seqpos *markpos = NULL; //fuer testen

  /*
  printf("xdropbelowscore = %d\n", lo->xdropbelowscore);
  printf("scores:\n");
  printf("mat  = %d\n", lo->arbitscores.mat);
  printf("mis  = %d\n", lo->arbitscores.mis);
  printf("ins  = %d\n", lo->arbitscores.ins);
  printf("del  = %d\n", lo->arbitscores.del);
  */

  for (repeatcounter = 0; repeatcounter < lo->repeatinfo.repeats.nextfreeRepeat;
       repeatcounter++)
  {
    //printf("\n\nAlignments of repeat nr. = %u :\n", repeatcounter+1);
   
    repeatptr = &(lo->repeatinfo.repeats.spaceRepeat[repeatcounter]);
    alilen = ((Seqpos)lo->repeatinfo.lmax) - repeatptr->len;
    
    /**** left (reverse) xdrop alignment ****/
    INITARRAY (&fronts, Myfrontvalue);
    if(alilen <= repeatptr->pos1)
    {
      evalxdroparbitscoresleft(&lo->arbitscores, 
	                       &xdropbest_left, 
			       &fronts, 
			       suffixarray->encseq,
			       suffixarray->encseq,
			       repeatptr->pos1, 
			       repeatptr->pos1 + repeatptr->offset,
			       //seq + repeatptr->pos1, 
			       //seq + repeatptr->pos1 + repeatptr->offset,
			       (int) alilen,
			       (int) alilen, 
			       (Xdropscore)lo->xdropbelowscore,
			       env);
    }
    else
    {
      evalxdroparbitscoresleft(&lo->arbitscores, 
	                       &xdropbest_left, 
			       &fronts, 
			       suffixarray->encseq,
			       suffixarray->encseq,
			       repeatptr->pos1, 
			       repeatptr->pos1 + repeatptr->offset,
			       //seq + repeatptr->pos1, 
			       //seq + repeatptr->pos1 + repeatptr->offset,
			       (int) repeatptr->pos1,
			       (int) (repeatptr->pos1 + repeatptr->offset), 
			       (Xdropscore)lo->xdropbelowscore,
			       env);
    }
    FREEARRAY (&fronts, Myfrontvalue);

    /**** right xdrop alignment ****/
    INITARRAY (&fronts, Myfrontvalue);
    totallength = getencseqtotallength(suffixarray->encseq);
    if(alilen <= totallength - (repeatptr->pos1 + repeatptr->offset +
		                repeatptr->len) )
    {
      evalxdroparbitscoresright (&lo->arbitscores, 
                                 &xdropbest_right, 
				 &fronts, 
				 //seq + repeatptr->pos1 + repeatptr->len,
				 //seq + repeatptr->pos1 + repeatptr->offset +
				 suffixarray->encseq,
				 suffixarray->encseq,
				 repeatptr->pos1 + repeatptr->len,
				 repeatptr->pos1 + repeatptr->offset +
				 repeatptr->len, 
				 (int) alilen, 
				 (int) alilen,
				 lo->xdropbelowscore,
				 env);
    } 
    else
    {
      evalxdroparbitscoresright(&lo->arbitscores, 
	                        &xdropbest_right, 
				&fronts, 
				//seq + repeatptr->pos1 + repeatptr->len,
				//seq + repeatptr->pos1 + repeatptr->offset +
				suffixarray->encseq,
				suffixarray->encseq,
				repeatptr->pos1 + repeatptr->len,
				repeatptr->pos1 + repeatptr->offset +
				repeatptr->len, 
				(int) (totallength - 
				(repeatptr->pos1 + repeatptr->len)), 
				(int) (totallength - 
				(repeatptr->pos1 + repeatptr->offset +
				 repeatptr->len)),
				lo->xdropbelowscore,
				env);
    }
    FREEARRAY (&fronts, Myfrontvalue);
 
    GETNEXTFREEINARRAY(boundaries,
	               &lo->arrayLTRboundaries,
		       LTRboundaries,
		       5);
    INITBOUNDARIES(boundaries);
// test 
    if( boundaries->contignumber == (unsigned long)0)
    {
      offset = 0;
    }
    else
    {
      markpos = calculatemarkpositions(suffixarray->encseq, 
                                     numofdbsequences, 
				     env);
      if(markpos == NULL)
      { 
        return -1;
      }
      offset = markpos[boundaries->contignumber-1];
      FREESPACE(markpos);
      //multiseqoffset = 
      //virtualtree->multiseq.markpos.spaceUint[boundaries->contignumber-1]+1;
    }
/* test
    printf("contig number: %lu\n",
               boundaries->contignumber);
    printf("offset to contig startpos: " FormatSeqpos "\n",
               PRINTSeqposcast(offset)); //brauche markpos array
    
    printf("Boundaries from vmatmaxoutdynamic:\n");
    
    printf("boundaries->leftLTR_5 abs. = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1));
    printf("boundaries->rightLTR_5 abs. = " FormatSeqpos "\n", 
	      PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset));
    
    printf("boundaries->leftLTR_5  = " FormatSeqpos "\n",
              PRINTSeqposcast(repeatptr->pos1 - offset));
    printf("boundaries->leftLTR_3  = " FormatSeqpos "\n", 
	      PRINTSeqposcast(repeatptr->pos1 + repeatptr->len - 1 
	                - offset));
    printf("boundaries->rightLTR_5 = " FormatSeqpos "\n", 
	      PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset 
	                - offset));
    printf("boundaries->rightLTR_3 = " FormatSeqpos "\n", 
	      PRINTSeqposcast(repeatptr->pos1 + repeatptr->offset + 
	                repeatptr->len - 1 - offset));
*/
    /* store new boundaries-positions in boundaries */
    adjustboundariesfromXdropextension( 
	           xdropbest_left, 
		   xdropbest_right,
	           repeatptr->pos1,  /*seed1 startpos*/
		   repeatptr->pos1 + repeatptr->offset, /*seed2 startpos*/
                   repeatptr->pos1 + repeatptr->len - 1, /*seed1 endpos*/
                   repeatptr->pos1 + repeatptr->offset + repeatptr->len - 1,
		                                         /*seed2 endpos*/  
		   boundaries,
		   offset); //muss noch uerbergeben werden
/*  
    // if search for motif and/or TSD
    if( motif->allowedmismatches < (Ushort)4 || minlengthTSD > (unsigned int) 1)
    {
      if( findcorrectboundaries(
            minlengthTSD,
	    maxlengthTSD,
	    boundaries,
	    virtualtree, 
	    vicinityforcorrectboundaries,
	    motif) != 0 )
      {
        return (int) -1;
      }

      // if search for TSDs and (not) motif
      if( boundaries->tsd && 
	  (motif->allowedmismatches >= (Ushort)4 ||
	   (boundaries->motif_near_tsd && boundaries->motif_far_tsd))
        )
      {
	// predicted as full LTR-pair, keep it
      }
      else
      {
	// if search for motif only (and not TSD)
	if( minlengthTSD <= (unsigned int) 1 &&
	    boundaries->motif_near_tsd &&
	    boundaries->motif_far_tsd )
	{
	  // predicted as full LTR-pair, keep it
	}
	else
	{
	  // delete this LTR-pair candidate
	  arrayLTRboundaries->nextfreeLTRboundaries--;
	  continue;
	}
      }
    }
#ifdef DEBUG
    shownewboundariesifadjusted(boundaries, multiseqoffset); 
#endif

    // check length and distance constraints again
    if( checklengthanddistanceconstraints(boundaries,
					  repeatinfo) != 0)
    {
      // delete this LTR-pair candidate
      arrayLTRboundaries->nextfreeLTRboundaries--;
      continue;
    }

    // check similarity from candidate pair
    if( determinesimilarityfromcandidate(
	         boundaries,
                 virtualtree) != 0)
    {
      return (int) -1;
    }
    if( boundaries->similarity < similaritythreshold )
    {
      // delete this LTR-pair candidate
      arrayLTRboundaries->nextfreeLTRboundaries--;
    }
*/
  }

  return 0;
}
