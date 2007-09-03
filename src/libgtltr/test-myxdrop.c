#include <stdio.h>
#include <stdbool.h>

#include "libgtcore/str.h"
#include "libgtmatch/alphadef.h"

#include "myxdrop.h"

int testmyxdrop(int argc, char** argv, Env *env)
{
  Str *useq, 
      *vseq;
  unsigned long ulen, vlen;
  bool showfilledmatrix = True;
  Myxdropbest xdropbest,
	      xdropbest_vglkurtz;
  ArrayMyfrontvalue fronts;
  int distance_right = 0;

  useq = str_new(env);
  vseq = str_new(env);

  str_set(useq, "FREIZEIT", env);
  str_set(vseq, "ZEITGEIST", env);

  /** userdefined options **/
  /* arbitrary match, mismatch and indel/del scores */
  Myxdropscore xdropbelowscore;   //userdefined
  Arbitraryscores arbitscores;
  arbitscores.mat  = (int) 2,
    arbitscores.mis  = (int) -1,
    arbitscores.ins  = (int) -2,
    arbitscores.del  = (int) -2;
  arbitscores.gcd  = (int) 1; // set only for initialization
  printf("scores:\n");
  printf("mat  = %ld\n", arbitscores.mat);
  printf("mis  = %ld\n", arbitscores.mis);
  printf("ins  = %ld\n", arbitscores.ins);
  printf("del  = %ld\n", arbitscores.del);

  //CHECKARGNUM (4, "seq1 seq2 xdopbelowscore");

  xdropbelowscore = (int) 5;

  ulen = str_length(useq);
  vlen = str_length(vseq);

  printf("alignment sequence useq:\n%s\n", str_get(useq));
  printf("alignment sequence vseq:\n%s\n", str_get(vseq));

  INITARRAY (&fronts, Myfrontvalue);

  evalxdroparbitscoresright
    (&arbitscores, &xdropbest, &fronts, useq, vseq, ulen, vlen,
     xdropbelowscore, &distance_right);

  if (showfilledmatrix)
  {
    if (showmatrix (&fronts, distance_right, useq, vseq, ulen, vlen) != 0)
    {
      STANDARDMESSAGE;
    }
  }

  // compare with kurtz implementation
  evaleditxdropright(&xdropbest_vglkurtz, useq, ulen, vseq, vlen, xdropbelowscore);
  if(xdropbest.score != xdropbest_vglkurtz.score
      || xdropbest.ivalue != xdropbest_vglkurtz.ivalue
      || xdropbest.jvalue != xdropbest_vglkurtz.jvalue
    )
  {
    fprintf(stderr, "comparison with kurtz implementation: Failure, best (prefix-) scores are different.\n");
    fprintf(stderr, "my score:  %ld             kurtz score: %ld\n", (Showsint) xdropbest.score,
	(Showsint) xdropbest_vglkurtz.score);
    fprintf(stderr, "my best prefix: (%lu,%lu)  kurtz best prefix: (%lu,%lu)\n",
	(Showuint) xdropbest.ivalue, (Showuint) xdropbest.jvalue,
	(Showuint) xdropbest_vglkurtz.ivalue, (Showuint) xdropbest_vglkurtz.ivalue);
    return EXIT_FAILURE;
  }
  else
  {
    fprintf(stderr, "comparison with kurtz implementation: OK, best (prefix-) scores are equal.\n");
    fprintf(stderr, "my score:  %ld             kurtz score: %ld\n", (Showsint) xdropbest.score,
	(Showsint) xdropbest_vglkurtz.score);
    fprintf(stderr, "my best prefix: (%lu,%lu)  kurtz best prefix: (%lu,%lu)\n",
	(Showuint) xdropbest.ivalue, (Showuint) xdropbest.jvalue,
	(Showuint) xdropbest_vglkurtz.ivalue, (Showuint) xdropbest_vglkurtz.ivalue);
  }


  FREEARRAY (&fronts, Frontvalue);
  FREESPACE (useq);
  FREESPACE (vseq);

  checkspaceleak ();

  return EXIT_SUCCESS;
}
