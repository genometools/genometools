#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "types.h"

#include "guessprot.pr"

#define FPRINTUSAGE\
        fprintf(stderr,"Usage: %s <filename>\n",argv[0])

int main(int argc,const char *argv[])
{
  if (argc != 2)
  {
    FPRINTUSAGE;
    return EXIT_FAILURE;
  }
  if (guessifproteinsequencestream(argv[1]))
  {
    return 1;
  } else
  {
    return 0;
  }
}
