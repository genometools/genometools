/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/bitpackstring.h>

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 100000000,
};

static inline void
timePrefixPrint(FILE *fp)
{
  time_t t = time(0);
  fprintf(fp, "%.24s", asctime(localtime(&t)));
}


int
main(int argc, char *argv[])
{
  BitElem *bitStore;
  size_t i, numRnd, offsetStart, offset;
  int  numBits;
  unsigned long seedval;
  struct timeval startTime, endTime;
  if(argc < 2)
  {
    struct timeval seed;
    gettimeofday(&seed, NULL);
/*     srandom(seedval = seed.tv_sec + seed.tv_usec); */
    srandom(seedval = 1182618578);
  }
  else
  {
    srandom(seedval = strtol(argv[1], NULL, 0));
  }
  fprintf(stderr, "seedval = %lu\n", seedval);
  if(argc < 3)
  {
    numRnd = 100000000;
  }
  else
  {
    numRnd = strtol(argv[2], NULL, 0);
  }
  offset = offsetStart = 0;
  numBits = 10;
  if(!(bitStore = 
       malloc(bitElemsAllocSize(sizeof(uint32_t)
                                *CHAR_BIT*numRnd) * sizeof(BitElem))))
  {
    perror("Storage allocations failed");
    return EXIT_FAILURE;
  }
  gettimeofday(&startTime, NULL);
  for(i = 0; i < numRnd; ++i)
  {
    uint32_t v = random();
    bsStoreUInt(bitStore, offset, numBits, v);
    offset += numBits;
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Storing %lu random numbers at %d bits precision sequentially"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }
  offset = offsetStart;
  gettimeofday(&startTime, NULL);  
  for(i = 0; i < numRnd; ++i)
  {
    bsGetUInt(bitStore, offset, numBits);
    offset += numBits;
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Retrieving %lu random numbers at %d bits precision sequentially,"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }
  gettimeofday(&startTime, NULL);  
  for(i = 0; i < numRnd; ++i)
  {
    bsGetUInt(bitStore, numBits*(random()%numRnd), numBits);
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Retrieving %lu random numbers at %d bits precision randomly,"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }

  free(bitStore);
  return EXIT_SUCCESS;
}

