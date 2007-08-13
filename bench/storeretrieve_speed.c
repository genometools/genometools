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

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */
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

extern uint32_t
retrieveArrayElem(uint32_t a[], size_t numElem);
extern void
storeArrayElem(uint32_t a[], size_t numElem, uint32_t val);


int
main(int argc, char *argv[])
{
  uint32_t *numStore; /*< used for random ints read back */
  size_t i, numRnd;
  int numBits = sizeof(numStore[0]) * CHAR_BIT;
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
  if(!(numStore = malloc(sizeof(uint32_t) *numRnd)))
  {
    perror("Storage allocations failed");
    return EXIT_FAILURE;
  }
  gettimeofday(&startTime, NULL);
  for(i = 0; i < numRnd; ++i)
  {
    uint32_t v = random();
    storeArrayElem(numStore, i, v);
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Storing %lu random numbers at %d bits precision sequentially"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }
  gettimeofday(&startTime, NULL);  
  for(i = 0; i < numRnd; ++i)
  {
    retrieveArrayElem(numStore, i);
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Retrieving %lu random numbers at %d bits precision sequentially"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }
  gettimeofday(&startTime, NULL);  
  for(i = 0; i < numRnd; ++i)
  {
    retrieveArrayElem(numStore, random()%numRnd);
  }
  gettimeofday(&endTime, NULL);
  {
    double timeDiff = ((double)endTime.tv_sec + (double)endTime.tv_usec * 1e-6)
      - ((double)startTime.tv_sec + (double)startTime.tv_usec * 1e-6);
    printf("Retrieving %lu random numbers at %d bits precision randomly"
           " required %f seconds\n", (unsigned long)numRnd, numBits, timeDiff);
  }
  free(numStore);
  return EXIT_SUCCESS;
}

