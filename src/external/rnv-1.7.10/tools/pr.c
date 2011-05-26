/* $Id: pr.c 153 2003-12-24 10:52:24Z dvd $ */

#include <stdlib.h>
#include <stdio.h>

int main(int argc,char **argv) {
  int *primes;
  int N,i,step;
  primes=calloc(N=atoi(*(++argv)),sizeof(int));
  for(i=0;i!=N;++i) primes[i]=1;
  for(step=2;step!=256;++step) {
    i=step+step;
    for(;;) {
      if(i>=N) break;
      primes[i]=0;
      i+=step;
    }
  }
  for(i=0;i!=N;++i) if(primes[i]) printf("%x\n",i);
}
