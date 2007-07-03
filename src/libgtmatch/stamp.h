#ifndef STAMP_H
#define STAMP_H

#define STAMP\
        printf("STAMP(%d,%s)\n",__LINE__,__FILE__);\
        (void) fflush(stdout)

#ifndef STAMP
#define STAMP
#endif

#endif
