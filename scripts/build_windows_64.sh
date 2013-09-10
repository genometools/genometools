#!/bin/sh -ex

make SYSTEM=Windows                                             \
     MACHINE=i686                                               \
     64bit=yes                                                  \
     CC=x86_64-w64-mingw32-gcc                                  \
     AR=x86_64-w64-mingw32-ar                                   \
     fpic=no                                                    \
     cairo=no                                                   \
     sharedlib=no                                               \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter -DSQLITE_MALLOCSIZE=_msize' \
     $*
