#!/bin/sh -ex

# build manuals for distributions first
make cleanup
make cairo=no manuals

# Windows 32-bit
make clean
make SYSTEM=Windows                                             \
     MACHINE=i686                                               \
     32bit=yes                                                  \
     CC=i686-w64-mingw32-gcc                                    \
     AR=i686-w64-mingw32-ar                                     \
     fpic=no                                                    \
     cairo=no                                                   \
     sharedlib=no                                               \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter -DSQLITE_MALLOCSIZE=_msize' \
     dist $*

# Windows 64-bit
make clean
make SYSTEM=Windows                                             \
     MACHINE=i686                                               \
     64bit=yes                                                  \
     CC=x86_64-w64-mingw32-gcc                                  \
     AR=x86_64-w64-mingw32-ar                                   \
     fpic=no                                                    \
     cairo=no                                                   \
     sharedlib=no                                               \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter -DSQLITE_MALLOCSIZE=_msize' \
     dist $*
