#!/bin/sh -ex

make CC=i686-w64-mingw32-gcc                                    \
     AR=i686-w64-mingw32-ar                                     \
     fpic=no                                                    \
     cairo=no                                                   \
     with-sqlite=no                                             \
     sharedlib=no                                               \
     SYSTEM=Windows                                             \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter' \
     $*
