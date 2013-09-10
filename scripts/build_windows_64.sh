#!/bin/sh -ex

make CC=x86_64-w64-mingw32-gcc                                  \
     AR=x86_64-w64-mingw32-ar                                   \
     fpic=no                                                    \
     cairo=no                                                   \
     with-sqlite=no                                             \
     sharedlib=no                                               \
     64bit=yes                                                  \
     SYSTEM=Windows                                             \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter' \
     $*
