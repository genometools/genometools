#!/bin/sh -ex

make CC=x86_64-w64-mingw32-gcc                                  \
     AR=x86_64-w64-mingw32-ar                                   \
     fpic=no                                                    \
     cairo=no                                                   \
     curses=no                                                  \
     with-sqlite=no                                             \
     sharedlib=no                                               \
     SYSTEM=Windows                                             \
     CFLAGS='-Wno-error=attributes -Wno-error=unused-parameter' \
     $*
