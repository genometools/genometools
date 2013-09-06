#!/bin/sh -ex

# build manuals for distributions first
make cleanup
make manuals

# Darwin
# https://launchpad.net/~flosoft/+archive/cross-apple
make clean
make SYSTEM=Darwin                     \
     MACHINE=i386                      \
     CC=i686-apple-darwin10-gcc        \
     AR=i686-apple-darwin10-ar         \
     RANLIB=i686-apple-darwin10-ranlib \
     STRIP=i686-apple-darwin10-strip   \
     cairo=no                          \
     curses=no                         \
     dist $*
