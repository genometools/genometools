#!/bin/sh -ex

make cleanup
make CC=i686-w64-mingw32-gcc AR=i686-w64-mingw32-ar fpic=no
