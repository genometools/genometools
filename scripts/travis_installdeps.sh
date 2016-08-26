#!/bin/bash

echo $TRAVIS_OS_NAME
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  apt-get update -q
  apt-get install ncbi-blast+ gcc-multilib -y
  apt-get install binutils-mingw-w64-i686 gcc-mingw-w64-i686 -y
  apt-get install binutils-mingw-w64-x86-64 gcc-mingw-w64-x86-64 -y
else
  brew update
  brew install pango cairo
  brew install homebrew/science/blast
fi
