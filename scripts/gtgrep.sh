#!/bin/sh

egrep "$@" src/match/*\
           src/ltr/*\
           src/mgth/*\
           src/tools/*.[ch]\
           src/*.[ch]\
           scripts/*\
           testsuite/*.rb
