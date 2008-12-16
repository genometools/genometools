#!/bin/sh
find testdata ${GTTESTDATA} -type f \( -name '*.fsa.gz' -o\
                            -name '*.fna' -o\
                            -name '*.fsa' \)
