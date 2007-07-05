#!/bin/sh

egrep $* src/libgtmatch/* src/*.[ch] src/tools/*.[ch]
