#!/bin/sh

egrep "$@" src/libgtmatch/* src/libgtltr/* src/*.[ch] src/tools/*.[ch] scripts/*
