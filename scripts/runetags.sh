#! /bin/sh
find . -type f -name '*.[ch]' -print | \
  etags -