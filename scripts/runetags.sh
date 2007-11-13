#! /bin/sh
find . -type f \( -name '*.[ch]' -o -name '*.pr' \) -print | \
  etags -