#! /bin/sh
find src/match -type f \( -name '*.[ch]' -o -name '*.pr' \) -print | \
  etags -o TAGS.match -
find . -type f \( -name '*.[ch]' -o -name '*.pr' \) -print | \
  etags -
