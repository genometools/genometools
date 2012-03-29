#!/bin/sh

gen-radixsort-ip.rb > src/core/radixsort-ip-ulong.inc
gen-radixsort-ip.rb --ulongpair > src/core/radixsort-ip-ulongpair.inc
