#!/bin/sh

LIBDIR=$WORKVSTREE/lib/x86_64-unknown-linux-gnu/32bit

extractsyms.sh ${LIBDIR}/libkurtz.a\
               ${LIBDIR}/libkurtz-basic.a\
               ${LIBDIR}/libmkvtree.a\
               ${LIBDIR}/libvmengine.a | sort -u > VMNAMES

extractsyms.sh lib/libgenometools.a | sort -u > GTNAMES

echo "========= the folllowing name clashes have been identified ======"
comm -1 -2 GTNAMES VMNAMES
