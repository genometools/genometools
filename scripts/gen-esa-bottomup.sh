#!/bin/sh

TEMPLATE=src/match/esa-bottomup
SC=scripts/gen-esa-bottomup.rb

${SC} --key spmsk --nobranch > ${TEMPLATE}-spmsk.inc
${SC} --key spmeq > ${TEMPLATE}-spmeq.inc
${SC} --key spmvar > ${TEMPLATE}-spmvar.inc
${SC} --key rdjcv --reader --absolute > ${TEMPLATE}-rdjcv.inc
${SC} --key rdjce --reader --absolute > ${TEMPLATE}-rdjce.inc
${SC} --key shulen --reader --absolute > ${TEMPLATE}-shulen.inc
${SC} --key maxpairs --reader --absolute > ${TEMPLATE}-maxpairs.inc
${SC} --key errfind --reader --absolute > ${TEMPLATE}-errfind.inc
