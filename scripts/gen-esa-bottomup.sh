#!/bin/sh

TEMPLATE=src/match/esa-bottomup
SC=scripts/gen-esa-bottomup.rb

${SC} --key maxpairs --reader --absolute \
                     --no_process_lcpinterval > ${TEMPLATE}-maxpairs.inc
${SC} --key spmsk --no_process_branchingedge > ${TEMPLATE}-spmsk.inc
${SC} --key rdjcv --reader --absolute \
                  --no_process_lcpinterval > ${TEMPLATE}-rdjcv.inc
${SC} --key rdjce --reader --absolute \
                  --no_process_lcpinterval > ${TEMPLATE}-rdjce.inc
${SC} --key errfind --reader --absolute > ${TEMPLATE}-errfind.inc
${SC} --key spmeq --no_process_lcpinterval > ${TEMPLATE}-spmeq.inc
${SC} --key spmvar > ${TEMPLATE}-spmvar.inc
${SC} --key shulen --reader --absolute \
                   --no_process_lcpinterval > ${TEMPLATE}-shulen.inc
${SC} --key shulen --gtlcpvaluetypeset --absolute --no_process_lastvalue \
                   --no_process_lcpinterval \
                   --withlastfrompreviousbucket \
                   --additionaluint32bucket \
                   --no_declarations > ${TEMPLATE}-shulen-RAM.inc
