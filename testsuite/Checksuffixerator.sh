#!/bin/sh -ex

Checkall.sh Random-Small.fna Random.fna Atinsert.fna RandomN.fna TTT-small.fna
Cmpdbfile.sh Random-Small.fna
Cmpdbfile.sh Random.fna
Cmpdbfile.sh Atinsert.fna Random.fna
Cmpdbfile.sh TTT-small.fna
Checkmapped.sh -parts 1 -pl 3 ${SWK} ${SW}
Checkmapped.sh -parts 2 -pl 3 ${SWK} ${SW}
Checkmapped.sh -parts 3 -pl 3 ${SWK} ${SW}
Cmpdbfile.sh ${ATK} ${AT} ${GRUMBACH}/*.fna
Checkmapped.sh -parts 1 -pl 3 Random-Small.fna
Checkmapped.sh -parts 1 -pl 3 Random.fna
Checkmapped.sh -parts 1 -pl 3 RandomN.fna
Checkmapped.sh -parts 1 -pl 3 RandomN.fna Random.fna
Checkmapped.sh -parts 1 -pl 7 -smap TransDNA ${AT}
Checkmapped.sh -parts 1 -pl 3 Atinsert.fna RandomN.fna Random.fna
Checkmapped.sh -parts 1 -pl 10 Random.fna Atinsert.fna
Checkmapped.sh -parts 2 -pl 10 Random.fna Atinsert.fna
Checkmapped.sh -parts 3 -pl 10 Random.fna Atinsert.fna
Checkmapped.sh -parts 1 -pl 8 ${ATK} ${AT} ${GRUMBACH}/*.fna
Checkmapped.sh -parts 2 -pl 8 ${ATK} ${AT} ${GRUMBACH}/*.fna
Checkmapped.sh -parts 3 -pl 8 ${ATK} ${AT} ${GRUMBACH}/*.fna
Checkmapped.sh -parts 1 -pl 3 -smap TransProt11 ${SWK} ${SW}

