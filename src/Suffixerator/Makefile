CC=ccache gcc
LD=${CC}
VPATH=LIB
GTHDIR=${HOME}/genometools
GTHLIB=${GTHDIR}/lib
LIBGTHCOREBZ2=${GTHLIB}/libgtcore.a
INCLUDEDIR=-I${GTHDIR}/include -I${VPATH}
WARN=-Wstrict-prototypes -Wundef -Wshadow -Wcast-align -Wsign-compare -Wnested-externs -Wall -Winline -Werror
CFLAGS=-m32 -g -O3 ${WARN} ${INCLUDEDIR} -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SPLINTFLAGS=-DDEBUG -f Splintoptions ${INCLUDEDIR}
LDFLAGS=-m32
LDLIBS=-lm -lz -lbz2
SPLINTALL=${subst .c,.splint,${wildcard *.c} ${wildcard LIB/*.c}}
LIBARCHIVE=mylib.a
LIBFILES=${subst LIB/,,${subst .c,.o,${wildcard LIB/*.c}}}

SUFFIXERATOROBJ=suffixerator-mn.o ${LIBARCHIVE}

GUESSPROTOBJ=guessprot-mn.o ${LIBARCHIVE}

SFXMAPOBJ=sfxmap-mn.o ${LIBARCHIVE}

all:suffixerator.x guessprot.x sfxmap.x splintall

libcfiles:
	echo ${LIBFILES}

${LIBARCHIVE}:${LIBFILES}
	ar ruv $@ ${LIBFILES}
	ranlib.sh $@

suffixerator.x:${SUFFIXERATOROBJ} ${LIBGTHCOREBZ2}
	$(LD) $(LDFLAGS) ${SUFFIXERATOROBJ} ${LIBGTHCOREBZ2} ${LDLIBS} -o $@

guessprot.x:${GUESSPROTOBJ}
	$(LD) $(LDFLAGS) ${GUESSPROTOBJ} ${LDLIBS} -o $@

sfxmap.x:${SFXMAPOBJ} ${LIBGTHCOREBZ2}
	 $(LD) $(LDFLAGS) ${SFXMAPOBJ} ${LIBGTHCOREBZ2} ${LDLIBS} -o $@

%.o:%.c
	${CC} ${CFLAGS} -c $< -o $@

%.prepro:%.c
	${CC} -E -g3 ${CFLAGS} -DDEBUG -c $< -o $@
	indent $@
	rm -f ${@}~

%.pr:%.c
	skproto.x $< > $@

%.splint:%.c
	splint ${SPLINTFLAGS} $<
	touch $*.splint

splintall:${SPLINTALL}

splintclean:
	rm -f *.splint LIB/*.splint

cflagsstring:
	@echo "${CFLAGS}"

-include Dependencies.mf

check:all Random.fna Atinsert.fna RandomN.fna
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

Random.fna:
	echo ">" > $@
	oneFE.rand.x "acn" 10000
	fold --bytes --spaces --width=70 -s R10000acn >> $@
	echo "" >> $@
	rm R10000acn

RandomN.fna:
	echo ">" > $@
	oneFE.rand.x "n" 10000
	fold --bytes --spaces --width=70 -s R10000n >> $@
	echo "" >> $@

Atinsert.fna:
	head -n 200 ${AT} | tr G N > Atinsert.fna

clean:cleanpp splintclean
	rm -f *.o *.x *.splint R10000acn R10000n *.prepro TMP.* mylib.a

cleanpp:
	cleanpp.sh
