CC=gcc
LD=${CC}
#define GTDIR appropriately
CURL=${PROJECT}/curl
CURLLIBS=-L${CURL}/lib -lcurl
CURLINCLUDE=-I${CURL}/include
INCLUDES=-I${GTDIR}/src ${CURLINCLUDE} -I.
CFLAGS=-m32 -O3 -Wall -Werror -g ${INCLUDES}
LDFLAGS=-L${GTDIR}/lib -m32

LDCURLLIBS=${CURLLIBS}

LDLIBS=-lgtcore\
       -lgtext\
       -lexpat\
       -lbz2\
       -lgz\
       -lm

OBJECTS=metagenome.o\
        metagenomethreader.o\
        mg_reverse.o\
        mg_codon.o\
        mg_xmlparser.o\
        mg_combinedscore.o\
        mg_computepath.o\
        mg_compute_gene_prediction.o\
        mg_outputwriter.o

OBJECTSCURL=mg_curl.o

mgth: ${OBJECTS}
	${LD} ${OBJECTS} ${LDFLAGS} ${LDLIBS} -o $@

curl: ${OBJECTSCURL} ${OBJECTS}
	${LD} ${OBJECTSCURL} ${OBJECTS} ${LDFLAGS} ${LDCURLLIBS} ${LDLIBS} -o mgth

splint:
	splint -f ${GTDIR}/testdata/Splintoptions -DNDEBUG ${INCLUDES} *.c

clean:
	rm -f *.o

cleanup:clean
	rm -f mgth
