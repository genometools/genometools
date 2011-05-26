AUTOMAKE_OPTIONS = 1.10.1 foreign dist-bzip2 dist-zip dist-lzma

bin_PROGRAMS = rnv arx rvp xsdck
check_PROGRAMS = rnvtest
noinst_LIBRARIES = librnv1.a librnv2.a

TESTS = rnvtest

man_MANS = rnv.1



librnv1_a_SOURCES = \
	ary.h ary.c \
	drv.h drv.c \
	er.h er.c \
	erbit.h \
	ht.h ht.c \
	m.h m.c \
	rn.h rn.c \
	rnc.h rnc.c \
	rnd.h rnd.c \
	rnl.h rnl.c \
	rnv.h rnv.c \
	rnx.h rnx.c \
	rx.h rx.c \
	s.h s.c \
	sc.h sc.c \
	u.h u.c \
	xmlc.h xmlc.c \
	xsd.h xsd.c \
	xsd_tm.h xsd_tm.c



librnv2_a_SOURCES = \
	dxl.h dxl.c \
	dsl.h dsl.c



rnv_CPPFLAGS = -DRNV_VERSION=\"@VERSION@\"
rnv_LDADD = librnv1.a librnv2.a

rnv_SOURCES = \
	ll.h \
	xcl.c



arx_CPPFLAGS = -DARX_VERSION=\"@VERSION@\"
arx_LDADD = librnv1.a

arx_SOURCES = \
	arx.c



rvp_CPPFLAGS = -DRVP_VERSION=\"@VERSION@\"
rvp_LDADD = librnv1.a librnv2.a

rvp_SOURCES = \
	rvp.c



xsdck_LDADD = librnv1.a

xsdck_SOURCES = \
	xsdck.c



rnvtest_LDADD = librnv1.a

rnvtest_SOURCES = \
	test.c



EXTRA_DIST = \
	scm/dsl.scm \
	scm/M.scm \
	scm/rx.scm \
	scm/rx-ranges.scm \
	scm/s48.scm \
	scm/spat.scm \
	scm/u.scm \
	scm/xml-ranges.scm \
	\
	tools/addr-spec.rnc \
	tools/addr-spec-dsl.rnc \
	tools/arx.conf \
	tools/b64gen.c \
	tools/gn.sed \
	tools/pr.c \
	tools/rgen.pl \
	tools/rnv.vim \
	tools/rvp.pl \
	tools/rvp.py \
	tools/xck \
	tools/xsd2rnc.xsl \
	tools/xslt-dsl.rnc \
	\
	readme.txt \
	readme32.txt \
	rnv.1 \
	src.txt \
	\
	build_vms.com \
	\
	rx_cls_ranges.c \
	rx_cls_u.c

## NOTE: The two source files above are included by rx.c

