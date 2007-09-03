default:
	./configure
	$(MAKE)

distclean:
	./configure --without-man-pages
	$(MAKE) $@

normal reentrant demos demos_r clean install_lib install_bin install_inc \
 install_man install:
	./configure
	$(MAKE) $@
