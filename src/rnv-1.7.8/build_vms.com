$! Replace these two strings to indicate where your expat install is located.
$	expat_headers = "my_disk:[ref.c_include]"
$	expat_olb = "my_disk:[olb]expat.olb"
$
$	create vms.h
$	deck
#ifndef VMS_H
#define VMS_H
#define UNISTD_H <unistd.h>
#define EXPAT_H "expat.h"
#define RNV_VERSION "1.7.7"
#define ARX_VERSION "1.7.7"
#define RVP_VERSION "1.7.7"
#endif /* VMS_H */
$	eod
$
$	cc := cc/first_include=vms.h/incl='expat_headers'
$	modules = "XCL,RNV,ARY,DRV,DSL,DXL,ER,HT,M,RN,RNC,RND,RNL,RNX," + -
	   	  "RX,RX_CLS_RANGES,RX_CLS_U,S,SC,U,XMLC,XSD,XSD_TM"
$
$	library/create/object rnv.olb
$	count = 0
$loop:
$	module = f$element (count, ",", modules)
$	if module .eqs. ","
$	then
$	    goto end_loop
$	endif
$	cc 'module'
$	library/object/insert rnv.olb 'module'
$	delete/nolog 'module'.obj;*
$	count = count + 1
$	goto loop
$end_loop:
$	link/exe=rnv.exe rnv/lib/include=xcl,'expat_olb'/lib
$
$! Now for the supporting cast...
$	cc test
$	link test,rnv/lib
$	cc arx
$	link arx,rnv/lib,'expat_olb'/lib
$	cc rvp
$	link rvp,rnv/lib
$	purge/nolog *.exe
$	purge/nolog *.olb
$	delete/nolog vms.h;*
$	exit
