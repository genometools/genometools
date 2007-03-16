/* $Id: ll.h,v 1.12 2004/03/13 13:28:02 dvd Exp $ */

#ifndef LL_H
#define LL_H 1

/* all limits that can affect speed or memory consumption;
 prefixes correspond to module names
 */

#define RN_LEN_P 1024
#define RN_PRIME_P 0x3fd
#define RN_LIM_P (4*RN_LEN_P)
#define RN_LEN_NC 256
#define RN_PRIME_NC 0xfb
#define RN_LEN_S 256

#define SC_LEN 64

#define RND_LEN_F 1024

#define DRV_LEN_DTL 4
#define DRV_LEN_M 4096
#define DRV_PRIME_M 0xffd
#define DRV_LIM_M (8*DRV_LEN_M)

#define RNX_LEN_EXP 16
#define RNX_LIM_EXP 64

#define XCL_LEN_T 1024
#define XCL_LIM_T 16384

#define RX_LEN_P 256
#define RX_PRIME_P 0xfb
#define RX_LIM_P (4*RX_LEN_P)
#define RX_LEN_R 32
#define RX_PRIME_R 0x1f
#define RX_LEN_2 RX_PRIME_R
#define RX_PRIME_2 RX_PRIME_R
#define RX_LEN_M 1024
#define RX_PRIME_M 0x3fd
#define RX_LIM_M (8*RX_LEN_M)

#endif
