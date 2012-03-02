/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifndef TAB_DIFFCOVER_H
#define TAB_DIFFCOVER_H

static Diffvalue differencecovertab[] = {
  /*
     1
   */ UScast (0),
  /*
     2
   */ UScast (0), UScast (1),
  /*
     4
   */ UScast (0), UScast (1), UScast (2),
  /*
     8
   */ UScast (0), UScast (1), UScast (2), UScast (4),
  /*
     16
   */ UScast (0), UScast (1), UScast (2), UScast (5), UScast (8),
  /*
     32
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (7),
    UScast (11), UScast (19),
  /*
     64
   */ UScast (0), UScast (1), UScast (2), UScast (5),
    UScast (14), UScast (16), UScast (34), UScast (42),
    UScast (59),
  /*
     128
   */ UScast (0), UScast (1), UScast (3), UScast (7),
    UScast (17), UScast (40), UScast (55), UScast (64),
    UScast (75), UScast (85), UScast (104), UScast (109),
    UScast (117),
  /*
     256
   */ UScast (0), UScast (1), UScast (3), UScast (7),
    UScast (12), UScast (20), UScast (30), UScast (44),
    UScast (65), UScast (80), UScast (89), UScast (96),
    UScast (114), UScast (122), UScast (128), UScast (150),
    UScast (196), UScast (197), UScast (201), UScast (219),
  /*
     512
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (9), UScast (18), UScast (27), UScast (36),
    UScast (45), UScast (64), UScast (83), UScast (102),
    UScast (121), UScast (140), UScast (159), UScast (178),
    UScast (197), UScast (216), UScast (226), UScast (236),
    UScast (246), UScast (256), UScast (266), UScast (267),
    UScast (268), UScast (269), UScast (270),
  /*
     1024
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (13), UScast (26),
    UScast (39), UScast (52), UScast (65), UScast (78),
    UScast (91), UScast (118), UScast (145), UScast (172),
    UScast (199), UScast (226), UScast (253), UScast (280),
    UScast (307), UScast (334), UScast (361), UScast (388),
    UScast (415), UScast (442), UScast (456), UScast (470),
    UScast (484), UScast (498), UScast (512), UScast (526),
    UScast (540), UScast (541), UScast (542), UScast (543),
    UScast (544), UScast (545), UScast (546),
  /*
     2048
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (7), UScast (8), UScast (9),
    UScast (19), UScast (38), UScast (57), UScast (76),
    UScast (95), UScast (114), UScast (133), UScast (152),
    UScast (171), UScast (190), UScast (229), UScast (268),
    UScast (307), UScast (346), UScast (385), UScast (424),
    UScast (463), UScast (502), UScast (541), UScast (580),
    UScast (619), UScast (658), UScast (697), UScast (736),
    UScast (775), UScast (814), UScast (853), UScast (892),
    UScast (931), UScast (951), UScast (971), UScast (991),
    UScast (1011), UScast (1031), UScast (1051), UScast (1071),
    UScast (1091), UScast (1111), UScast (1131), UScast (1132),
    UScast (1133), UScast (1134), UScast (1135), UScast (1136),
    UScast (1137), UScast (1138), UScast (1139), UScast (1140),
  /*
     4096
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (7), UScast (8), UScast (9),
    UScast (10), UScast (11), UScast (12), UScast (13),
    UScast (27), UScast (54), UScast (81), UScast (108),
    UScast (135), UScast (162), UScast (189), UScast (216),
    UScast (243), UScast (270), UScast (297), UScast (324),
    UScast (351), UScast (378), UScast (433), UScast (488),
    UScast (543), UScast (598), UScast (653), UScast (708),
    UScast (763), UScast (818), UScast (873), UScast (928),
    UScast (983), UScast (1038), UScast (1093), UScast (1148),
    UScast (1203), UScast (1258), UScast (1313), UScast (1368),
    UScast (1423), UScast (1478), UScast (1533), UScast (1588),
    UScast (1643), UScast (1698), UScast (1753), UScast (1808),
    UScast (1863), UScast (1891), UScast (1919), UScast (1947),
    UScast (1975), UScast (2003), UScast (2031), UScast (2059),
    UScast (2087), UScast (2115), UScast (2143), UScast (2171),
    UScast (2199), UScast (2227), UScast (2255), UScast (2256),
    UScast (2257), UScast (2258), UScast (2259), UScast (2260),
    UScast (2261), UScast (2262), UScast (2263), UScast (2264),
    UScast (2265), UScast (2266), UScast (2267), UScast (2268),
  /*
     8192
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (7), UScast (8), UScast (9),
    UScast (10), UScast (11), UScast (12), UScast (13),
    UScast (14), UScast (15), UScast (16), UScast (17),
    UScast (18), UScast (37), UScast (74), UScast (111),
    UScast (148), UScast (185), UScast (222), UScast (259),
    UScast (296), UScast (333), UScast (370), UScast (407),
    UScast (444), UScast (481), UScast (518), UScast (555),
    UScast (592), UScast (629), UScast (666), UScast (703),
    UScast (778), UScast (853), UScast (928), UScast (1003),
    UScast (1078), UScast (1153), UScast (1228), UScast (1303),
    UScast (1378), UScast (1453), UScast (1528), UScast (1603),
    UScast (1678), UScast (1753), UScast (1828), UScast (1903),
    UScast (1978), UScast (2053), UScast (2128), UScast (2203),
    UScast (2278), UScast (2353), UScast (2428), UScast (2503),
    UScast (2578), UScast (2653), UScast (2728), UScast (2803),
    UScast (2878), UScast (2953), UScast (3028), UScast (3103),
    UScast (3178), UScast (3253), UScast (3328), UScast (3403),
    UScast (3478), UScast (3516), UScast (3554), UScast (3592),
    UScast (3630), UScast (3668), UScast (3706), UScast (3744),
    UScast (3782), UScast (3820), UScast (3858), UScast (3896),
    UScast (3934), UScast (3972), UScast (4010), UScast (4048),
    UScast (4086), UScast (4124), UScast (4162), UScast (4200),
    UScast (4201), UScast (4202), UScast (4203), UScast (4204),
    UScast (4205), UScast (4206), UScast (4207), UScast (4208),
    UScast (4209), UScast (4210), UScast (4211), UScast (4212),
    UScast (4213), UScast (4214), UScast (4215), UScast (4216),
    UScast (4217), UScast (4218),
  /*
     16384
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (7), UScast (8), UScast (9),
    UScast (10), UScast (11), UScast (12), UScast (13),
    UScast (14), UScast (15), UScast (16), UScast (17),
    UScast (18), UScast (19), UScast (20), UScast (21),
    UScast (22), UScast (23), UScast (24), UScast (25),
    UScast (26), UScast (53), UScast (106), UScast (159),
    UScast (212), UScast (265), UScast (318), UScast (371),
    UScast (424), UScast (477), UScast (530), UScast (583),
    UScast (636), UScast (689), UScast (742), UScast (795),
    UScast (848), UScast (901), UScast (954), UScast (1007),
    UScast (1060), UScast (1113), UScast (1166), UScast (1219),
    UScast (1272), UScast (1325), UScast (1378), UScast (1431),
    UScast (1538), UScast (1645), UScast (1752), UScast (1859),
    UScast (1966), UScast (2073), UScast (2180), UScast (2287),
    UScast (2394), UScast (2501), UScast (2608), UScast (2715),
    UScast (2822), UScast (2929), UScast (3036), UScast (3143),
    UScast (3250), UScast (3357), UScast (3464), UScast (3571),
    UScast (3678), UScast (3785), UScast (3892), UScast (3999),
    UScast (4106), UScast (4213), UScast (4320), UScast (4427),
    UScast (4534), UScast (4641), UScast (4748), UScast (4855),
    UScast (4962), UScast (5069), UScast (5176), UScast (5283),
    UScast (5390), UScast (5497), UScast (5604), UScast (5711),
    UScast (5818), UScast (5925), UScast (6032), UScast (6139),
    UScast (6246), UScast (6353), UScast (6460), UScast (6567),
    UScast (6674), UScast (6781), UScast (6888), UScast (6995),
    UScast (7102), UScast (7156), UScast (7210), UScast (7264),
    UScast (7318), UScast (7372), UScast (7426), UScast (7480),
    UScast (7534), UScast (7588), UScast (7642), UScast (7696),
    UScast (7750), UScast (7804), UScast (7858), UScast (7912),
    UScast (7966), UScast (8020), UScast (8074), UScast (8128),
    UScast (8182), UScast (8236), UScast (8290), UScast (8344),
    UScast (8398), UScast (8452), UScast (8506), UScast (8560),
    UScast (8561), UScast (8562), UScast (8563), UScast (8564),
    UScast (8565), UScast (8566), UScast (8567), UScast (8568),
    UScast (8569), UScast (8570), UScast (8571), UScast (8572),
    UScast (8573), UScast (8574), UScast (8575), UScast (8576),
    UScast (8577), UScast (8578), UScast (8579), UScast (8580),
    UScast (8581), UScast (8582), UScast (8583), UScast (8584),
    UScast (8585), UScast (8586),
  /*
     32768
   */ UScast (0), UScast (1), UScast (2), UScast (3), UScast (4),
    UScast (5), UScast (6), UScast (7), UScast (8), UScast (9),
    UScast (10), UScast (11), UScast (12), UScast (13),
    UScast (14), UScast (15), UScast (16), UScast (17),
    UScast (18), UScast (19), UScast (20), UScast (21),
    UScast (22), UScast (23), UScast (24), UScast (25),
    UScast (26), UScast (27), UScast (28), UScast (29),
    UScast (30), UScast (31), UScast (32), UScast (33),
    UScast (34), UScast (35), UScast (36), UScast (37),
    UScast (75), UScast (150), UScast (225), UScast (300),
    UScast (375), UScast (450), UScast (525), UScast (600),
    UScast (675), UScast (750), UScast (825), UScast (900),
    UScast (975), UScast (1050), UScast (1125), UScast (1200),
    UScast (1275), UScast (1350), UScast (1425), UScast (1500),
    UScast (1575), UScast (1650), UScast (1725), UScast (1800),
    UScast (1875), UScast (1950), UScast (2025), UScast (2100),
    UScast (2175), UScast (2250), UScast (2325), UScast (2400),
    UScast (2475), UScast (2550), UScast (2625), UScast (2700),
    UScast (2775), UScast (2850), UScast (3001), UScast (3152),
    UScast (3303), UScast (3454), UScast (3605), UScast (3756),
    UScast (3907), UScast (4058), UScast (4209), UScast (4360),
    UScast (4511), UScast (4662), UScast (4813), UScast (4964),
    UScast (5115), UScast (5266), UScast (5417), UScast (5568),
    UScast (5719), UScast (5870), UScast (6021), UScast (6172),
    UScast (6323), UScast (6474), UScast (6625), UScast (6776),
    UScast (6927), UScast (7078), UScast (7229), UScast (7380),
    UScast (7531), UScast (7682), UScast (7833), UScast (7984),
    UScast (8135), UScast (8286), UScast (8437), UScast (8588),
    UScast (8739), UScast (8890), UScast (9041), UScast (9192),
    UScast (9343), UScast (9494), UScast (9645), UScast (9796),
    UScast (9947), UScast (10098), UScast (10249),
    UScast (10400), UScast (10551), UScast (10702),
    UScast (10853), UScast (11004), UScast (11155),
    UScast (11306), UScast (11457), UScast (11608),
    UScast (11759), UScast (11910), UScast (12061),
    UScast (12212), UScast (12363), UScast (12514),
    UScast (12665), UScast (12816), UScast (12967),
    UScast (13118), UScast (13269), UScast (13420),
    UScast (13571), UScast (13722), UScast (13873),
    UScast (14024), UScast (14175), UScast (14251),
    UScast (14327), UScast (14403), UScast (14479),
    UScast (14555), UScast (14631), UScast (14707),
    UScast (14783), UScast (14859), UScast (14935),
    UScast (15011), UScast (15087), UScast (15163),
    UScast (15239), UScast (15315), UScast (15391),
    UScast (15467), UScast (15543), UScast (15619),
    UScast (15695), UScast (15771), UScast (15847),
    UScast (15923), UScast (15999), UScast (16075),
    UScast (16151), UScast (16227), UScast (16303),
    UScast (16379), UScast (16455), UScast (16531),
    UScast (16607), UScast (16683), UScast (16759),
    UScast (16835), UScast (16911), UScast (16987),
    UScast (17063), UScast (17064), UScast (17065),
    UScast (17066), UScast (17067), UScast (17068),
    UScast (17069), UScast (17070), UScast (17071),
    UScast (17072), UScast (17073), UScast (17074),
    UScast (17075), UScast (17076), UScast (17077),
    UScast (17078), UScast (17079), UScast (17080),
    UScast (17081), UScast (17082), UScast (17083),
    UScast (17084), UScast (17085), UScast (17086),
    UScast (17087), UScast (17088), UScast (17089),
    UScast (17090), UScast (17091), UScast (17092),
    UScast (17093), UScast (17094), UScast (17095),
    UScast (17096), UScast (17097), UScast (17098),
    UScast (17099), UScast (17100)
};

static Diffrank differencecoversizes[] = {
    UCcast (1), /* 2^0 */
    UCcast (2), /* 2^1 */
    UCcast (3), /* 2^2 */
    UCcast (4), /* 2^3 */
    UCcast (5), /* 2^4 */
    UCcast (7), /* 2^5 */
    UCcast (9), /* 2^6 */
    UCcast (13), /* 2^7 */
    UCcast (20), /* 2^8 */
    UCcast (28), /* 2^9 */
    UCcast (40), /* 2^10 */
    UCcast (58), /* 2^11 */
    UCcast (82), /* 2^12 */
    UCcast (112), /* 2^13 */
    UCcast (160), /* 2^14 */
    UCcast (226) /* 2^15 */
};

#endif
