/* $Id: rx_cls_u.c,v 1.2 2003/12/21 18:57:34 dvd Exp $ */

#define CLS_U_  0
#define CLS_U_C  1
#define CLS_U_Cc  2
#define CLS_U_Cf  3
#define CLS_U_Co  4
#define CLS_U_IsAlphabeticPresentationForms  5
#define CLS_U_IsArabic  6
#define CLS_U_IsArabicPresentationForms_A  7
#define CLS_U_IsArabicPresentationForms_B  8
#define CLS_U_IsArmenian  9
#define CLS_U_IsArrows  10
#define CLS_U_IsBasicLatin  11
#define CLS_U_IsBengali  12
#define CLS_U_IsBlockElements  13
#define CLS_U_IsBopomofo  14
#define CLS_U_IsBopomofoExtended  15
#define CLS_U_IsBoxDrawing  16
#define CLS_U_IsBraillePatterns  17
#define CLS_U_IsByzantineMusicalSymbols  18
#define CLS_U_IsCJKCompatibility  19
#define CLS_U_IsCJKCompatibilityForms  20
#define CLS_U_IsCJKCompatibilityIdeographs  21
#define CLS_U_IsCJKCompatibilityIdeographsSupplement  22
#define CLS_U_IsCJKRadicalsSupplement  23
#define CLS_U_IsCJKSymbolsandPunctuation  24
#define CLS_U_IsCJKUnifiedIdeographs  25
#define CLS_U_IsCJKUnifiedIdeographsExtensionA  26
#define CLS_U_IsCJKUnifiedIdeographsExtensionB  27
#define CLS_U_IsCherokee  28
#define CLS_U_IsCombiningDiacriticalMarks  29
#define CLS_U_IsCombiningHalfMarks  30
#define CLS_U_IsCombiningMarksforSymbols  31
#define CLS_U_IsControlPictures  32
#define CLS_U_IsCurrencySymbols  33
#define CLS_U_IsCyrillic  34
#define CLS_U_IsDeseret  35
#define CLS_U_IsDevanagari  36
#define CLS_U_IsDingbats  37
#define CLS_U_IsEnclosedAlphanumerics  38
#define CLS_U_IsEnclosedCJKLettersandMonths  39
#define CLS_U_IsEthiopic  40
#define CLS_U_IsGeneralPunctuation  41
#define CLS_U_IsGeometricShapes  42
#define CLS_U_IsGeorgian  43
#define CLS_U_IsGothic  44
#define CLS_U_IsGreek  45
#define CLS_U_IsGreekExtended  46
#define CLS_U_IsGujarati  47
#define CLS_U_IsGurmukhi  48
#define CLS_U_IsHalfwidthandFullwidthForms  49
#define CLS_U_IsHangulCompatibilityJamo  50
#define CLS_U_IsHangulJamo  51
#define CLS_U_IsHangulSyllables  52
#define CLS_U_IsHebrew  53
#define CLS_U_IsHiragana  54
#define CLS_U_IsIPAExtensions  55
#define CLS_U_IsIdeographicDescriptionCharacters  56
#define CLS_U_IsKanbun  57
#define CLS_U_IsKangxiRadicals  58
#define CLS_U_IsKannada  59
#define CLS_U_IsKatakana  60
#define CLS_U_IsKhmer  61
#define CLS_U_IsLao  62
#define CLS_U_IsLatin_1Supplement  63
#define CLS_U_IsLatinExtended_A  64
#define CLS_U_IsLatinExtended_B  65
#define CLS_U_IsLatinExtendedAdditional  66
#define CLS_U_IsLetterlikeSymbols  67
#define CLS_U_IsMalayalam  68
#define CLS_U_IsMathematicalAlphanumericSymbols  69
#define CLS_U_IsMathematicalOperators  70
#define CLS_U_IsMiscellaneousSymbols  71
#define CLS_U_IsMiscellaneousTechnical  72
#define CLS_U_IsMongolian  73
#define CLS_U_IsMusicalSymbols  74
#define CLS_U_IsMyanmar  75
#define CLS_U_IsNumberForms  76
#define CLS_U_IsOgham  77
#define CLS_U_IsOldItalic  78
#define CLS_U_IsOpticalCharacterRecognition  79
#define CLS_U_IsOriya  80
#define CLS_U_IsPrivateUse  81
#define CLS_U_IsRunic  82
#define CLS_U_IsSinhala  83
#define CLS_U_IsSmallFormVariants  84
#define CLS_U_IsSpacingModifierLetters  85
#define CLS_U_IsSpecials  86
#define CLS_U_IsSuperscriptsandSubscripts  87
#define CLS_U_IsSyriac  88
#define CLS_U_IsTags  89
#define CLS_U_IsTamil  90
#define CLS_U_IsTelugu  91
#define CLS_U_IsThaana  92
#define CLS_U_IsThai  93
#define CLS_U_IsTibetan  94
#define CLS_U_IsUnifiedCanadianAboriginalSyllabics  95
#define CLS_U_IsYiRadicals  96
#define CLS_U_IsYiSyllables  97
#define CLS_U_L  98
#define CLS_U_Ll  99
#define CLS_U_Lm  100
#define CLS_U_Lo  101
#define CLS_U_Lt  102
#define CLS_U_Lu  103
#define CLS_U_M  104
#define CLS_U_Mc  105
#define CLS_U_Me  106
#define CLS_U_Mn  107
#define CLS_U_N  108
#define CLS_U_Nd  109
#define CLS_U_Nl  110
#define CLS_U_No  111
#define CLS_U_P  112
#define CLS_U_Pc  113
#define CLS_U_Pd  114
#define CLS_U_Pe  115
#define CLS_U_Pf  116
#define CLS_U_Pi  117
#define CLS_U_Po  118
#define CLS_U_Ps  119
#define CLS_U_S  120
#define CLS_U_Sc  121
#define CLS_U_Sk  122
#define CLS_U_Sm  123
#define CLS_U_So  124
#define CLS_U_Z  125
#define CLS_U_Zl  126
#define CLS_U_Zp  127
#define CLS_U_Zs  128
#define NUM_CLS_U 129
static char *clstab[NUM_CLS_U]={"",
  "C",
  "Cc",
  "Cf",
  "Co",
  "IsAlphabeticPresentationForms",
  "IsArabic",
  "IsArabicPresentationForms-A",
  "IsArabicPresentationForms-B",
  "IsArmenian",
  "IsArrows",
  "IsBasicLatin",
  "IsBengali",
  "IsBlockElements",
  "IsBopomofo",
  "IsBopomofoExtended",
  "IsBoxDrawing",
  "IsBraillePatterns",
  "IsByzantineMusicalSymbols",
  "IsCJKCompatibility",
  "IsCJKCompatibilityForms",
  "IsCJKCompatibilityIdeographs",
  "IsCJKCompatibilityIdeographsSupplement",
  "IsCJKRadicalsSupplement",
  "IsCJKSymbolsandPunctuation",
  "IsCJKUnifiedIdeographs",
  "IsCJKUnifiedIdeographsExtensionA",
  "IsCJKUnifiedIdeographsExtensionB",
  "IsCherokee",
  "IsCombiningDiacriticalMarks",
  "IsCombiningHalfMarks",
  "IsCombiningMarksforSymbols",
  "IsControlPictures",
  "IsCurrencySymbols",
  "IsCyrillic",
  "IsDeseret",
  "IsDevanagari",
  "IsDingbats",
  "IsEnclosedAlphanumerics",
  "IsEnclosedCJKLettersandMonths",
  "IsEthiopic",
  "IsGeneralPunctuation",
  "IsGeometricShapes",
  "IsGeorgian",
  "IsGothic",
  "IsGreek",
  "IsGreekExtended",
  "IsGujarati",
  "IsGurmukhi",
  "IsHalfwidthandFullwidthForms",
  "IsHangulCompatibilityJamo",
  "IsHangulJamo",
  "IsHangulSyllables",
  "IsHebrew",
  "IsHiragana",
  "IsIPAExtensions",
  "IsIdeographicDescriptionCharacters",
  "IsKanbun",
  "IsKangxiRadicals",
  "IsKannada",
  "IsKatakana",
  "IsKhmer",
  "IsLao",
  "IsLatin-1Supplement",
  "IsLatinExtended-A",
  "IsLatinExtended-B",
  "IsLatinExtendedAdditional",
  "IsLetterlikeSymbols",
  "IsMalayalam",
  "IsMathematicalAlphanumericSymbols",
  "IsMathematicalOperators",
  "IsMiscellaneousSymbols",
  "IsMiscellaneousTechnical",
  "IsMongolian",
  "IsMusicalSymbols",
  "IsMyanmar",
  "IsNumberForms",
  "IsOgham",
  "IsOldItalic",
  "IsOpticalCharacterRecognition",
  "IsOriya",
  "IsPrivateUse",
  "IsRunic",
  "IsSinhala",
  "IsSmallFormVariants",
  "IsSpacingModifierLetters",
  "IsSpecials",
  "IsSuperscriptsandSubscripts",
  "IsSyriac",
  "IsTags",
  "IsTamil",
  "IsTelugu",
  "IsThaana",
  "IsThai",
  "IsTibetan",
  "IsUnifiedCanadianAboriginalSyllabics",
  "IsYiRadicals",
  "IsYiSyllables",
  "L",
  "Ll",
  "Lm",
  "Lo",
  "Lt",
  "Lu",
  "M",
  "Mc",
  "Me",
  "Mn",
  "N",
  "Nd",
  "Nl",
  "No",
  "P",
  "Pc",
  "Pd",
  "Pe",
  "Pf",
  "Pi",
  "Po",
  "Ps",
  "S",
  "Sc",
  "Sk",
  "Sm",
  "So",
  "Z",
  "Zl",
  "Zp",
  "Zs"
  };
