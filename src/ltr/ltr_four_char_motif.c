/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "ltr/ltr_four_char_motif.h"

/* test the motif and encode the characters by using alpha */
int gt_ltr_four_char_motif_encode (GtLTRFourCharMotif *motif,
                                   const GtEncseq *encseq,
                                   GtError *err)
{
  const GtUchar *symbolmap;
  GtUchar c_tab[UCHAR_MAX+1];
  unsigned int i;

  symbolmap = gt_alphabet_symbolmap(gt_encseq_alphabet(encseq));
  if ( symbolmap[(unsigned int)motif->firstleft] == (GtUchar) UNDEFCHAR)
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondleft] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondleft);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->firstright] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->firstright);
    return -1;
  }
  if ( symbolmap[(unsigned int)motif->secondright] == (GtUchar) UNDEFCHAR )
  {
    gt_error_set(err,"Illegal nucleotide character %c "
                      "as argument to option -motif", motif->secondright);
    return -1;
  }

  for (i=0; i<=(unsigned int) UCHAR_MAX; i++)
  {
    c_tab[i] = (GtUchar) UNDEFCHAR;
  }
  /* define complementary symbols */
  c_tab[symbolmap['a']] = symbolmap['t'];
  c_tab[symbolmap['c']] = symbolmap['g'];
  c_tab[symbolmap['g']] = symbolmap['c'];
  c_tab[symbolmap['t']] = symbolmap['a'];

  /* if motif is not palindromic */
  if ( (c_tab[symbolmap[(unsigned int)motif->firstleft]] !=
       c_tab[c_tab[symbolmap[(unsigned int)motif->secondright]]])
           ||
      (c_tab[symbolmap[(unsigned int)motif->secondleft]] !=
       c_tab[c_tab[symbolmap[(unsigned int)motif->firstright]]]) )
  {
    gt_error_set(err, "Illegal motif, motif not palindromic");
    return -1;
  }

  /* encode the symbols */
  motif->firstleft = symbolmap[(unsigned int)motif->firstleft];
  motif->secondleft = symbolmap[(unsigned int)motif->secondleft];
  motif->firstright = symbolmap[(unsigned int)motif->firstright];
  motif->secondright = symbolmap[(unsigned int)motif->secondright];

  return 0;
}
