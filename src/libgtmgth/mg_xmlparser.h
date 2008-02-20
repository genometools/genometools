
/*
  Copyright (c) 2007 David Schmitz-Huebsch <dschmitz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <expat.h>
#ifndef CURLDEF
#include <metagenomethreader.h>
#endif
#ifdef CURLDEF
#include <curl.h>
#endif

/* Funktion, in der die Struktur 'gefuellt' und der XML-Parser
   initialisiert und gestartet wird
   Parameter:  Struktur, die als Data den Expat-Funktionen uebergeben wird;
               XML-Filename; Error-Variable
   Returnwert: Fehlercode (int) */
int mg_xmlparser(ParseStruct *, GenFile *, Error *);

/* Expat-Funktion zur Verarbeitung des Textes zwischen XML-Tags
   Parameter: void-Zeiger auf die (Nutzer)Daten - hier: parsestruct-Struktur;
              Zeiger auf das erste Zeichen nach einem XML-Tag;
              int-Wert der Laenge des Textes zwischen den XML-Tags
   Returnwert: void */
void text(void *, const XML_Char *, int);
