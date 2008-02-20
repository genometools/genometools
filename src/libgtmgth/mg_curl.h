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

#ifndef MG_CURL_H
#define MG_CURL_H

#include "libgtmgth/metagenomethreader.h"

#ifndef S_SPLINT_S
#include <curl/curl.h>
#include <curl/types.h>
#include <curl/easy.h>
#endif

/* Expat Hilfs-Struktur zum Abspeichern der Hit-DNA-Sequenz und der
   Hit-DNA-Laenge */
typedef struct
{
  char *memory;
  size_t size;
} MemoryStruct;

/* Expat-Hilfsfunktion zur dyn. Speicherreservierung
   Parameter:  void-Zeiger auf den Anfang des Speicherbereichs,
               Groesse des momentan reservierten Speicherbereichs
   Returnwert: void */
void *myrealloc(void *, size_t);

/* Expat-Hilfsfunktion zum Abspeichern der empfangenen Daten
   Parameter:  void-Zeiger auf den Anfang des Speicherbereichs,
               Anzahl zu speichernder Elemente, Groesse des zu
               speichernden Datentyps, void-Zeiger auf die zu
               speichernden Daten
   Returnwert: Groesse des neu allokierten Speicherbereichs*/
size_t WriteMemoryCallback(void *, size_t, size_t, void *);

#endif
