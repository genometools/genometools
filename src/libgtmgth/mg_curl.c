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

#ifdef CURLDEF

#include "metagenomethreader.h"

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

/* Expat-Hilfsfunktion zum Abspeichern der empfangenen Daten
   Parameter:  void-Zeiger auf den Anfang des Speicherbereichs,
               Anzahl zu speichernder Elemente, Groesse des zu
               speichernden Datentyps, void-Zeiger auf die zu
               speichernden Daten
   Returnwert: Groesse des neu allokierten Speicherbereichs*/
size_t WriteMemoryCallback(void *, size_t, size_t, void *);

size_t WriteMemoryCallback(void *ptr,
                           size_t size, size_t nmemb, void *data)
{
  size_t realsize = size * nmemb;
  MemoryStruct *mem = (MemoryStruct *) data;

  mem->memory = (char *) realloc(mem->memory, mem->size + realsize + 1);
  if (mem->memory)
  {
    memcpy(&(mem->memory[mem->size]), ptr, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;
  }
  return realsize;
}

int mg_curl(ParseStruct *parsestruct_ptr,
            unsigned long hit_counter, Error * err)
{
  int had_err = 0,
    curl_errornr = 0;

  /* Laenge der aus dem XML-File stammenden Hit-DNA-Sequnez */
  unsigned long seq_len;
  long numb_from = 0, numb_to = 0, numb_diff = 0;

  Str *seq_var,
   *http_adr;

  MemoryStruct memorystruct;

  /* char-Zeiger auf die HTTP-Adresse des cgi-Skriptes efetch von NCBI */
  char *http_adr_ptr,
   *seq_pos;                           /* char-Zeiger, wird benutzt zum
                                          Auslesen der Sequenzinformation
                                          aus dem XML-File, welche
                                          Ergebnis der efetch-Anfrage ist */
  const char *curlerror;

  /* Curl-Handle */
  CURL *curl_handle;

  /* char-Zeiger auf die Daten ist NULL */
  memorystruct.memory = NULL;
  /* noch keine Daten eingetragen bzw. abgespeichert */
  memorystruct.size = 0;

  /* Zwischenspeicher fuer die Sequnezinformation, da die StrArray-Klasse
     keine Funktion zum begrenzten Einfuegen eines Strings zur Verfuegung
     stellt; setzen des ersten Teils der HTTP-Adresse */
  seq_var = str_new();
  http_adr =
    str_new_cstr
    ("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=");

  /* Check der Umgebungsvariablen */
  error_check(err);

  curl_global_init(CURL_GLOBAL_ALL);

  /* initialisieren der curl-session */
  curl_handle = curl_easy_init();

  /* Zusammensetzen der http-Adresse durch Anhaengen der query-GI-Nummer,
     des Hit-from, des Hit-to Wertes und des Rueckgabetyps an den ersten
     Teil der HTTP-Adresse */
  str_append_str(http_adr, ARGUMENTSSTRUCT(curl_fcgi_db));
  str_append_cstr(http_adr, "&id=gi|");
  str_append_str(http_adr, parsestruct_ptr->hit_gi_nr_tmp);
  str_append_cstr(http_adr, "&seq_start=");
  str_append_cstr(http_adr,
                  strarray_get(MATRIXSTRUCT(hit_from), hit_counter));
  str_append_cstr(http_adr, "&seq_stop=");
  str_append_cstr(http_adr,
                  strarray_get(MATRIXSTRUCT(hit_to), hit_counter));
  str_append_cstr(http_adr, "&retmode=xml");

  /* char-Zeiger wird benoetigt, da curl_easy_setopt als 3. Parameter
     einen char-Zeiger erwartet */
  http_adr_ptr = str_get(http_adr);

  /* festlegen, welche HTTP-Adresse aufgerufen werden soll */
  curl_easy_setopt(curl_handle, CURLOPT_URL, http_adr_ptr);

  /* die empfangenen Daten werden an die Funktion WriteMemoryCallback
     gesendet, wo Speicherplatz reserviert und die Daten in diesen
     Speicherbereich kopiert werden */
  curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION,
                   WriteMemoryCallback);

  /* Die Daten werden in die Struktur eingetragen */
  curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *) &memorystruct);

  /* setzen des user-agent field, da einige Server diesen voraussetzen */
  curl_easy_setopt(curl_handle, CURLOPT_USERAGENT, "libcurl-agent/1.0");

  /* Anfrage wird ausgefuehrt */
  curl_errornr = curl_easy_perform(curl_handle);
  curlerror = curl_easy_strerror(curl_errornr);

  if (curl_errornr)
  {
    error_set(err,
              "an error occurred during curl-processing (error-code %d):\
               \"%s\"", curl_errornr, curlerror);
    had_err = -1;
  }

  if (!had_err)
  {
    /* Die Hit-DNA steht zwischen dem <GBSeq_sequence> und dem
       </GBSeq_sequence> XML-Tag, Zeiger auf das < Zeichen von
       <GBSeq_sequence> */
    seq_pos = strstr(memorystruct.memory, "<GBSeq_sequence>");

    if (!seq_pos)
    {
      error_set(err,
                "an error occurred while retrieving sequence-information\
                 with the following request: \"%s\"", http_adr_ptr);
      had_err = -1;
    }

    if (!had_err)
    {
      /* seq_pos+16 zeigt auf das erste Zeichen der Sequence; gezaehlt
         wird die Laenge bis zum naechsten Zeichen das kein g,a,c oder t
         etc. ist */
      assert(seq_pos != NULL);
      seq_len = strspn(seq_pos + 16, "gactrymkswhbvdnu");

      numb_from = atol(strarray_get(MATRIXSTRUCT(hit_from), hit_counter));
      numb_to = atol(strarray_get(MATRIXSTRUCT(hit_to), hit_counter));

      numb_diff = numb_to - numb_from +1;

      if (numb_diff == seq_len)
      {
        /* seq_len Zeichen werden in die Hilfsvariable seq_var kopiert */
        str_append_cstr_nt(seq_var, seq_pos + 16, seq_len);

        /* Die Sequenz in seq_var wird in das StrArray hit_dna kopiert */
        strarray_add_cstr(MATRIXSTRUCT(hit_dna), str_get(seq_var));

        /* das Hit-Sequenz-File wird geschrieben; die erste Zeile eines
           Eintrages ist die Hit-GI-Def, an die durch ein Leerzeichen
           getrennt die Hsp-Num des jeweiligen Hits angehaengt wird */
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, ">%s ",
                        strarray_get(MATRIXSTRUCT(hit_num), hit_counter));
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                        strarray_get(MATRIXSTRUCT(hit_from), hit_counter));
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                        strarray_get(MATRIXSTRUCT(hit_to), hit_counter));
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                        strarray_get(MATRIXSTRUCT(fasta_row), hit_counter));
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, "%s\n",
                        strarray_get(MATRIXSTRUCT(hit_gi_def), hit_counter));

        /* nach dem GI-Def Eintrag folgt in der naechsten Zeile die Sequenz */
        genfile_xprintf(parsestruct_ptr->fp_blasthit_file, "%s\n",
                        str_get(seq_var));
      }
      else
      {
        PARSESTRUCT(gi_flag) = 1;
      }
    }
  }
  /* cleanup des Curl-Handle */
  curl_easy_cleanup(curl_handle);

  if (memorystruct.memory)
    free(memorystruct.memory);

  str_delete(http_adr);
  str_delete(seq_var);

  return had_err;
}

#endif
