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

#include <ctype.h>
#include <expat.h>
#include "mg_xmlparser.h"
#include "core/unused_api.h"
#include "metagenomethreader.h"

#ifdef CURLDEF
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

static size_t WriteMemoryCallback(void *ptr,
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
            unsigned long hit_counter, GtError * err)
{
  int had_err = 0,
    curl_errornr = 0;

  /* Laenge der aus dem XML-File stammenden Hit-DNA-Sequenz */
  unsigned long seq_len;
  long numb_from = 0, numb_to = 0, numb_diff = 0;

  GtStr *seq_var,
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

  /* Zwischenspeicher fuer die Sequnezinformation, da die GtStrArray-Klasse
     keine Funktion zum begrenzten Einfuegen eines Strings zur Verfuegung
     stellt; setzen des ersten Teils der HTTP-Adresse */
  seq_var = gt_str_new();
  http_adr =
    gt_str_new_cstr
    ("http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=");

  /* Check der Umgebungsvariablen */
  gt_error_check(err);

  curl_global_init(CURL_GLOBAL_ALL);

  /* initialisieren der curl-session */
  curl_handle = curl_easy_init();

  /* Zusammensetzen der http-Adresse durch Anhaengen der query-GI-Nummer,
     des Hit-from, des Hit-to Wertes und des Rueckgabetyps an den ersten
     Teil der HTTP-Adresse */
  gt_str_append_str(http_adr, ARGUMENTSSTRUCT(curl_fcgi_db));
  gt_str_append_cstr(http_adr, "&id=gi|");
  gt_str_append_str(http_adr, parsestruct_ptr->hit_gi_nr_tmp);
  gt_str_append_cstr(http_adr, "&seq_start=");
  gt_str_append_cstr(http_adr,
                  gt_str_array_get(MATRIXSTRUCT(hit_from), hit_counter));
  gt_str_append_cstr(http_adr, "&seq_stop=");
  gt_str_append_cstr(http_adr,
                  gt_str_array_get(MATRIXSTRUCT(hit_to), hit_counter));
  gt_str_append_cstr(http_adr, "&retmode=xml");

  /* char-Zeiger wird benoetigt, da curl_easy_setopt als 3. Parameter
     einen char-Zeiger erwartet */
  http_adr_ptr = gt_str_get(http_adr);

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
    gt_error_set(err,
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
      gt_error_set(err,
                "an error occurred while retrieving sequence-information\
                 with the following request: \"%s\"", http_adr_ptr);
      had_err = -1;
    }

    if (!had_err)
    {
      /* seq_pos+16 zeigt auf das erste Zeichen der Sequence; gezaehlt
         wird die Laenge bis zum naechsten Zeichen das kein g,a,c oder t
         etc. ist */
      gt_assert(seq_pos != NULL);
      seq_len = strspn(seq_pos + 16, "gactrymkswhbvdnu");

      numb_from = atol(gt_str_array_get(MATRIXSTRUCT(hit_from), hit_counter));
      numb_to = atol(gt_str_array_get(MATRIXSTRUCT(hit_to), hit_counter));

      numb_diff = numb_to - numb_from +1;

      if (numb_diff == seq_len)
      {
        /* seq_len Zeichen werden in die Hilfsvariable seq_var kopiert */
        gt_str_append_cstr_nt(seq_var, seq_pos + 16, seq_len);

        /* Die Sequenz in seq_var wird in das GtStrArray hit_dna kopiert */
        gt_str_array_add_cstr(MATRIXSTRUCT(hit_dna), gt_str_get(seq_var));

        /* das Hit-Sequenz-File wird geschrieben; die erste Zeile eines
           Eintrages ist die Hit-GI-Def, an die durch ein Leerzeichen
           getrennt die Hsp-Num des jeweiligen Hits angehaengt wird */
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, ">%s ",
                           gt_str_array_get(MATRIXSTRUCT(hit_num),
                                            hit_counter));
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                           gt_str_array_get(MATRIXSTRUCT(hit_from),
                                            hit_counter));
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                           gt_str_array_get(MATRIXSTRUCT(hit_to),
                                            hit_counter));
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, "%s ",
                           gt_str_array_get(MATRIXSTRUCT(fasta_row),
                                            hit_counter));
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, "%s\n",
                           gt_str_array_get(MATRIXSTRUCT(hit_gi_def),
                                            hit_counter));

        /* nach dem GI-Def Eintrag folgt in der naechsten Zeile die Sequenz */
        gt_file_xprintf(parsestruct_ptr->fp_blasthit_file, "%s\n",
                        gt_str_get(seq_var));
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

  gt_str_delete(http_adr);
  gt_str_delete(seq_var);

  return had_err;
}

#endif

/* Expat-Funktion zur Verarbeitung des Textes zwischen XML-Tags
   Parameter: void-Zeiger auf die (Nutzer)Daten - hier: parsestruct-Struktur;
              Zeiger auf das erste Zeichen nach einem XML-Tag;
              int-Wert der Laenge des Textes zwischen den XML-Tags
   Returnwert: void */
static void textElement(void *, const XML_Char *, int);

/* Expat-Funktion zur Behandlung oeffnender XML-Tags
   Parameter: void-Zeiger auf die (Nutzer)Daten - hier: parsestruct-Struktur;
              Name des aktuellen XML-Tags;
              Expat XML-Attribute (nicht verwendet)
   Returnwert: void */
static void XMLCALL startElement(void *, const char *, const char **);

/* Expat-Funktion zur Behandlung schliessender XML-Tags
   Parameter: void-Zeiger auf die (Nutzer)Daten - hier: parsestruct-Struktur;
              Name des aktuellen XML-Tags
   Returnwert: void */
static void XMLCALL endElement(void *, const char *);

/* Funktion zum Setzen von Flags, bei Uebereinstimmung von gesuchtem und
   aktuell oeffnendem XML-Tag-Namen
   Parameter: das zu setzende Flag; der Counter, der das zu betrachtende
              XML-Tag bestimmt; ein Zeichen zur Unterscheidung von Query-,
              Hit- und Hsp-XML-Tags
   Returnwert: void */
static void flag_setting(unsigned short *, unsigned short *, char);

/* Funktion zum Loeschen von Flags, bei Uebereinstimmung von gesuchtem und
   aktuell schliessendem XML-Tag-Namen
   Parameter: das zu setzende Flag; der counter, der das zu betrachtende
              XML-Tag bestimmt; ein Zeichen zur Unterscheidung von Query-,
              Hit- und Hsp-XML-Tags
   Returnwert: void */
static void flag_delete(unsigned short *, unsigned short *, char);

/* Funktion, die den Stand der Counter ueberprueft und bei Bedarf auf
   UNSET setzt; stimmt der Counter mit den definierten Grenzen ueberein
   wird er auf UNSET gesetzt
   Parameter: der counter, der das zu betrachtende XML-Tag bestimmt; ein
              Zeichen zur Unterscheidung von Query-, Hit- und Hsp-XML-Tags
   Returnwert: void */
static void check_counter(unsigned short *, char);

int mg_xmlparser(ParseStruct *parsestruct_ptr, GtFile * fp_xmlfile,
                 GtError * err)
{
  int had_err = 0;

  /* Expat XML-Error setzen */
  enum XML_Error error;

  /* Puffer zum zeilenweisen Einlesen des XML-Files, XML_Parser
     deklarieren */
  GtStr *buf;
  XML_Parser parser;

  /* Check Umgebungsvariablen */
  gt_error_check(err);

  /* Puffer-String anlegen */
  buf = gt_str_new();

  /* XML-Parser wird initialisiert */
  parser = XML_ParserCreate(NULL);
  /* Die Struktur parsestruct wird als UserData gesetzt */
  XML_SetUserData(parser, parsestruct_ptr);
  /* Start- und Endelement Handler werden gesetzt */
  XML_SetElementHandler(parser, startElement, endElement);
  /* Text-Handler setzen */
  XML_SetCharacterDataHandler(parser, textElement);

  while ((gt_str_read_next_line_generic(buf, fp_xmlfile) != EOF) && !had_err)
  {
    PARSESTRUCT(xml_linenumber)++;
    if (PARSESTRUCT(had_err))
    {
      had_err = -1;
    }

    if ((XML_Parse(parser, gt_str_get(buf), gt_str_length(buf), false) ==
         XML_STATUS_ERROR) && !had_err)
    {
      error = XML_GetErrorCode(parser);
      gt_error_set(err,
                "an error occurred parsing line %lu of file \"%s\": %s",
                PARSESTRUCT(xml_linenumber), gt_str_get(PARSESTRUCT(xmlfile)),
                XML_ErrorString(error));

      had_err = -1;
    }
    /* reset line buffer */
    gt_str_reset(buf);
  }

  /* finish parsing */
  if ((XML_Parse(parser, NULL, 0, true) == XML_STATUS_ERROR) && !had_err)
  {
    error = XML_GetErrorCode(parser);
    gt_error_set(err,
              "an error occurred while finishing the parsing of file\
               \"%s\": %s",
              gt_str_get(PARSESTRUCT(xmlfile)), XML_ErrorString(error));

    had_err = -1;
  }

  if (PARSESTRUCT(xml_tag_flag) && !(!had_err) && PARSESTRUCT(giexp_flag))
  {
    gt_str_array_delete(MATRIXSTRUCT(hit_gi_nr));
    gt_str_array_delete(MATRIXSTRUCT(hit_num));
    gt_str_array_delete(MATRIXSTRUCT(hit_dna));
    gt_str_array_delete(MATRIXSTRUCT(hit_gi_def));
    gt_str_array_delete(MATRIXSTRUCT(hit_acc));
    gt_str_array_delete(MATRIXSTRUCT(fasta_row));
    gt_str_array_delete(MATRIXSTRUCT(hit_from));
    gt_str_array_delete(MATRIXSTRUCT(hit_to));
    gt_str_array_delete(MATRIXSTRUCT(hsp_qseq));
    gt_str_array_delete(MATRIXSTRUCT(hsp_hseq));
    gt_str_array_delete(PARSESTRUCT(query_frame_tmp));
    gt_str_array_delete(PARSESTRUCT(hit_frame_tmp));
  }

  /* Freigeben des XML-Parser und loeschen des Puffer-Strings */
  XML_ParserFree(parser);
  gt_str_delete(buf);

  return had_err;
}

static void XMLCALL startElement(void *data, const char *name,
                                 GT_UNUSED const char **atts)
{
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  if (!PARSESTRUCT(had_err))
  {
    /* Vergleich des aktuellen XML-Tags mit dem aktuellen Tag-Namen im
       jew. GtArray bei Uebereinstimmung wird das jeweilige Flag gesetzt */
    if (strcmp
        (name,
         gt_str_array_get(PARSESTRUCT(query_array),
                      XMLPARSERSTRUCT(query_array_index_start))) == 0)
    {
      flag_setting(&XMLPARSERSTRUCT(query_array_index_start),
                   &PARSESTRUCT(def_flag), 'q');
    }
    else
      if (strcmp
          (name,
           gt_str_array_get(PARSESTRUCT(hit_array),
                        XMLPARSERSTRUCT(hit_array_index_start))) == 0)
    {
      flag_setting(&XMLPARSERSTRUCT(hit_array_index_start),
                   &PARSESTRUCT(hit_flag), 'h');
    }
    else
      if (strcmp
          (name,
           gt_str_array_get(PARSESTRUCT(hit_hsp_array),
                        XMLPARSERSTRUCT(hit_hsp_array_index_start))) == 0)
    {
      flag_setting(&XMLPARSERSTRUCT(hit_hsp_array_index_start),
                   &PARSESTRUCT(hit_hsp_flag), 't');
    }
  }
}

static void XMLCALL endElement(void *data, const char *name)
{
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  if (!PARSESTRUCT(had_err))
  {
    GtError *err = PARSESTRUCT(err);

    /* Laenge der GI-Nr in einem GI-Def-XML-Eintrag */
    unsigned short gi_len = 0;

    /* Temp-Variablen zum Zwischenspeichern der Query-Start- bzw.
       -End-Werte sowie der Frame-Informationen */
    unsigned long ulong_numb_buf = 0,
      query_nr = 0, **query_nr_p;
    long numb_buf = 0;

    /* Zeiger auf die erste Zahl der GI-Nr in einem Hit-ID-XML-Eintrag */
    const char *gi_ptr = NULL;

    /* Check Umgebungsvariablen */
    gt_error_check(err);

    /* wurde das Iteration_hits XML-Tag erreicht wird die mg_combinedscore
       Methode aufgerufen und diverse Variablen fuer einen neuen Eintrag
       resetet bzw. geloescht */
    if (strcmp(name, gt_str_get(PARSESTRUCT(xml_tag))) == 0
                             && PARSESTRUCT(giexp_flag))
    {
      if (XMLPARSERSTRUCT(hit_counter) > 0)
      {
        PARSESTRUCT(had_err) =
          mg_combinedscore(parsestruct_ptr, XMLPARSERSTRUCT(hit_counter),
                           err);

        /* Zaehler der Hits pro Hit-GI-Nr */
        XMLPARSERSTRUCT(hit_counter) = 0;
        PARSESTRUCT(gi_flag) = 0;

        gt_array_reset(MATRIXSTRUCT(query_from));
        gt_array_reset(MATRIXSTRUCT(query_to));
        gt_array_reset(MATRIXSTRUCT(hit_frame));
        gt_array_reset(MATRIXSTRUCT(query_frame));
      }

      /* Schreiben der schliessenden XML-Tags nach der Hit-Bearbeitung */
      if (ARGUMENTSSTRUCT(outputfile_format) == 3)
      {
        /* Abschluss des Iteration-Bereichs im XML-File */
        mg_outputwriter(parsestruct_ptr, NULL, NULL, NULL, 'x', err);
      }

      gt_str_array_delete(MATRIXSTRUCT(hit_gi_nr));
      gt_str_array_delete(MATRIXSTRUCT(hit_gi_def));
      gt_str_array_delete(MATRIXSTRUCT(hit_acc));
      gt_str_array_delete(MATRIXSTRUCT(fasta_row));
      gt_str_array_delete(MATRIXSTRUCT(hit_num));

      gt_str_array_delete(MATRIXSTRUCT(hit_dna));
      gt_str_array_delete(MATRIXSTRUCT(hit_from));
      gt_str_array_delete(MATRIXSTRUCT(hit_to));
      gt_str_array_delete(MATRIXSTRUCT(hsp_qseq));
      gt_str_array_delete(MATRIXSTRUCT(hsp_hseq));
      gt_str_array_delete(PARSESTRUCT(query_frame_tmp));
      gt_str_array_delete(PARSESTRUCT(hit_frame_tmp));

      PARSESTRUCT(xml_tag_flag) = MGTH_UNSET;
    }
    /* nur wenn ein Flag gesetzt ist und Text eingelesen wurde, ist eine
       weitere Bearbeitung notwendig */
    if ((PARSESTRUCT(def_flag) == 1 || PARSESTRUCT(hit_flag) == 1
         || PARSESTRUCT(hit_hsp_flag) == 1))
    {
      /* Anhand der Query-GI-Def kann mittels Hash-Table und
         GtBioseq-Struktur die Query-DNA ausgelesen werden */
      if (strcmp(name, gt_str_array_get(PARSESTRUCT(query_array), 0)) == 0
                                    && PARSESTRUCT(giexp_flag))
      {
        /* Query-DNA-Strings fuer die Sequenz und die Definition werden
           geloescht, bevor der neue Eintrag hineinkopiert wird */
        gt_str_reset(MATRIXSTRUCT(query_dna));
        gt_str_reset(MATRIXSTRUCT(query_def));

        /* Abspeichern der Query-Def */
        gt_str_set(MATRIXSTRUCT(query_def), gt_str_get(PARSESTRUCT(buf_ptr)));

        /* Erzeugen von StringArrays zur Aufnahme der
           Hit-DNA-Sequenz-Informationen */
        MATRIXSTRUCT(hit_gi_nr) = gt_str_array_new();
        MATRIXSTRUCT(hit_num) = gt_str_array_new();
        MATRIXSTRUCT(hit_dna) = gt_str_array_new();
        MATRIXSTRUCT(hit_gi_def) = gt_str_array_new();
        MATRIXSTRUCT(hit_acc) = gt_str_array_new();
        MATRIXSTRUCT(fasta_row) = gt_str_array_new();
        MATRIXSTRUCT(hit_from) = gt_str_array_new();
        MATRIXSTRUCT(hit_to) = gt_str_array_new();
        MATRIXSTRUCT(hsp_qseq) = gt_str_array_new();
        MATRIXSTRUCT(hsp_hseq) = gt_str_array_new();
        PARSESTRUCT(query_frame_tmp) = gt_str_array_new();
        PARSESTRUCT(hit_frame_tmp) = gt_str_array_new();

        /* Fuer den Fall eines Parse-Fehlers - Das Flag zeigt an, dass die
           Strings angelegt und nicht geloescht wurden */
        PARSESTRUCT(xml_tag_flag) = MGTH_SET;

        /* Auslesen der Eintrags-Nr aus der Hashtabelle */
        if ((query_nr_p = gt_cstr_nofree_ulp_gt_hashmap_get(
               PARSESTRUCT(queryhash), gt_str_get(PARSESTRUCT(buf_ptr)))))
        {
          char *seq;
          query_nr = **query_nr_p;
          seq = gt_bioseq_get_sequence(PARSESTRUCT(queryseq), query_nr);
          /* Abspeichern der zur Query-Def passenden Query-DNA */
          gt_str_append_cstr_nt(MATRIXSTRUCT(query_dna), seq,
                                gt_bioseq_get_sequence_length(PARSESTRUCT
                                                              (queryseq),
                                                              query_nr));
          gt_free(seq);
          mg_outputwriter(parsestruct_ptr, NULL, NULL, NULL, 'q', err);
        }
        else
        {
          gt_error_set(err,
                    "query-dna entry in xml-file does not exist in\
                     query-hash. wrong query-dna file?");
          PARSESTRUCT(had_err) = -1;
        }
      }
      /* Hit_id liefert die GI-Nummer, die fuer die efetch-NCBI-Abfrage
         benoetigt wird */
      else
        if ((strcmp(name, gt_str_array_get(PARSESTRUCT(hit_array), 0)) == 0))
      {
        /* loeschen des alten hit_gi_nr Eintrages */
        gt_str_reset(PARSESTRUCT(hit_gi_nr_tmp));
        gt_str_reset(PARSESTRUCT(fasta_row));

        gt_str_set(PARSESTRUCT(fasta_row), gt_str_get(PARSESTRUCT(buf_ptr)));

        /* gi_ptr zeigt auf die erste Zahl der Hit-GI-Nummer */
        gi_ptr = strchr(gt_str_get(PARSESTRUCT(buf_ptr)), '|');

        gt_assert(gi_ptr != NULL);
        gi_ptr = gi_ptr + 1;

        if (!isalpha(*gi_ptr))
        {
          /* bestimmen der GI-Nummer Laenge; strspn liest solange Zeichen,
             bis ein Zeichen auftritt, dass nicht im angegebenen String
             auftaucht - hier: '|' */
          gi_len = strspn(gi_ptr + 1, "0123456789") + 1;

          /* Die GI-Nummer der Laenge gi_len wird abgespeichert */
          gt_str_append_cstr_nt(PARSESTRUCT(hit_gi_nr_tmp), gi_ptr, gi_len);
        }
        else
        {
          gt_error_set(err,
                    "incorrect gi-hit-number in xmlfile - required format\
                     is gi|[0-9]");
          PARSESTRUCT(had_err) = -1;
        }
      }
      /* einlesen der hit_gi_def */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_array), 1)) == 0
                                         && PARSESTRUCT(giexp_flag))
      {
        gt_str_set(PARSESTRUCT(gi_def_tmp), gt_str_get(PARSESTRUCT(buf_ptr)));
      }
      /* einlesen der hit_acc_nr */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_array), 2)) == 0
                                         && PARSESTRUCT(giexp_flag))
      {
        gt_str_set(PARSESTRUCT(gi_acc_tmp), gt_str_get(PARSESTRUCT(buf_ptr)));
      }
      /* abspeichern der Hit-GI-Def, Hit-ACC-Nr. und der Hsp-num -
         Kombination ist Bestandteil der eindeutigen Hit-FASTA-Zeile im
         Hit-File und dient der Identifizierung der Hit-DNA-Sequenzen im
         XML-File */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 0)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* Die GI-Nummer wird entsprechend ihrer Laenge gi_len
           abgespeichert */
        gt_str_array_add_cstr(MATRIXSTRUCT(hit_gi_nr),
                          gt_str_get(PARSESTRUCT(hit_gi_nr_tmp)));
        /* einlesen der aktuellen hit_gi_def */
        gt_str_array_add_cstr(MATRIXSTRUCT(hit_gi_def),
                          gt_str_get(PARSESTRUCT(gi_def_tmp)));
        /* einlesen der aktuellen hit_accession-number */
        gt_str_array_add_cstr(MATRIXSTRUCT(hit_acc),
                          gt_str_get(PARSESTRUCT(gi_acc_tmp)));
        gt_str_array_add_cstr(MATRIXSTRUCT(fasta_row),
                          gt_str_get(PARSESTRUCT(fasta_row)));
        /* einlesen der aktuellen hit_num */
        gt_str_array_add_cstr(MATRIXSTRUCT(hit_num),
                          gt_str_get(PARSESTRUCT(buf_ptr)));
      }
      /* Der Query-Start-Wert wird gespeichert */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 1)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* Der Query-from Wert wird als Long-Wert gespeichert, dazu
           zunaechst Umwandlung des Strings mittels atol */
        ulong_numb_buf = atol(gt_str_get(PARSESTRUCT(buf_ptr)));
        /* Der Query-from Wert wird zum GtArray query_from hinzugefuegt */
        gt_array_add_elem(MATRIXSTRUCT(query_from), &ulong_numb_buf,
                       sizeof (unsigned long));
      }
      /* Query-Stop-Wert wird gespeichert/Bearbeitung siehe
         Query-Start-Wert */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 2)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        ulong_numb_buf = atol(gt_str_get(PARSESTRUCT(buf_ptr)));
        gt_array_add_elem(MATRIXSTRUCT(query_to), &ulong_numb_buf,
                       sizeof (unsigned long));
      }
      /* Hit-from XML-Tag */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 3)) ==
               0)
      {
        if (PARSESTRUCT(giexp_flag))
        {
          /* Speichern des Hit-from Wertes in der matrix_info Struktur
             innerhalb der parsestruct-Struktur */
          gt_str_array_add_cstr(MATRIXSTRUCT(hit_from),
                               gt_str_get(PARSESTRUCT(buf_ptr)));
        }
        else
        {
          gt_file_xprintf(HITFILEOUT, "%s ",
                             gt_str_get(PARSESTRUCT(hit_gi_nr_tmp)));
          gt_file_xprintf(HITFILEOUT, "%s ",
                             gt_str_get(PARSESTRUCT(buf_ptr)));
        }
      }
      /* Hit-to XML-Tag - Bearbeitung siehe Hit-from-Tag */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 4)) ==
               0)
      {
        if (PARSESTRUCT(giexp_flag))
        {
          /* Speichern des Hit-from Wertes in der matrix_info Struktur
             innerhalb der parsestruct-Struktur */
          gt_str_array_add_cstr(MATRIXSTRUCT(hit_to),
                            gt_str_get(PARSESTRUCT(buf_ptr)));
        }
        else
        {
          gt_file_xprintf(HITFILEOUT, "%s \n",
                             gt_str_get(PARSESTRUCT(buf_ptr)));
        }
      }
      /* Query-Frame XML-Tag; bei der Berechnung der Combined-Scores
         bilden die Hits Cluster entsprechend der Query-Frames;
         Speicherung der Frames erfolgt als Long-Wert */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 5)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* abspeichern des query_frames fuer die FASTA-File Zeile im
           Hit-File als String */
        gt_str_array_add_cstr(PARSESTRUCT(query_frame_tmp),
                          gt_str_get(PARSESTRUCT(buf_ptr)));

        /* abspeichern des query_frames fuer spaetere Sequenzberechnnugen
           als long Value */
        numb_buf = atol(gt_str_get(PARSESTRUCT(buf_ptr)));
        gt_array_add_elem(MATRIXSTRUCT(query_frame), &numb_buf, sizeof (long));
      }
      /* Hit-Frame XML-Tag/Bearbeitung siehe Query-Frame XML-Tag als
         String */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 6)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* abspeichern des hit_frames fuer die FASTA-File Zeile im
           Hit-File */
        gt_str_array_add_cstr(PARSESTRUCT(hit_frame_tmp),
                          gt_str_get(PARSESTRUCT(buf_ptr)));

        /* abspeichern des hit_frames fuer spaetere Sequenzberechnnugen
           als long Value */
        numb_buf = atol(gt_str_get(PARSESTRUCT(buf_ptr)));
        gt_array_add_elem(MATRIXSTRUCT(hit_frame), &numb_buf, sizeof (long));

        /* Wenn ein Hit-FASTA-File vorliegt existiert eine GtBioseq-Struktur
           und eine Hashtabelle, ueber die die Hit-Sequenz-Informationen
           eingelesen werden */

        unsigned long hit_nr = 0;

        GtStr *hit_tmp;
        GtStr *hit_dna_tmp;

        hit_tmp = gt_str_new();
        hit_dna_tmp = gt_str_new();

        /* Die Fasta-Zeile im Hit-File besteht aus der Hit-Gi-Def, der
           Accession-Nr, der Hit-Hsp-Nr, dem Hit-From und Hit-To Wert
           sowie des Query- und Hit-Frames getrennt durch ein
           Leerzeichen; Zeile muss eindeutig sein - durch diese
           Kombination ist dies gewaehrleistet */

        gt_str_set(hit_tmp,
                gt_str_array_get(MATRIXSTRUCT(hit_gi_nr),
                             XMLPARSERSTRUCT(hit_counter)));
        gt_str_append_cstr(hit_tmp, " ");
        gt_str_append_cstr(hit_tmp,
                        gt_str_array_get(MATRIXSTRUCT(hit_from),
                                     XMLPARSERSTRUCT(hit_counter)));
        gt_str_append_cstr(hit_tmp, " ");
        gt_str_append_cstr(hit_tmp,
                        gt_str_array_get(MATRIXSTRUCT(hit_to),
                                     XMLPARSERSTRUCT(hit_counter)));
        gt_str_append_cstr(hit_tmp, " ");
        gt_str_append_cstr(hit_tmp,
                        gt_str_array_get(MATRIXSTRUCT(fasta_row),
                                     XMLPARSERSTRUCT(hit_counter)));
        gt_str_append_cstr(hit_tmp, " ");
        gt_str_append_cstr(hit_tmp,
                        gt_str_array_get(MATRIXSTRUCT(hit_gi_def),
                                     XMLPARSERSTRUCT(hit_counter)));

        /* Hit-Hashtabelle enthaelt den konstruierten Eintrag */
        if ((query_nr_p = gt_cstr_nofree_ulp_gt_hashmap_get(
               PARSESTRUCT(hithash), gt_str_get(hit_tmp))))
        {
          char *seq;
          /* Positionsbestimmung des Eintrages in der GtBioseq-Struktur */
          hit_nr = **query_nr_p;
          seq = gt_bioseq_get_sequence(PARSESTRUCT(hitseq), hit_nr);
          /* auslesen der Sequenzinformation */
          gt_str_append_cstr_nt(hit_dna_tmp, seq,
                                gt_bioseq_get_sequence_length(PARSESTRUCT
                                                              (hitseq),
                                                              hit_nr));
          gt_free(seq);
          /* abspeichern der Hit-DNA in der Matrix-Info Struktur */
          gt_str_array_add_cstr(MATRIXSTRUCT(hit_dna), gt_str_get(hit_dna_tmp));
        }
        /* Falls kein Eintrag in der Hashtabelle gefunden wurde,
           Fehlercode setzen - falsche Hit-DNA-Datei als Parameter bei
           Programmaufruf angegeben ? */
        else
        {
          /* Hit-GI-Nr nicht in der Hashtabelle, also ueber cURL von
             z.B. NCBI nachladen */
#ifdef CURLDEF
          PARSESTRUCT(had_err) =
            mg_curl(parsestruct_ptr, XMLPARSERSTRUCT(hit_counter), err);
#else
          PARSESTRUCT(gi_flag) = 1;
#endif
        }
        gt_str_delete(hit_tmp);
        gt_str_delete(hit_dna_tmp);
      }
      /* Einlesen der translatierten Query-DNA-Sequenz */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 7)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        gt_str_array_add_cstr(MATRIXSTRUCT(hsp_qseq),
                          gt_str_get(PARSESTRUCT(buf_ptr)));
      }
      /* Einlesen der translatierten Hit-DNA-Sequenz */
      else if (strcmp(name, gt_str_array_get(PARSESTRUCT(hit_hsp_array), 8)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        gt_str_array_add_cstr(MATRIXSTRUCT(hsp_hseq),
                          gt_str_get(PARSESTRUCT(buf_ptr)));
        /* Zaehler fuer die Hits pro GI-Nummer wird um 1 erhoeht */
        XMLPARSERSTRUCT(hit_counter)++;

        if (PARSESTRUCT(gi_flag))
        {
           PARSESTRUCT(gi_flag) = 0;

           gt_str_array_set_size(MATRIXSTRUCT(hit_gi_nr),
                             gt_str_array_size(MATRIXSTRUCT(hit_gi_nr))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hit_gi_def),
                             gt_str_array_size(MATRIXSTRUCT(hit_gi_def))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hit_acc),
                             gt_str_array_size(MATRIXSTRUCT(hit_acc))-1);
           gt_str_array_set_size(MATRIXSTRUCT(fasta_row),
                             gt_str_array_size(MATRIXSTRUCT(fasta_row))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hit_num),
                             gt_str_array_size(MATRIXSTRUCT(hit_num))-1);
           gt_array_set_size(MATRIXSTRUCT(query_from),
                          gt_array_size(MATRIXSTRUCT(query_from))-1);
           gt_array_set_size(MATRIXSTRUCT(query_to),
                          gt_array_size(MATRIXSTRUCT(query_to))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hit_from),
                             gt_str_array_size(MATRIXSTRUCT(hit_from))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hit_to),
                             gt_str_array_size(MATRIXSTRUCT(hit_to))-1);
           gt_array_set_size(MATRIXSTRUCT(query_frame),
                          gt_array_size(MATRIXSTRUCT(query_frame))-1);
           gt_array_set_size(MATRIXSTRUCT(hit_frame),
                          gt_array_size(MATRIXSTRUCT(hit_frame))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hsp_qseq),
                             gt_str_array_size(MATRIXSTRUCT(hsp_qseq))-1);
           gt_str_array_set_size(MATRIXSTRUCT(hsp_hseq),
                             gt_str_array_size(MATRIXSTRUCT(hsp_hseq))-1);

           XMLPARSERSTRUCT(hit_counter)--;
        }
      }

      /* Flagberechnungen */
      if (strcmp
          (name,
           gt_str_array_get(PARSESTRUCT(query_array),
                        XMLPARSERSTRUCT(query_array_index_end))) == 0
          && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(query_array_index_end),
                    &PARSESTRUCT(def_flag), 'q');
      }
      else
        if (strcmp
            (name,
             gt_str_array_get(PARSESTRUCT(hit_array),
                          XMLPARSERSTRUCT(hit_array_index_end))) == 0
            && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(hit_array_index_end),
                    &PARSESTRUCT(hit_flag), 'h');
      }
      else
        if (strcmp
            (name,
             gt_str_array_get(PARSESTRUCT(hit_hsp_array),
                          XMLPARSERSTRUCT(hit_hsp_array_index_end))) == 0
            && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(hit_hsp_array_index_end),
                    &PARSESTRUCT(hit_hsp_flag), 't');
      }

      /* Der Lesepuffer fuer den Text zwischen 2 XML-Tags wird resetet */
      gt_str_reset(PARSESTRUCT(buf_ptr));
    }
  }
}

static void textElement(void *data, const XML_Char *txt_element, int len)
{
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  if (!PARSESTRUCT(had_err))
  {
    /* falls ein Flag gesetzt ist (relevanter XML-Tag Zwischenbereich),
       wird mit der Bearbeitung der Textpassage begonnen; dazu wird der
       "aktuelle Text" txt_element an den bereits eingelesenen Text angehaengt
     */
    if (PARSESTRUCT(hit_flag) == MGTH_SET || PARSESTRUCT(def_flag) == MGTH_SET
        || PARSESTRUCT(hit_hsp_flag) == MGTH_SET)
    {
      gt_str_append_cstr_nt(PARSESTRUCT(buf_ptr), txt_element, len);
    }
  }
}

static void flag_setting(unsigned short *counter_fct, unsigned short *flag,
                         char flag_sign)
{
  /* das uebergebene Flag wird auf SET gesetzt und der Counter um eine
     Stelle auf den naechsten Tag-Namen in der Struktur weitergesetzt */
  *flag = MGTH_SET;
  (*counter_fct)++;

  /* ueberpruefen, ob der Counter die Grenzen zur naechsten Kategorie
     erreicht hat */
  check_counter(counter_fct, flag_sign);
}

static void flag_delete(unsigned short *counter_fct, unsigned short *flag,
                        char flag_sign)
{
  /* das uebergebene Flag wird auf UNSET gesetzt und der Counter um eine
     Stelle auf den naechsten Tag-Namen in der Struktur weitergesetzt */
  *flag = MGTH_UNSET;
  (*counter_fct)++;

  /* ueberpruefen, ob der Counter die Grenzen zur neachsten Kategorie
     erreicht hat */
  check_counter(counter_fct, flag_sign);
}

static void check_counter(unsigned short *counter_check,
                          char flag_sign_fct)
{
  /* ist das Ende der aktuellen XML-Tag Kategorie erreicht, werden die
     Counter wieder auf den Anfang (0) zurueckgesetzt; das flag_sign zeigt
     an,um welchen counter es sich aktuell handelt */
  if (*counter_check == MGTH_QUERY_SIZE && flag_sign_fct == 'q')
  {
    *counter_check = 0;
  }
  if (*counter_check == MGTH_HIT_SIZE && flag_sign_fct == 'h')
  {
    *counter_check = 0;
  }
  if (*counter_check == MGTH_HIT_TO_QUERY_SIZE && flag_sign_fct == 't')
  {
    *counter_check = 0;
  }
}
