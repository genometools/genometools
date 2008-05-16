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
#include "libgtcore/unused.h"
#include "metagenomethreader.h"

#ifdef CURLDEF
#include <curl.h>
#endif

/* Expat-Funktion zur Verarbeitung des Textes zwischen XML-Tags
   Parameter: void-Zeiger auf die (Nutzer)Daten - hier: parsestruct-Struktur;
              Zeiger auf das erste Zeichen nach einem XML-Tag;
              int-Wert der Laenge des Textes zwischen den XML-Tags
   Returnwert: void */
void textElement(void *, const XML_Char *, int);

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

int mg_xmlparser(ParseStruct *parsestruct_ptr, GenFile * fp_xmlfile,
                 Error * err)
{
  int had_err = 0;

  /* Expat XML-Error setzen */
  enum XML_Error error;

  /* Puffer zum zeilenweisen Einlesen des XML-Files, XML_Parser
     deklarieren */
  Str *buf;
  XML_Parser parser;

  /* Check Umgebungsvariablen */
  error_check(err);

  /* Puffer-String anlegen */
  buf = str_new();

  /* XML-Parser wird initialisiert */
  parser = XML_ParserCreate(NULL);
  /* Die Struktur parsestruct wird als UserData gesetzt */
  XML_SetUserData(parser, parsestruct_ptr);
  /* Start- und Endelement Handler werden gesetzt */
  XML_SetElementHandler(parser, startElement, endElement);
  /* Text-Handler setzen */
  XML_SetCharacterDataHandler(parser, textElement);

  while ((str_read_next_line_generic(buf, fp_xmlfile) != EOF) && !had_err)
  {
    PARSESTRUCT(xml_linenumber)++;
    if (PARSESTRUCT(had_err))
    {
      had_err = -1;
    }

    if ((XML_Parse(parser, str_get(buf), str_length(buf), false) ==
         XML_STATUS_ERROR) && !had_err)
    {
      error = XML_GetErrorCode(parser);
      error_set(err,
                "an error occurred parsing line %lu of file \"%s\": %s",
                PARSESTRUCT(xml_linenumber), str_get(PARSESTRUCT(xmlfile)),
                XML_ErrorString(error));

      had_err = -1;
    }
    /* reset line buffer */
    str_reset(buf);
  }

  /* finish parsing */
  if ((XML_Parse(parser, NULL, 0, true) == XML_STATUS_ERROR) && !had_err)
  {
    error = XML_GetErrorCode(parser);
    error_set(err,
              "an error occurred while finishing the parsing of file\
               \"%s\": %s",
              str_get(PARSESTRUCT(xmlfile)), XML_ErrorString(error));

    had_err = -1;
  }

  if (PARSESTRUCT(xml_tag_flag) && !(!had_err) && PARSESTRUCT(giexp_flag))
  {
    strarray_delete(MATRIXSTRUCT(hit_gi_nr));
    strarray_delete(MATRIXSTRUCT(hit_num));
    strarray_delete(MATRIXSTRUCT(hit_dna));
    strarray_delete(MATRIXSTRUCT(hit_gi_def));
    strarray_delete(MATRIXSTRUCT(hit_acc));
    strarray_delete(MATRIXSTRUCT(fasta_row));
    strarray_delete(MATRIXSTRUCT(hit_from));
    strarray_delete(MATRIXSTRUCT(hit_to));
    strarray_delete(MATRIXSTRUCT(hsp_qseq));
    strarray_delete(MATRIXSTRUCT(hsp_hseq));
    strarray_delete(PARSESTRUCT(query_frame_tmp));
    strarray_delete(PARSESTRUCT(hit_frame_tmp));
  }

  /* Freigeben des XML-Parser und loeschen des Puffer-Strings */
  XML_ParserFree(parser);
  str_delete(buf);

  return had_err;
}

static void XMLCALL startElement(void *data, const char *name,
                                 UNUSED const char **atts)
{
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  if (!PARSESTRUCT(had_err))
  {
    /* Vergleich des aktuellen XML-Tags mit dem aktuellen Tag-Namen im
       jew. Array bei Uebereinstimmung wird das jeweilige Flag gesetzt */
    if (strcmp
        (name,
         strarray_get(PARSESTRUCT(query_array),
                      XMLPARSERSTRUCT(query_array_index_start))) == 0)
    {
      flag_setting(&XMLPARSERSTRUCT(query_array_index_start),
                   &PARSESTRUCT(def_flag), 'q');
    }
    else
      if (strcmp
          (name,
           strarray_get(PARSESTRUCT(hit_array),
                        XMLPARSERSTRUCT(hit_array_index_start))) == 0)
    {
      flag_setting(&XMLPARSERSTRUCT(hit_array_index_start),
                   &PARSESTRUCT(hit_flag), 'h');
    }
    else
      if (strcmp
          (name,
           strarray_get(PARSESTRUCT(hit_hsp_array),
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
    Error *err = PARSESTRUCT(err);

    /* Laenge der GI-Nr in einem GI-Def-XML-Eintrag */
    unsigned short gi_len = 0;

    /* Temp-Variablen zum Zwischenspeichern der Query-Start- bzw.
       -End-Werte sowie der Frame-Informationen */
    unsigned long ulong_numb_buf = 0,
      query_nr = 0;
    long numb_buf = 0;

    /* Zeiger auf die erste Zahl der GI-Nr in einem Hit-ID-XML-Eintrag */
    const char *gi_ptr = NULL;

    /* Check Umgebungsvariablen */
    error_check(err);

    /* wurde das Iteration_hits XML-Tag erreicht wird die mg_combinedscore
       Methode aufgerufen und diverse Variablen fuer einen neuen Eintrag
       resetet bzw. geloescht */
    if (strcmp(name, str_get(PARSESTRUCT(xml_tag))) == 0
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

        array_reset(MATRIXSTRUCT(query_from));
        array_reset(MATRIXSTRUCT(query_to));
        array_reset(MATRIXSTRUCT(hit_frame));
        array_reset(MATRIXSTRUCT(query_frame));
      }

      /* Schreiben der schliessenden XML-Tags nach der Hit-Bearbeitung */
      if (ARGUMENTSSTRUCT(outputfile_format) == 3)
      {
        /* Abschluss des Iteration-Bereichs im XML-File */
        mg_outputwriter(parsestruct_ptr, NULL, NULL, NULL, 'x', err);
      }

      strarray_delete(MATRIXSTRUCT(hit_gi_nr));
      strarray_delete(MATRIXSTRUCT(hit_gi_def));
      strarray_delete(MATRIXSTRUCT(hit_acc));
      strarray_delete(MATRIXSTRUCT(fasta_row));
      strarray_delete(MATRIXSTRUCT(hit_num));

      strarray_delete(MATRIXSTRUCT(hit_dna));
      strarray_delete(MATRIXSTRUCT(hit_from));
      strarray_delete(MATRIXSTRUCT(hit_to));
      strarray_delete(MATRIXSTRUCT(hsp_qseq));
      strarray_delete(MATRIXSTRUCT(hsp_hseq));
      strarray_delete(PARSESTRUCT(query_frame_tmp));
      strarray_delete(PARSESTRUCT(hit_frame_tmp));

      PARSESTRUCT(xml_tag_flag) = UNSET;
    }
    /* nur wenn ein Flag gesetzt ist und Text eingelesen wurde, ist eine
       weitere Bearbeitung notwendig */
    if ((PARSESTRUCT(def_flag) == 1 || PARSESTRUCT(hit_flag) == 1
         || PARSESTRUCT(hit_hsp_flag) == 1))
    {
      /* Anhand der Query-GI-Def kann mittels Hash-Table und
         Bioseq-Struktur die Query-DNA ausgelesen werden */
      if (strcmp(name, strarray_get(PARSESTRUCT(query_array), 0)) == 0
                                    && PARSESTRUCT(giexp_flag))
      {
        /* Query-DNA-Strings fuer die Sequenz und die Definition werden
           geloescht, bevor der neue Eintrag hineinkopiert wird */
        str_reset(MATRIXSTRUCT(query_dna));
        str_reset(MATRIXSTRUCT(query_def));

        /* Abspeichern der Query-Def */
        str_set(MATRIXSTRUCT(query_def), str_get(PARSESTRUCT(buf_ptr)));

        /* Erzeugen von StringArrays zur Aufnahme der
           Hit-DNA-Sequenz-Informationen */
        MATRIXSTRUCT(hit_gi_nr) = strarray_new();
        MATRIXSTRUCT(hit_num) = strarray_new();
        MATRIXSTRUCT(hit_dna) = strarray_new();
        MATRIXSTRUCT(hit_gi_def) = strarray_new();
        MATRIXSTRUCT(hit_acc) = strarray_new();
        MATRIXSTRUCT(fasta_row) = strarray_new();
        MATRIXSTRUCT(hit_from) = strarray_new();
        MATRIXSTRUCT(hit_to) = strarray_new();
        MATRIXSTRUCT(hsp_qseq) = strarray_new();
        MATRIXSTRUCT(hsp_hseq) = strarray_new();
        PARSESTRUCT(query_frame_tmp) = strarray_new();
        PARSESTRUCT(hit_frame_tmp) = strarray_new();

        /* Fuer den Fall eines Parse-Fehlers - Das Flag zeigt an, dass die
           Strings angelegt und nicht geloescht wurden */
        PARSESTRUCT(xml_tag_flag) = SET;

        if (hashtable_get
            (PARSESTRUCT(queryhash), str_get(PARSESTRUCT(buf_ptr))))
        {
          /* Auslesen der Eintrags-Nr aus der Hashtabelle */
          query_nr =
            *(unsigned long *) hashtable_get(PARSESTRUCT(queryhash),
                                             str_get(PARSESTRUCT
                                                     (buf_ptr)));

          /* Abspeichern der zur Query-Def passenden Query-DNA */
          str_append_cstr_nt(MATRIXSTRUCT(query_dna),
                             bioseq_get_sequence(PARSESTRUCT(queryseq),
                                                 query_nr),
                             bioseq_get_sequence_length(PARSESTRUCT
                                                        (queryseq),
                                                        query_nr));

          mg_outputwriter(parsestruct_ptr, NULL, NULL, NULL, 'q', err);
        }
        else
        {
          error_set(err,
                    "query-dna entry in xml-file does not exist in\
                     query-hash. wrong query-dna file?");
          PARSESTRUCT(had_err) = -1;
        }
      }
      /* Hit_id liefert die GI-Nummer, die fuer die efetch-NCBI-Abfrage
         benoetigt wird */
      else
        if ((strcmp(name, strarray_get(PARSESTRUCT(hit_array), 0)) == 0))
      {
        /* loeschen des alten hit_gi_nr Eintrages */
        str_reset(PARSESTRUCT(hit_gi_nr_tmp));
        str_reset(PARSESTRUCT(fasta_row));

        str_set(PARSESTRUCT(fasta_row), str_get(PARSESTRUCT(buf_ptr)));

        /* gi_ptr zeigt auf die erste Zahl der Hit-GI-Nummer */
        gi_ptr = strchr(str_get(PARSESTRUCT(buf_ptr)), '|');

        assert(gi_ptr != NULL);
        gi_ptr = gi_ptr + 1;

        if (!isalpha(*gi_ptr))
        {
          /* bestimmen der GI-Nummer Laenge; strspn liest solange Zeichen,
             bis ein Zeichen auftritt, dass nicht im angegebenen String
             auftaucht - hier: '|' */
          gi_len = strspn(gi_ptr + 1, "0123456789") + 1;

          /* Die GI-Nummer der Laenge gi_len wird abgespeichert */
          str_append_cstr_nt(PARSESTRUCT(hit_gi_nr_tmp), gi_ptr, gi_len);
        }
        else
        {
          error_set(err,
                    "incorrect gi-hit-number in xmlfile - required format\
                     is gi|[0-9]");
          PARSESTRUCT(had_err) = -1;
        }
      }
      /* einlesen der hit_gi_def */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_array), 1)) == 0
                                         && PARSESTRUCT(giexp_flag))
      {
        str_set(PARSESTRUCT(gi_def_tmp), str_get(PARSESTRUCT(buf_ptr)));
      }
      /* einlesen der hit_acc_nr */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_array), 2)) == 0
                                         && PARSESTRUCT(giexp_flag))
      {
        str_set(PARSESTRUCT(gi_acc_tmp), str_get(PARSESTRUCT(buf_ptr)));
      }
      /* abspeichern der Hit-GI-Def, Hit-ACC-Nr. und der Hsp-num -
         Kombination ist Bestandteil der eindeutigen Hit-FASTA-Zeile im
         Hit-File und dient der Identifizierung der Hit-DNA-Sequenzen im
         XML-File */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 0)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* Die GI-Nummer wird entsprechend ihrer Laenge gi_len
           abgespeichert */
        strarray_add_cstr(MATRIXSTRUCT(hit_gi_nr),
                          str_get(PARSESTRUCT(hit_gi_nr_tmp)));
        /* einlesen der aktuellen hit_gi_def */
        strarray_add_cstr(MATRIXSTRUCT(hit_gi_def),
                          str_get(PARSESTRUCT(gi_def_tmp)));
        /* einlesen der aktuellen hit_accession-number */
        strarray_add_cstr(MATRIXSTRUCT(hit_acc),
                          str_get(PARSESTRUCT(gi_acc_tmp)));
        strarray_add_cstr(MATRIXSTRUCT(fasta_row),
                          str_get(PARSESTRUCT(fasta_row)));
        /* einlesen der aktuellen hit_num */
        strarray_add_cstr(MATRIXSTRUCT(hit_num),
                          str_get(PARSESTRUCT(buf_ptr)));
      }
      /* Der Query-Start-Wert wird gespeichert */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 1)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* Der Query-from Wert wird als Long-Wert gespeichert, dazu
           zunaechst Umwandlung des Strings mittels atol */
        ulong_numb_buf = atol(str_get(PARSESTRUCT(buf_ptr)));
        /* Der Query-from Wert wird zum Array query_from hinzugefuegt */
        array_add_elem(MATRIXSTRUCT(query_from), &ulong_numb_buf,
                       sizeof (unsigned long));
      }
      /* Query-Stop-Wert wird gespeichert/Bearbeitung siehe
         Query-Start-Wert */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 2)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        ulong_numb_buf = atol(str_get(PARSESTRUCT(buf_ptr)));
        array_add_elem(MATRIXSTRUCT(query_to), &ulong_numb_buf,
                       sizeof (unsigned long));
      }
      /* Hit-from XML-Tag */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 3)) ==
               0)
      {
        if (PARSESTRUCT(giexp_flag))
        {
          /* Speichern des Hit-from Wertes in der matrix_info Struktur
             innerhalb der parsestruct-Struktur */
          strarray_add_cstr(MATRIXSTRUCT(hit_from),
                            str_get(PARSESTRUCT(buf_ptr)));
        }
        else
        {
          genfile_xprintf(HITFILEOUT, "%s ",
                          str_get(PARSESTRUCT(hit_gi_nr_tmp)));
          genfile_xprintf(HITFILEOUT, "%s ", str_get(PARSESTRUCT(buf_ptr)));
        }
      }
      /* Hit-to XML-Tag - Bearbeitung siehe Hit-from-Tag */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 4)) ==
               0)
      {
        if (PARSESTRUCT(giexp_flag))
        {
          /* Speichern des Hit-from Wertes in der matrix_info Struktur
             innerhalb der parsestruct-Struktur */
          strarray_add_cstr(MATRIXSTRUCT(hit_to),
                            str_get(PARSESTRUCT(buf_ptr)));
        }
        else
        {
          genfile_xprintf(HITFILEOUT, "%s \n", str_get(PARSESTRUCT(buf_ptr)));
        }
      }
      /* Query-Frame XML-Tag; bei der Berechnung der Combined-Scores
         bilden die Hits Cluster entsprechend der Query-Frames;
         Speicherung der Frames erfolgt als Long-Wert */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 5)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* abspeichern des query_frames fuer die FASTA-File Zeile im
           Hit-File als String */
        strarray_add_cstr(PARSESTRUCT(query_frame_tmp),
                          str_get(PARSESTRUCT(buf_ptr)));

        /* abspeichern des query_frames fuer spaetere Sequenzberechnnugen
           als long Value */
        numb_buf = atol(str_get(PARSESTRUCT(buf_ptr)));
        array_add_elem(MATRIXSTRUCT(query_frame), &numb_buf, sizeof (long));
      }
      /* Hit-Frame XML-Tag/Bearbeitung siehe Query-Frame XML-Tag als
         String */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 6)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        /* abspeichern des hit_frames fuer die FASTA-File Zeile im
           Hit-File */
        strarray_add_cstr(PARSESTRUCT(hit_frame_tmp),
                          str_get(PARSESTRUCT(buf_ptr)));

        /* abspeichern des hit_frames fuer spaetere Sequenzberechnnugen
           als long Value */
        numb_buf = atol(str_get(PARSESTRUCT(buf_ptr)));
        array_add_elem(MATRIXSTRUCT(hit_frame), &numb_buf, sizeof (long));

        /* Wenn ein Hit-FASTA-File vorliegt existiert eine Bioseq-Struktur
           und eine Hashtabelle, ueber die die Hit-Sequenz-Informationen
           eingelesen werden */

        unsigned long hit_nr = 0;

        Str *hit_tmp;
        Str *hit_dna_tmp;

        hit_tmp = str_new();
        hit_dna_tmp = str_new();

        /* Die Fasta-Zeile im Hit-File besteht aus der Hit-Gi-Def, der
           Accession-Nr, der Hit-Hsp-Nr, dem Hit-From und Hit-To Wert
           sowie des Query- und Hit-Frames getrennt durch ein
           Leerzeichen; Zeile muss eindeutig sein - durch diese
           Kombination ist dies gewaehrleistet */

        str_set(hit_tmp,
                strarray_get(MATRIXSTRUCT(hit_gi_nr),
                             XMLPARSERSTRUCT(hit_counter)));
        str_append_cstr(hit_tmp, " ");
        str_append_cstr(hit_tmp,
                        strarray_get(MATRIXSTRUCT(hit_from),
                                     XMLPARSERSTRUCT(hit_counter)));
        str_append_cstr(hit_tmp, " ");
        str_append_cstr(hit_tmp,
                        strarray_get(MATRIXSTRUCT(hit_to),
                                     XMLPARSERSTRUCT(hit_counter)));
        str_append_cstr(hit_tmp, " ");
        str_append_cstr(hit_tmp,
                        strarray_get(MATRIXSTRUCT(fasta_row),
                                     XMLPARSERSTRUCT(hit_counter)));
        str_append_cstr(hit_tmp, " ");
        str_append_cstr(hit_tmp,
                        strarray_get(MATRIXSTRUCT(hit_gi_def),
                                     XMLPARSERSTRUCT(hit_counter)));

        /* Hit-Hashtabelle enthaelt den konstruierten Eintrag */
        if (hashtable_get(PARSESTRUCT(hithash), str_get(hit_tmp)))
        {
          /* Positionsbestimmung des Eintrages in der Bioseq-Struktur */
          hit_nr =
            *(unsigned long *) hashtable_get(PARSESTRUCT(hithash),
                                             str_get(hit_tmp));

          /* auslesen der Sequenzinformation */
          str_append_cstr_nt(hit_dna_tmp,
                             bioseq_get_sequence(PARSESTRUCT(hitseq),
                                                 hit_nr),
                             bioseq_get_sequence_length(PARSESTRUCT
                                                        (hitseq),
                                                        hit_nr));
          /* abspeichern der Hit-DNA in der Matrix-Info Struktur */
          strarray_add_cstr(MATRIXSTRUCT(hit_dna), str_get(hit_dna_tmp));
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
#endif
#ifndef CURLDEF
          PARSESTRUCT(gi_flag) = 1;
#endif
        }
        str_delete(hit_tmp);
        str_delete(hit_dna_tmp);
      }
      /* Einlesen der translatierten Query-DNA-Sequenz */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 7)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        strarray_add_cstr(MATRIXSTRUCT(hsp_qseq),
                          str_get(PARSESTRUCT(buf_ptr)));
      }
      /* Einlesen der translatierten Hit-DNA-Sequenz */
      else if (strcmp(name, strarray_get(PARSESTRUCT(hit_hsp_array), 8)) ==
               0 && PARSESTRUCT(giexp_flag))
      {
        strarray_add_cstr(MATRIXSTRUCT(hsp_hseq),
                          str_get(PARSESTRUCT(buf_ptr)));
        /* Zaehler fuer die Hits pro GI-Nummer wird um 1 erhoeht */
        XMLPARSERSTRUCT(hit_counter)++;

        if (PARSESTRUCT(gi_flag))
        {
           PARSESTRUCT(gi_flag) = 0;

           strarray_set_size(MATRIXSTRUCT(hit_gi_nr),
                             strarray_size(MATRIXSTRUCT(hit_gi_nr))-1);
           strarray_set_size(MATRIXSTRUCT(hit_gi_def),
                             strarray_size(MATRIXSTRUCT(hit_gi_def))-1);
           strarray_set_size(MATRIXSTRUCT(hit_acc),
                             strarray_size(MATRIXSTRUCT(hit_acc))-1);
           strarray_set_size(MATRIXSTRUCT(fasta_row),
                             strarray_size(MATRIXSTRUCT(fasta_row))-1);
           strarray_set_size(MATRIXSTRUCT(hit_num),
                             strarray_size(MATRIXSTRUCT(hit_num))-1);
           array_set_size(MATRIXSTRUCT(query_from),
                          array_size(MATRIXSTRUCT(query_from))-1);
           array_set_size(MATRIXSTRUCT(query_to),
                          array_size(MATRIXSTRUCT(query_to))-1);
           strarray_set_size(MATRIXSTRUCT(hit_from),
                             strarray_size(MATRIXSTRUCT(hit_from))-1);
           strarray_set_size(MATRIXSTRUCT(hit_to),
                             strarray_size(MATRIXSTRUCT(hit_to))-1);
           array_set_size(MATRIXSTRUCT(query_frame),
                          array_size(MATRIXSTRUCT(query_frame))-1);
           array_set_size(MATRIXSTRUCT(hit_frame),
                          array_size(MATRIXSTRUCT(hit_frame))-1);
           strarray_set_size(MATRIXSTRUCT(hsp_qseq),
                             strarray_size(MATRIXSTRUCT(hsp_qseq))-1);
           strarray_set_size(MATRIXSTRUCT(hsp_hseq),
                             strarray_size(MATRIXSTRUCT(hsp_hseq))-1);

           XMLPARSERSTRUCT(hit_counter)--;
        }
      }

      /* Flagberechnungen */
      if (strcmp
          (name,
           strarray_get(PARSESTRUCT(query_array),
                        XMLPARSERSTRUCT(query_array_index_end))) == 0
          && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(query_array_index_end),
                    &PARSESTRUCT(def_flag), 'q');
      }
      else
        if (strcmp
            (name,
             strarray_get(PARSESTRUCT(hit_array),
                          XMLPARSERSTRUCT(hit_array_index_end))) == 0
            && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(hit_array_index_end),
                    &PARSESTRUCT(hit_flag), 'h');
      }
      else
        if (strcmp
            (name,
             strarray_get(PARSESTRUCT(hit_hsp_array),
                          XMLPARSERSTRUCT(hit_hsp_array_index_end))) == 0
            && !PARSESTRUCT(had_err))
      {
        flag_delete(&XMLPARSERSTRUCT(hit_hsp_array_index_end),
                    &PARSESTRUCT(hit_hsp_flag), 't');
      }

      /* Der Lesepuffer fuer den Text zwischen 2 XML-Tags wird resetet */
      str_reset(PARSESTRUCT(buf_ptr));
    }
  }
}

void textElement(void *data, const XML_Char *txt_element, int len)
{
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  if (!PARSESTRUCT(had_err))
  {
    /* falls ein Flag gesetzt ist (relevanter XML-Tag Zwischenbereich),
       wird mit der Bearbeitung der Textpassage begonnen; dazu wird der
       "aktuelle Text" txt_element an den bereits eingelesenen Text angehaengt
     */
    if (PARSESTRUCT(hit_flag) == SET || PARSESTRUCT(def_flag) == SET
        || PARSESTRUCT(hit_hsp_flag) == SET)
    {
      str_append_cstr_nt(PARSESTRUCT(buf_ptr), txt_element, len);
    }
  }
}

static void flag_setting(unsigned short *counter_fct, unsigned short *flag,
                         char flag_sign)
{
  /* das uebergebene Flag wird auf SET gesetzt und der Counter um eine
     Stelle auf den naechsten Tag-Namen in der Struktur weitergesetzt */
  *flag = SET;
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
  *flag = UNSET;
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
  if (*counter_check == QUERY_SIZE && flag_sign_fct == 'q')
  {
    *counter_check = 0;
  }
  if (*counter_check == HIT_SIZE && flag_sign_fct == 'h')
  {
    *counter_check = 0;
  }
  if (*counter_check == HIT_TO_QUERY_SIZE && flag_sign_fct == 't')
  {
    *counter_check = 0;
  }
}
