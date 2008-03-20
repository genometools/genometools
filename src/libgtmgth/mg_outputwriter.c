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

#include "libgtcore/unused.h"
#include "mg_outputwriter.h"

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des txt-Dokuments aufruft
   Parameter:  Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
               die HitInformation-Struktur, die RegionStruct-Struktur,
               das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_txt(ParseStruct *,
                             CombinedScoreMatrixEntry **,
                             HitInformation *,
                             RegionStruct **, char, Error *);

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des html-Dokuments aufruft
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_html(ParseStruct *,
                              CombinedScoreMatrixEntry **,
                              HitInformation *,
                              RegionStruct **, char, Error *);

/* Funktion, die nacheinander die erforderlichen Funktionen zur Erstellung
   des xml-Dokuments aufruft
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void outputwriter_xml(ParseStruct *,
                             CombinedScoreMatrixEntry **,
                             HitInformation *,
                             RegionStruct **, char, Error *);

/* Funktion zum Schreiben des Output-Header der txt-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_txt(const ParseStruct *);

/* Funktion zum Schreiben des Output-Header der html-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_html(const ParseStruct *);

/* Funktion zum Schreiben des Output-Header der xml-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_header_xml(const ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der txt-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_txt(const ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der html-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_html(const ParseStruct *);

/* Funktion zum Schreiben des Query-DNA Abschnittes der xml-Datei
   Parameter:  Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_querydna_xml(const ParseStruct *);

/* Funktion zur Berechnung der Ausgabemenge fuer den Hit-Information
   Abschnittes;
   Parameter: Zeiger auf ParseStruct-Struktur, CombinedScore-Matrix,
              die HitInformation-Struktur, die RegionStruct-Struktur,
              das Char-Zeichen, um welchen Bereich es sich handelt
   Returnwert: void */
static void output_hitdna(ParseStruct *,
                          CombinedScoreMatrixEntry **,
                          HitInformation *, RegionStruct **, Error *);

/* Funktion zum Schreiben des Coding-DNA-Abschnittes
   Parameter:  Zeiger auf ParseStruct, Zeiger auf die DNA-Seq,
               Str der Protein-Seq, Ausgabe-Format
   Returnwert: void */
static void print_codingheader(const ParseStruct *, const char *, Str *);

/* Funktion zum Schreiben des Hit-Information-Abschnittes
   Parameter: Zeiger auf die ParseStruct, Zeiger auf die
              HitInformation-Struktur, Zeiger auf die kodierende-Seq,
              Ausgabeformat, Sequenzposition
   Returnwert: void */
static void print_hitinformation(const ParseStruct *,
                                 const HitInformation *, unsigned long);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der Text-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_txt(const ParseStruct *);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der HTML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_html(const ParseStruct *);

/* Funktion zum Schreiben der Metagenomethreader-Statistik der XML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_xml(const ParseStruct *);

/* Funktion zum Schreiben des Output-Footer der HTML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_footer_html(const ParseStruct *);

/* Funktion zum Schreiben des Output-Footer der XML-Datei
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_footer_xml(const ParseStruct *);

/* Funktion zur Berechnung der AS-Sequenz
   Parameter: Zeiger auf die DNA-Sequnez, Str fuer die Proteinsequenz,
              from-Wert, to-Wert, aktueller Frame
   Returnwert: void */
static int as_coding(const ParseStruct *,
                     char *,
                     Str *,
                     unsigned long,
                     unsigned long, unsigned short, Error *);

static int newmemory_hash(void *, void *, void *, Error *);

/* Funktion zum Schreiben des Statistic-Headers
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_statistics_header(const ParseStruct *);

static short check_startcodon(const ParseStruct *, const char *);

/* Funktion zum Abschluss eines Iterations-Bereich im XML-File
   Parameter: Zeiger auf die ParseStruct-Struktur
   Returnwert: void */
static void output_close_iteration_xml(const ParseStruct *);

/* Zeitstruktur liefert die Struktur fuer das aktuelle Datum in der Form
   dd.mm.yyyy */
struct tm *today(void)
{
  static struct tm *tmnow;

  time_t tnow;

  (void) time(&tnow);
  tmnow = localtime(&tnow);

  return tmnow;
}

void mg_outputwriter(ParseStruct *parsestruct_ptr,
                     CombinedScoreMatrixEntry **combinedscore_matrix,
                     HitInformation *hit_information,
                     RegionStruct **regionmatrix, char type, Error * err)
{
  error_check(err);

  /* je nach dem per Option angegebenen Output-Fileformat wird die
     entsprechende Ausgabefunktion aufgerufen */
  switch (ARGUMENTSSTRUCT(outputfile_format))
  {
      /* txt */
    case 1:
      outputwriter_txt(parsestruct_ptr, combinedscore_matrix,
                       hit_information, regionmatrix, type, err);
      break;
      /* html */
    case 2:
      outputwriter_html(parsestruct_ptr, combinedscore_matrix,
                        hit_information, regionmatrix, type, err);
      break;
      /* xml */
    case 3:
      outputwriter_xml(parsestruct_ptr, combinedscore_matrix,
                       hit_information, regionmatrix, type, err);
      break;
  }
}

static void outputwriter_txt(ParseStruct *parsestruct_ptr,
                      CombinedScoreMatrixEntry **combinedscore_matrix,
                      HitInformation *hit_information,
                      RegionStruct **regionmatrix, char type, Error * err)
{
  error_check(err);

  /* Aufruf der jeweiligen Funtion zum erstellen der Bereiche Header(t),
     Query-DNA(q) und Coding-DNA(h) */
  if (type == 't')
    output_header_txt(parsestruct_ptr);
  else if (type == 'q')
    output_querydna_txt(parsestruct_ptr);
  else if (type == 'h')
    output_hitdna(parsestruct_ptr, combinedscore_matrix, hit_information,
                  regionmatrix, err);
  else if (type == 'v')
    output_statistics_header(parsestruct_ptr);
  else if (type == 's')
    output_statistics_txt(parsestruct_ptr);
}

static void outputwriter_html(ParseStruct *parsestruct_ptr,
                       CombinedScoreMatrixEntry **combinedscore_matrix,
                       HitInformation *hit_information,
                       RegionStruct **regionmatrix, char type, Error * err)
{
  error_check(err);

  /* Aufruf der jeweiligen Funtion zum erstellen der Bereiche Header(t),
     Query-DNA(q) und Coding-DNA(h) und des Footer(f) */
  if (type == 't')
    output_header_html(parsestruct_ptr);
  else if (type == 'q')
    output_querydna_html(parsestruct_ptr);
  else if (type == 'h')
    output_hitdna(parsestruct_ptr, combinedscore_matrix, hit_information,
                  regionmatrix, err);
  else if (type == 's')
    output_statistics_html(parsestruct_ptr);
  else if (type == 'v')
    output_statistics_header(parsestruct_ptr);
  else if (type == 'f')
    output_footer_html(parsestruct_ptr);
}

static void outputwriter_xml(ParseStruct *parsestruct_ptr,
                      CombinedScoreMatrixEntry **combinedscore_matrix,
                      HitInformation *hit_information,
                      RegionStruct **regionmatrix, char type, Error * err)
{
  error_check(err);

  /* Aufruf der jeweiligen Funtion zum erstellen der Bereiche Header(t),
     Query-DNA(q) und Coding-DNA(h) und des Footer(f) */
  if (type == 't')
    output_header_xml(parsestruct_ptr);
  else if (type == 'q')
    output_querydna_xml(parsestruct_ptr);
  else if (type == 'h')
    output_hitdna(parsestruct_ptr, combinedscore_matrix, hit_information,
                  regionmatrix, err);
  else if (type == 's')
    output_statistics_xml(parsestruct_ptr);
  else if (type == 'v')
    output_statistics_header(parsestruct_ptr);
  else if (type == 'f')
    output_footer_xml(parsestruct_ptr);
  else if (type == 'x')
    output_close_iteration_xml(parsestruct_ptr);
}

static void output_header_txt(const ParseStruct *parsestruct_ptr)
{
  /* Aktuelle Datum in Variable speichern */
  struct tm *tmstamp;

  tmstamp = today();

  if (!ARGUMENTSSTRUCT(testmodus_mode))
  {
    /* Headerbereich schreiben inkl. Auflistung der Parametereinstellugen */
    genfile_xprintf(FILEPOINTEROUT,
                    "\nMetagenomethreader Result %d.%d.%d\n\n",
                    tmstamp->tm_mday, tmstamp->tm_mon + 1,
                    tmstamp->tm_year + 1900);
  }
  genfile_xprintf(FILEPOINTEROUT,
                  "\nParametereinstellungen\n Synonymic Value: %.4f\n ",
                  ARGUMENTSSTRUCT(synonomic_value));
  genfile_xprintf(FILEPOINTEROUT, "Nonsynonymic Value: %.4f\n ",
                  ARGUMENTSSTRUCT(nonsynonomic_value));
  genfile_xprintf(FILEPOINTEROUT, "Blasthit-End Value: %.4f\n ",
                  ARGUMENTSSTRUCT(blasthit_end_value));
  genfile_xprintf(FILEPOINTEROUT, "Query-Stopcodon-Value: %.4f\n ",
                  ARGUMENTSSTRUCT(stopcodon_queryseq));
  genfile_xprintf(FILEPOINTEROUT, "Hit-Stopcodon-Value: %.4f\n ",
                  ARGUMENTSSTRUCT(stopcodon_hitseq));
  genfile_xprintf(FILEPOINTEROUT, "Frameshift-Span: %.4f\n ",
                  ARGUMENTSSTRUCT(frameshift_span));
  genfile_xprintf(FILEPOINTEROUT, "Prediction-Span: %.4f\n ",
                  ARGUMENTSSTRUCT(prediction_span));
  genfile_xprintf(FILEPOINTEROUT, "Leavegene-Value: %.4f\n ",
                  ARGUMENTSSTRUCT(leavegene_value));
  genfile_xprintf(FILEPOINTEROUT, "Curl-DB: %s\n ",
                  str_get(ARGUMENTSSTRUCT(curl_fcgi_db)));
  genfile_xprintf(FILEPOINTEROUT, "Output-Filename: %s\n ",
                  str_get(ARGUMENTSSTRUCT(outputtextfile_name)));
  genfile_xprintf(FILEPOINTEROUT, "Output-Fileformat: %d\n ",
                  ARGUMENTSSTRUCT(outputfile_format));
  genfile_xprintf(FILEPOINTEROUT, "Hitfile (yes=1/no=0): %d\n ",
                  ARGUMENTSSTRUCT(hitfile_bool));
  genfile_xprintf(FILEPOINTEROUT, "Min Protein-Length (>=15): %lu\n ",
                  ARGUMENTSSTRUCT(min_as));
  genfile_xprintf(FILEPOINTEROUT, "Min Result-Percentage: %.4f\n ",
                  ARGUMENTSSTRUCT(percent_value));
  genfile_xprintf(FILEPOINTEROUT, "Extended-Modus (yes=1/no=0): %d\n ",
                  ARGUMENTSSTRUCT(extended_mode));
  genfile_xprintf(FILEPOINTEROUT, "Homology-Modus (yes=1/no=0): %d\n ",
                  ARGUMENTSSTRUCT(homology_mode));
  genfile_xprintf(FILEPOINTEROUT, "Codon-Modus (yes=1/no=0): %d\n\n",
                  ARGUMENTSSTRUCT(codon_mode));
}

static void output_header_html(const ParseStruct *parsestruct_ptr)
{
  /* Aktuelle Datum in Variable speichern */
  struct tm *tmstamp;

  tmstamp = today();

  /* Headerbereich schreiben inkl. Auflistung der Parametereinstellugen */
  genfile_xprintf(FILEPOINTEROUT,
                  "<!DOCTYPE html PUBLIC \"-/""/W3C/""/DTD XHTML 1.0 "
                  "Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/"
                  "xhtml1-transitional.dtd\">");
  genfile_xprintf(FILEPOINTEROUT,
                  "<html xmlns=\"http://www.w3.org/1999/xhtml\" "
                  "xml:lang=\"de\" lang=\"de\">");
  genfile_xprintf(FILEPOINTEROUT, "<head>");
  if (!ARGUMENTSSTRUCT(testmodus_mode))
  {
    genfile_xprintf(FILEPOINTEROUT,
                    "<title>Metagenomethreader Result %d.%d.%d</title>",
                    tmstamp->tm_mday, tmstamp->tm_mon + 1,
                    tmstamp->tm_year + 1900);
  }
    genfile_xprintf(FILEPOINTEROUT,
                    "<meta http-equiv=\"Content-type\" content=\"text/html; "
                    "charset=iso-8859-1\"/>");
    genfile_xprintf(FILEPOINTEROUT,
                    "<link rel=\"stylesheet\" type=\"text/css\" "
                    "href=\"styles.css\" media=\"all\"/>");
    genfile_xprintf(FILEPOINTEROUT, "</head>");
    genfile_xprintf(FILEPOINTEROUT, "<body>");
    genfile_xprintf(FILEPOINTEROUT,
                    "<table border=\"0\" width=\"800\" cellspacing=\"1\" "
                    "cellpadding=\"2\">");
  if (!ARGUMENTSSTRUCT(testmodus_mode))
  {
    genfile_xprintf(FILEPOINTEROUT,
                    "<tr><td width=\"200\"><font class=\"font_header\">"
                    "Metagenomethreader Result %d.%d.%d</font><br><br></td>"
                    "<td></td></tr>",
                    tmstamp->tm_mday, tmstamp->tm_mon + 1,
                    tmstamp->tm_year + 1900);
  }
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Parametereinstellungen</font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\"></font></td></tr>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Synonymic Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(synonomic_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Non-Synonymic Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(nonsynonomic_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Blast-Hit-End Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(blasthit_end_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Query Stop-Codon Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(stopcodon_queryseq));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Hit Stop-Codon Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(stopcodon_hitseq));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Frameshift-Span: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(frameshift_span));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Prediction-Span: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(prediction_span));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Leavegene-Value: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(leavegene_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">cURL-DB: "
                  "</font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%s</font></td></tr>",
                  str_get(ARGUMENTSSTRUCT(curl_fcgi_db)));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Output-Filename: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%s</font></td></tr>",
                  str_get(ARGUMENTSSTRUCT(outputtextfile_name)));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Output-Fileformat<br>(1/2/3): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%d</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(outputfile_format));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">Hitfile<br>"
                  "(yes=1/no=0): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%d</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(hitfile_bool));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Min-Protein-Length<br>(>=15): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%lu</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(min_as));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Min-Result-Percentage: </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"class\">%.4f</font></td></tr>",
                  ARGUMENTSSTRUCT(percent_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Extended-Modus<br>(yes=1/no=0): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%d</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(extended_mode));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Homology-Modus<br>(yes=1/no=0): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%d</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(homology_mode));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td width=\"200\"><font class=\"class\">"
                  "Codon-Modus<br>(1/2/3): </font></td>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<td valign=\"top\"><font class=\"class\">%d</font>"
                  "</td></tr>",
                  ARGUMENTSSTRUCT(codon_mode));
}

static void output_header_xml(const ParseStruct *parsestruct_ptr)
{
  /* Aktuelle Datum in Variable speichern */
  struct tm *tmstamp;

  tmstamp = today();

  /* Headerbereich schreiben inkl. Auflistung der Parametereinstellugen */
  genfile_xprintf(FILEPOINTEROUT, "<?xml version=\"1.0\"?>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "<!DOCTYPE BlastOutput PUBLIC \"-/""/NCBI/""/NCBI "
                  "BlastOutput/EN\" \"NCBI_BlastOutput.dtd\">\n");
  genfile_xprintf(FILEPOINTEROUT, "<MetagenomethreaderOutput>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "  <MetagenomethreaderOutput_title>Metagenomethreader"
                  "</MetagenomethreaderOutput_title>\n");

  if (!ARGUMENTSSTRUCT(testmodus_mode))
  {
    genfile_xprintf(FILEPOINTEROUT,
                    "  <MetagenomethreaderOutput_date>Result %d.%d.%d"
                    "</MetagenomethreaderOutput_date>\n",
                    tmstamp->tm_mday, tmstamp->tm_mon + 1,
                    tmstamp->tm_year + 1900);
  }
  genfile_xprintf(FILEPOINTEROUT, "  <MetagenomethreaderOutput_param>\n");
  genfile_xprintf(FILEPOINTEROUT, "    <Parameters>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_syn>%.4f</Parameters_syn>\n",
                  ARGUMENTSSTRUCT(synonomic_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_nonsyn>%.4f</Parameters_nonsyn>\n",
                  ARGUMENTSSTRUCT(nonsynonomic_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_blastend>%.4f"
                  "</Parameters_blastend>\n",
                  ARGUMENTSSTRUCT(blasthit_end_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_stopcodon-query>%.4f"
                  "</Parameters_stopcodon-query>\n",
                  ARGUMENTSSTRUCT(stopcodon_queryseq));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_stopcodon-hit>%.4f"
                  "</Parameters_stopcodon-hit>\n",
                  ARGUMENTSSTRUCT(stopcodon_hitseq));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_frameshift>%.4f"
                  "</Parameters_frameshift>\n",
                  ARGUMENTSSTRUCT(frameshift_span));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_prediction>%.4f"
                  "</Parameters_prediction>\n",
                  ARGUMENTSSTRUCT(prediction_span));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_leavegene>%.4f"
                  "</Parameters_leavegene>\n",
                  ARGUMENTSSTRUCT(leavegene_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_curl-db>%s"
                  "</Parameters_curl-db>\n",
                  str_get(ARGUMENTSSTRUCT(curl_fcgi_db)));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_output-file>%s"
                  "</Parameters_output-file>\n",
                  str_get(ARGUMENTSSTRUCT(outputtextfile_name)));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_output-format>%d"
                  "</Parameters_output-format>\n",
                  ARGUMENTSSTRUCT(outputfile_format));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_hitfile>%d</Parameters_hitfile>\n",
                  ARGUMENTSSTRUCT(hitfile_bool));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_min-as>%lu</Parameters_min-as>\n",
                  ARGUMENTSSTRUCT(min_as));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_min_resultpercentage>%.4f"
                  "</Parameters_min_resultpercentage>\n",
                  ARGUMENTSSTRUCT(percent_value));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_extended_modus>%d"
                  "</Parameters_extended_modus>\n",
                  ARGUMENTSSTRUCT(extended_mode));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_homology_modus>%d"
                  "</Parameters_homology_modus>\n",
                  ARGUMENTSSTRUCT(homology_mode));
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Parameters_codon_modus>%d"
                  "</Parameters_codon_modus>\n",
                  ARGUMENTSSTRUCT(codon_mode));
  genfile_xprintf(FILEPOINTEROUT, "    </Parameters>\n");
  genfile_xprintf(FILEPOINTEROUT, "  </MetagenomethreaderOutput_param>\n");
}

static void output_querydna_txt(const ParseStruct *parsestruct_ptr)
{
  /* schreiben des Query-DNA Headers inkl. Query-Def. und Query-Sequenz */
  genfile_xprintf(FILEPOINTEROUT, "Query-DNA-Entry-Section\n\n");
  genfile_xprintf(FILEPOINTEROUT, "Query-DNA-Def: %s\n",
                  str_get(MATRIXSTRUCT(query_def)));
  genfile_xprintf(FILEPOINTEROUT, "Query_DNA-Sequence:\n%s\n",
                  str_get(MATRIXSTRUCT(query_dna)));
  genfile_xprintf(FILEPOINTEROUT, "\nCoding-DNA-Entry-Section\n\n");
}

static void output_querydna_html(const ParseStruct *parsestruct_ptr)
{
  /* schreiben des Query-DNA Headers inkl. Query-Def. und Query-Sequenz */
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td colspan=\"2\"><font class=\"font_header\"><br>"
                  "<br>Query-DNA-Entry-Section<br><br></font></td></tr>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td><font class=\"class\">Query-DNA-Def</font></td>"
                  "<td><font class=\"class\">%s</font></td></tr>",
                  str_get(MATRIXSTRUCT(query_def)));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td colspan=\"2\"><font class=\"class\">"
                  "Query_DNA-Sequence</font></td></tr>");
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td colspan=\"2\"><font class=\"class\">%s</font>"
                  "</td></tr>",
                  str_get(MATRIXSTRUCT(query_dna)));
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td colspan=\"2\"><br><font class=\"class\">"
                  "Coding-DNA-Entry-Section</font></td></tr>");
}

static void output_querydna_xml(const ParseStruct *parsestruct_ptr)
{
  /* schreiben des Query-DNA Headers inkl. Query-Def. und Query-Sequenz */
  genfile_xprintf(FILEPOINTEROUT,
                  "  <MetagenomethreaderOutput_iterations>\n");
  genfile_xprintf(FILEPOINTEROUT, "  <Iteration>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "    <Iteration_query-def>%s</Iteration_query-def>\n",
                  str_get(MATRIXSTRUCT(query_def)));
  genfile_xprintf(FILEPOINTEROUT,
                  "    <Iteration_query-dna>%s</Iteration_query-dna>\n",
                  str_get(MATRIXSTRUCT(query_dna)));
  genfile_xprintf(FILEPOINTEROUT, "    <Iteration_hits>\n");
}

static void output_hitdna(ParseStruct *parsestruct_ptr,
                   CombinedScoreMatrixEntry **combinedscore_matrix,
                   HitInformation *hit_information,
                   RegionStruct **regionmatrix, Error * err)
{
  int had_err = 0;

  unsigned long arraysize,
    from = 0,
    to = 0,
    contig_seq_diff,
    hit_numbers = 0,
    current,
    *hit_ptr,
    hitcounter = 0,
    array_index,
    seq_index,
    hit_index,
    string_number,
    tmp_var,
    *hitsnum_tmp,
    hit_from,
    hit_to,
    hash_index = 0,
    *memory_tmp;

  unsigned short row_index;

  char *contig_seq,
   *contig_seq_ptr;

  Str *as_seq;

  error_check(err);

  hit_ptr = NULL;
  as_seq = str_new();
  hitcounter = strarray_size(hit_information->hit_gi);
  contig_seq_ptr = str_get(MATRIXSTRUCT(query_dna));

  /* in allen Leserahmen nach kodierenden Abschnitten suchen */
  for (row_index = 0; row_index < 7; row_index++)
  {
    if (row_index != 3)
    {
      /* Anzahl kodierender Abschnitte im aktuellen Leserahmen */
      arraysize = array_size(regionmatrix[row_index][0].from);

      for (array_index = 0; array_index < arraysize; array_index++)
      {
        /* Bereichsgrenzen des kodierenden Abschnittes - Werte bereits in
           0- bis X-Kodierung */
        from =
          *(unsigned long *) array_get(regionmatrix[row_index][0].from,
                                       array_index);
        to =
          *(unsigned long *) array_get(regionmatrix[row_index][0].to,
                                       array_index);

        /* nur wenn die Mindestlaenge in Anzahl AS erfuellt ist, erfolgt
           eine Ausgabe */
        if ((to - from + 1) / 3 > ARGUMENTSSTRUCT(min_as))
        {
          /* to-from+1+1 entspricht Anzahl der zu betrachtenden Basen +
             Stringende-Zeichen */
          contig_seq_diff = to - from + 2;

          contig_seq = ma_malloc(contig_seq_diff*sizeof (char));
          /* kopieren von contig_seq_diff-1 Zeichen */
          (void) snprintf(contig_seq, contig_seq_diff, "%s",
                          contig_seq_ptr + from);
          assert(contig_seq);

          had_err = as_coding(parsestruct_ptr,
                              contig_seq_ptr, as_seq, from, to, row_index, err);

          /* je nach Ausgabeformat die entsprechende Ausgabefunktion
             aufrufen */
          switch (ARGUMENTSSTRUCT(outputfile_format))
          {
            case 1:
              print_codingheader(parsestruct_ptr, contig_seq, as_seq);
              break;
            case 2:
              print_codingheader(parsestruct_ptr, contig_seq, as_seq);
              break;
            case 3:
              print_codingheader(parsestruct_ptr, contig_seq, as_seq);
              break;
          }

          /* Sequnezbereich wird genauer betrachtet */
          for (seq_index = from; seq_index < to + 1; seq_index++)
          {
            /* hit_numbers ist die Anzahl der beteiligten Hits an
               Sequenzposition seq_index */
            hit_numbers =
              array_size(combinedscore_matrix[row_index][seq_index].
                         hit_number);

            if (hit_ptr == NULL)
            {
              hit_ptr = ma_calloc(hitcounter, sizeof (unsigned long));
            }
            /* Abfragen aller Hit-Nummern fuer jede Position */
            for (hit_index = 0; hit_index < hit_numbers; hit_index++)
            {
              /* current entspricht der Hit-Nummer */
              current =
                *(unsigned long *)
                array_get(combinedscore_matrix[row_index][seq_index].
                          hit_number, hit_index);
              /* hit_ptr wird an der Stelle current auf 1 gesetzt - so
                 wird eine Mehrfachnennung vermieden */
              hit_ptr[current] |= 1;
            }
          }

          /* hitcounter ist max. Anzahl der Hits */
          for (seq_index = 0; seq_index < hitcounter; seq_index++)
          {

            assert(hit_ptr != NULL);
            /* Wenn an der Position seq_index eine 1 steht, kommt der
               entsprechende Hit(-Eintrag) in der Ergebnismenge vor und
               muss ausgegeben werden */
            if (hit_ptr[seq_index])
            {
              /* Die aktuelle Hit-Definition wird eingelesen */
              str_set(parsestruct_ptr->result_hits,
                      strarray_get(hit_information->hit_def, seq_index));

              hit_from =
                atoi(strarray_get(hit_information->hit_from, seq_index));
              hit_to =
                atoi(strarray_get(hit_information->hit_to, seq_index));

              /* ueberpruefen, ob der aktuelle Hit bereits erfasst wurde */
              if (!hashtable_get
                  (parsestruct_ptr->resulthits,
                   str_get(parsestruct_ptr->result_hits)))
              {
                /* Hit noch nicht erfasst, also in das Statistik-Array
                   einfuegen */
                strarray_add_cstr(HITSTRUCT(hits_statistic),
                                  str_get(parsestruct_ptr->result_hits));
                /* Position der Hit-Def. im Statistik-Array */
                string_number =
                  strarray_size(HITSTRUCT(hits_statistic)) - 1;

                /* Hitscounter: zaehlt , Hitsnumber: zaehlt das Auftreten
                   der Hit-Def. in der Ergebnismenge */
                HITSTRUCT(hitsnumber) =
                  HITSTRUCT(hitsnumber) + hit_to - hit_from + 1;

                /* Bei Bedarf speicher reallokieren */
                if (string_number + 1 > parsestruct_ptr->hits_memory - 20)
                {
                  hitsnum_tmp = HITSTRUCT(hitsnum);
                  memory_tmp = HITSTRUCT(memory);

                  array_reset(parsestruct_ptr->value_tmp);
                  (void) hashtable_foreach(parsestruct_ptr->resulthits,
                                           newmemory_hash, parsestruct_ptr,
                                           err);

                  /* Speichergroesse erhoehen */
                  parsestruct_ptr->hits_memory =
                    MEMORY_SIZE + parsestruct_ptr->hits_memory;

                  HITSTRUCT(hitsnum) =
                    ma_realloc(hitsnum_tmp,
                               parsestruct_ptr->hits_memory *
                               sizeof (unsigned long));
                  HITSTRUCT(memory) =
                    ma_realloc(memory_tmp,
                               parsestruct_ptr->hits_memory *
                               sizeof (unsigned long));

                  hashtable_reset(parsestruct_ptr->resulthits);

                  for (hash_index = 0;
                       hash_index < array_size(parsestruct_ptr->value_tmp);
                       hash_index++)
                  {
                    hashtable_add(parsestruct_ptr->resulthits,
                                  (char *)
                                  strarray_get(HITSTRUCT(hits_statistic),
                                               *(unsigned long *)
                                               array_get(parsestruct_ptr->
                                                         value_tmp,
                                                         hash_index)),
                                  HITSTRUCT(memory +
                                            *(unsigned long *)
                                            array_get(parsestruct_ptr->
                                                      value_tmp,
                                                      hash_index)));
                  }
                }

                tmp_var = string_number;

                *(HITSTRUCT(memory + tmp_var)) = tmp_var;
                *(HITSTRUCT(hitsnum + tmp_var)) = hit_to - hit_from + 1;

                hashtable_add(parsestruct_ptr->resulthits,
                              (char *)
                              strarray_get(HITSTRUCT(hits_statistic),
                                           string_number),
                              HITSTRUCT(memory + tmp_var));
              }
              else
              {
                HITSTRUCT(hitsnumber) =
                  HITSTRUCT(hitsnumber) + hit_to - hit_from + 1;
                tmp_var =
                  *(unsigned long *) hashtable_get(parsestruct_ptr->
                                                   resulthits,
                                                   str_get
                                                   (parsestruct_ptr->
                                                    result_hits));

                *(HITSTRUCT(hitsnum) + tmp_var) =
                  *(HITSTRUCT(hitsnum) + tmp_var) + hit_to - hit_from + 1;;
              }

              /* je nach Ausgabeformat die entsprechende
                 Hit-Information-Ausgabefunktion aufrufen */
              switch (ARGUMENTSSTRUCT(outputfile_format))
              {
                case 1:
                  print_hitinformation(parsestruct_ptr, hit_information,
                                       seq_index);
                  break;
                case 2:
                  print_hitinformation(parsestruct_ptr, hit_information,
                                       seq_index);
                  break;
                case 3:
                  print_hitinformation(parsestruct_ptr, hit_information,
                                       seq_index);
                  break;
              }
              str_reset(parsestruct_ptr->result_hits);
            }
          }

          switch (ARGUMENTSSTRUCT(outputfile_format))
          {
              /* Ausgabe Coding-DNA-Header - txt */
            case 1:
              genfile_xprintf(FILEPOINTEROUT, "\n\n");
              break;
              /* Ausgabe Coding-DNA-Header - html */
            case 2:
              break;
              /* Ausgabe Coding-DNA-Header - xml */
            case 3:
              genfile_xprintf(FILEPOINTEROUT, "        </Hit_infos>\n");
              genfile_xprintf(FILEPOINTEROUT, "      </Hit>\n");
              break;
          }

          ma_free(hit_ptr);
          ma_free(contig_seq);
          hit_ptr = NULL;
          str_reset(as_seq);
        }
      }
    }
  }
  str_delete(as_seq);
}

static void print_codingheader(const ParseStruct *parsestruct_ptr,
                               const char *contig_seq, Str * as_seq)
{
  assert(contig_seq);

  switch (ARGUMENTSSTRUCT(outputfile_format))
  {
      /* Ausgabe Coding-DNA-Header - txt */
    case 1:
      genfile_xprintf(FILEPOINTEROUT, "Coding-DNA: \n");
      genfile_xprintf(FILEPOINTEROUT, "%s\n", contig_seq);
      genfile_xprintf(FILEPOINTEROUT, "Protein-Seq: ");
      genfile_xprintf(FILEPOINTEROUT, "%s\n", str_get(as_seq));
      genfile_xprintf(FILEPOINTEROUT, "Hit-Information Section\n");
      break;
      /* Ausgabe Coding-DNA-Header - html */
    case 2:
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"font_header\"><br>"
                      "Coding-DNA</font></td></tr>");
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"class\">%s</font>"
                      "</td></tr>",
                      contig_seq);
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"font_header\">"
                      "Protein-Sequence</font></td></tr>");
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"class\">%s</font>"
                      "</td></tr>", str_get(as_seq));
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"font_header\">"
                      "Hit-Information Section</font></td></tr>");
      break;
      /* Ausgabe Coding-DNA-Header - xml */
    case 3:
      genfile_xprintf(FILEPOINTEROUT, "      <Hit>\n");
      genfile_xprintf(FILEPOINTEROUT, "        <Hit_dna>%s</Hit_dna>\n",
                      contig_seq);
      genfile_xprintf(FILEPOINTEROUT,
                      "        <Hit_protein-seq>%s</Hit_protein-seq>\n",
                      str_get(as_seq));
      genfile_xprintf(FILEPOINTEROUT, "        <Hit_infos>\n");
      break;
  }
}

static void print_hitinformation(const ParseStruct *parsestruct_ptr,
                                 const HitInformation *hit_information,
                                 unsigned long seq_index)
{
  /* je nach Ausgabeformat schreiben der Hit-Informationen */
  switch (ARGUMENTSSTRUCT(outputfile_format))
  {
      /* txt */
    case 1:
      genfile_xprintf(FILEPOINTEROUT, "gi-nr: gi|%s ",
                      strarray_get(hit_information->hit_gi, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "gi_def: %s ",
                      strarray_get(hit_information->hit_def, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "hsp_num: %s ",
                      strarray_get(hit_information->hit_hsp_nr,
                                   seq_index));
      genfile_xprintf(FILEPOINTEROUT, "from: %s ",
                      strarray_get(hit_information->hit_from, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "to: %s\n",
                      strarray_get(hit_information->hit_to, seq_index));
      break;
      /* html */
    case 2:
      genfile_xprintf(FILEPOINTEROUT,
                      "<tr><td colspan=\"2\"><font class=\"class\">gi-nr: ");
      genfile_xprintf(FILEPOINTEROUT,
                      "<a href=\"http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi"
                      "?db=nuccore&id=%s\">",
                      strarray_get(hit_information->hit_gi, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "gi|%s</a>  ",
                      strarray_get(hit_information->hit_gi, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "gi_def: %s ",
                      strarray_get(hit_information->hit_def, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "hsp_num: %s ",
                      strarray_get(hit_information->hit_hsp_nr,
                                   seq_index));
      genfile_xprintf(FILEPOINTEROUT, "from: %s ",
                      strarray_get(hit_information->hit_from, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "to: %s</font></td></tr>",
                      strarray_get(hit_information->hit_to, seq_index));
      break;
      /* xml */
    case 3:
      genfile_xprintf(FILEPOINTEROUT, "          <Infos>\n");
      genfile_xprintf(FILEPOINTEROUT,
                      "            <Infos_gi-nr>gi|%s</Infos_gi-nr>\n",
                      strarray_get(hit_information->hit_gi, seq_index));
      genfile_xprintf(FILEPOINTEROUT,
                      "            <Infos_gi-def>%s</Infos_gi-def>\n",
                      strarray_get(hit_information->hit_def, seq_index));
      genfile_xprintf(FILEPOINTEROUT,
                      "            <Infos_hsp-num>%s</Infos_hsp-num>\n",
                      strarray_get(hit_information->hit_hsp_nr,
                                   seq_index));
      genfile_xprintf(FILEPOINTEROUT,
                      "            <Infos_from>%s</Infos_from>\n",
                      strarray_get(hit_information->hit_from, seq_index));
      genfile_xprintf(FILEPOINTEROUT,
                      "            <Infos_to>%s</Infos_to>\n",
                      strarray_get(hit_information->hit_to, seq_index));
      genfile_xprintf(FILEPOINTEROUT, "          </Infos>\n");
      break;
  }
}

static void output_statistics_txt(const ParseStruct *parsestruct_ptr)
{
  /* schreiben des Query-DNA Headers inkl. Query-Def. und Query-Sequenz */
  genfile_xprintf(FILEPOINTEROUT, "%-8.4f   ",
                  ((double) *(HITSTRUCT(hitsnum) + HITSTRUCT(stat_pos)) /
                   (double) HITSTRUCT(hitsnumber)) * 100);
  genfile_xprintf(FILEPOINTEROUT, "%s\n",
                  strarray_get(HITSTRUCT(hits_statistic),
                               *(HITSTRUCT(memory) +
                                 HITSTRUCT(stat_pos))));
}

static void output_statistics_html(const ParseStruct *parsestruct_ptr)
{
  genfile_xprintf(FILEPOINTEROUT,
                  "<tr><td align=\"right\" width=\"50\">%-8.4f </td>",
                  ((float) *(HITSTRUCT(hitsnum) + HITSTRUCT(stat_pos)) /
                   (float) HITSTRUCT(hitsnumber)) * 100);
  genfile_xprintf(FILEPOINTEROUT,
                  "<td><font class=\"font_header\"> %s</font></td></tr>",
                  strarray_get(HITSTRUCT(hits_statistic),
                               *(HITSTRUCT(memory) +
                                 HITSTRUCT(stat_pos))));
}

static void output_statistics_xml(const ParseStruct *parsestruct_ptr)
{
  genfile_xprintf(FILEPOINTEROUT, "    <Statistics>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Statistics_percent>%-8.4f</Statistics_percent>\n",
                  ((double) *(HITSTRUCT(hitsnum) + HITSTRUCT(stat_pos)) /
                   (double) HITSTRUCT(hitsnumber)) * 100);
  genfile_xprintf(FILEPOINTEROUT,
                  "      <Statistics_gi-def>%s</Statistics_gi-def>\n",
                  strarray_get(HITSTRUCT(hits_statistic),
                               *(HITSTRUCT(memory) +
                                 HITSTRUCT(stat_pos))));
  genfile_xprintf(FILEPOINTEROUT, "    </Statistics>\n");
}

static void output_footer_html(const ParseStruct *parsestruct_ptr)
{
  /* Schreiben des HTML-Footers */
  genfile_xprintf(FILEPOINTEROUT, "</td></tr></table>");
  genfile_xprintf(FILEPOINTEROUT, "</table>\n</body>\n</html>");
}

static void output_footer_xml(const ParseStruct *parsestruct_ptr)
{
  /* Schreiben des HTML-Footers */
  genfile_xprintf(FILEPOINTEROUT,
                  "  </MetagenomethreaderOutput_statistics>\n");
  genfile_xprintf(FILEPOINTEROUT, "</MetagenomethreaderOutput>\n");
}

static int as_coding(const ParseStruct *parsestruct_ptr,
                     char *contig_seq,
                      Str *as_seq,
                      unsigned long from,
                      unsigned long to,
                      unsigned short current_row, Error * err)
{
  int had_err = 0;

  unsigned long startpoint = from,
    endpoint = to,
    startpoint_start,
    startpoint_atg,
    startpoint_safe,
    contig_len;

  char contig_triplet[3],
   *contig_seq_tri = NULL,
    contig_as;

  short current_frame,
    current_frame_tmp,
    found = 0,
    found_start = 0,
    found_end = 0,
    start_codon = 0;

  error_check(err);
  assert(contig_seq);

  contig_len = strlen(contig_seq);
  current_frame = get_current_frame(current_row);
  current_frame_tmp = current_frame;

  if (current_frame_tmp < 0)
  {
    current_frame_tmp *= -1;

    had_err = mg_reverse_complement(contig_seq, contig_len, err);

    startpoint = contig_len - 1 - to;
    endpoint = contig_len - from;
  }

  if (!had_err)
  {
    if (startpoint < 3)
    {
      startpoint = current_frame_tmp - 1;
      startpoint_start = startpoint;
    }
    else
    {
      startpoint -= (((startpoint) - current_frame_tmp) % 3);
      startpoint -= 1;
      startpoint_start = startpoint;
    }

    startpoint_safe = startpoint;

    while ((startpoint <= endpoint) && (startpoint <= contig_len - 3))
    {
      /* Aktuelles Query-Triplet */
      contig_triplet[0] = contig_seq[startpoint];
      contig_triplet[1] = contig_seq[startpoint + 1];
      contig_triplet[2] = contig_seq[startpoint + 2];
      /* XXX assert(contig_triplet); */

      /* Bestimmen der AS der jeweiligen Triplets */
      contig_as = mg_codon2amino(contig_triplet[0],
                                 contig_triplet[1], contig_triplet[2]);
      assert(contig_as);

      str_append_char(as_seq, contig_as);
      startpoint += 3;
    }

    if (ARGUMENTSSTRUCT(extended_mode))
    {
      Str *as_seq_start;

      as_seq_start = str_new();

      contig_seq_tri = ma_malloc(4*sizeof (char));

      /* DNA-Basen-Triplet einlesen */
      contig_seq_tri[0] = tolower(contig_seq[startpoint - 3]);
      contig_seq_tri[1] = tolower(contig_seq[startpoint - 2]);
      contig_seq_tri[2] = tolower(contig_seq[startpoint - 1]);
      contig_seq_tri[3] = '\0';
      assert(contig_seq_tri);

      found = check_stopcodon(contig_seq_tri);

      while (startpoint <= contig_len - 3 && !found_end && found)
      {
        /* DNA-Basen-Triplet einlesen */
        contig_seq_tri[0] = tolower(contig_seq[startpoint - 3]);
        contig_seq_tri[1] = tolower(contig_seq[startpoint - 2]);
        contig_seq_tri[2] = tolower(contig_seq[startpoint - 1]);
        contig_seq_tri[3] = '\0';
        assert(contig_seq_tri);

        found_end = check_stopcodon(contig_seq_tri);

        if (found_end)
        {
          /* Aktuelles Query-Triplet */
          contig_triplet[0] = contig_seq[startpoint];
          contig_triplet[1] = contig_seq[startpoint + 1];
          contig_triplet[2] = contig_seq[startpoint + 2];
          /* XXX assert(contig_triplet); */

          /* Bestimmen der AS der jeweiligen Triplets */
          contig_as = mg_codon2amino(contig_triplet[0],
                                     contig_triplet[1], contig_triplet[2]);
          assert(contig_as);

          str_append_char(as_seq, contig_as);
        }
        /* Startwert um 3 Basen weitersetzen */
        startpoint += 3;
      }

      /* DNA-Basen-Triplet einlesen */
      contig_seq_tri[0] = tolower(contig_seq[startpoint_start]);
      contig_seq_tri[1] = tolower(contig_seq[startpoint_start + 1]);
      contig_seq_tri[2] = tolower(contig_seq[startpoint_start + 2]);
      contig_seq_tri[3] = '\0';
      assert(contig_seq_tri);

      start_codon = check_startcodon(parsestruct_ptr, contig_seq_tri);

      found = 0;

      if (!start_codon)
      {
        while (startpoint_start > 2 && !found)
        {
          /* DNA-Basen-Triplet einlesen */
          contig_seq_tri[0] = tolower(contig_seq[startpoint_start - 3]);
          contig_seq_tri[1] =
            tolower(contig_seq[startpoint_start - 2]);
          contig_seq_tri[2] =
            tolower(contig_seq[startpoint_start - 1]);
          contig_seq_tri[3] = '\0';
          assert(contig_seq_tri);

          found = check_stopcodon(contig_seq_tri);

          startpoint_atg = startpoint_start;

          if (found || startpoint_start < 3)
          {
            while (startpoint_atg <= startpoint_safe - 2)
            {
              if (!found_start)
              {
                /* DNA-Basen-Triplet einlesen */
                contig_seq_tri[0] = tolower(contig_seq[startpoint_atg]);
                contig_seq_tri[1] =
                  tolower(contig_seq[startpoint_atg + 1]);
                contig_seq_tri[2] =
                  tolower(contig_seq[startpoint_atg + 2]);
                contig_seq_tri[3] = '\0';
                assert(contig_seq_tri);

                start_codon = check_startcodon(parsestruct_ptr,
                                               contig_seq_tri);

                /* ueberpruefen auf Start-Codon */
                if (start_codon)
                {
                  str_append_char(as_seq_start, 'M');
                  found_start = 1;
                }
              }
              else
              {
                /* Aktuelles Query-Triplet */
                contig_triplet[0] = contig_seq[startpoint_atg];
                contig_triplet[1] = contig_seq[startpoint_atg + 1];
                contig_triplet[2] = contig_seq[startpoint_atg + 2];
                /* XXX assert(contig_triplet); */

                /* Bestimmen der AS der jeweiligen Triplets */
                contig_as = mg_codon2amino(contig_triplet[0],
                                           contig_triplet[1],
                                           contig_triplet[2]);
                assert(contig_as);

                str_append_char(as_seq_start, contig_as);
              }
              /* Startwert um 3 Basen weitersetzen */
              startpoint_atg += 3;
            }
          }
          startpoint_start -= 3;
        }
        str_append_str(as_seq_start, as_seq);
        str_reset(as_seq);
        str_append_str(as_seq, as_seq_start);
      }
      str_reset(as_seq_start);
      str_delete(as_seq_start);
      ma_free(contig_seq_tri);
    }

    if (current_frame < 0)
    {
      had_err = mg_reverse_complement(contig_seq, contig_len, err);
    }
  }
  return had_err;
}

static int newmemory_hash(UNUSED void *key,
                          void *value, void *data,
                          UNUSED Error * err)
{
  /* Parsestruct-Struktur */
  ParseStruct *parsestruct_ptr = (ParseStruct *) data;

  /* Position des aktuell betrachteten Schluessels */
  HITSTRUCT(stat_pos) = *(unsigned long *) value;

  error_check(err);

  array_add(parsestruct_ptr->value_tmp, HITSTRUCT(stat_pos));

  return 0;
}

static void output_statistics_header(const ParseStruct *parsestruct_ptr)
{
  /* Statistik-Bereich-Header schreiben */
  if (ARGUMENTSSTRUCT(outputfile_format) == 3)
  {
    genfile_xprintf(FILEPOINTEROUT,
                    "  <MetagenomethreaderOutput_statistics>\n");
  }
  else if (ARGUMENTSSTRUCT(outputfile_format) == 2)
  {
    genfile_xprintf(FILEPOINTEROUT,
                    "<tr><td colspan=\"2\"><table cellspacing=\"1\" "
                    "cellpadding=\"3\">");
    genfile_xprintf(FILEPOINTEROUT,
                    "<tr><td colspan=\"2\"><font class=\"font_header\"><br><br>"
                    "Statistic-Section<br><br></font></td></tr>");
  }
  else
  {
    genfile_xprintf(FILEPOINTEROUT, "Statistic-Section\n\n");
  }
}

static short check_startcodon(const ParseStruct *parsestruct_ptr,
                              const char *contig_seq_ptrfct)
{
  unsigned short codon_status = 0;

  /* jeder if-Zweig ueberprueft entsprechend des
     Metagenomethreader-Arguments zur Verwendung alternativer Start-Codons
     auf Start-Codons */
  if (ARGUMENTSSTRUCT(codon_mode) == 2)
  {
    if (!strcmp(contig_seq_ptrfct, "atg")
        || !strcmp(contig_seq_ptrfct, "ctg")
        || !strcmp(contig_seq_ptrfct, "gtg")
        || !strcmp(contig_seq_ptrfct, "aug")
        || !strcmp(contig_seq_ptrfct, "cug")
        || !strcmp(contig_seq_ptrfct, "gug"))
    {
      codon_status = 1;
    }
  }
  else if (ARGUMENTSSTRUCT(codon_mode) == 3)
  {
    if (!strcmp(contig_seq_ptrfct, "atg")
        || !strcmp(contig_seq_ptrfct, "ctg")
        || !strcmp(contig_seq_ptrfct, "gtg")
        || !strcmp(contig_seq_ptrfct, "ttg")
        || !strcmp(contig_seq_ptrfct, "aug")
        || !strcmp(contig_seq_ptrfct, "cug")
        || !strcmp(contig_seq_ptrfct, "gug")
        || !strcmp(contig_seq_ptrfct, "uug"))
    {
      codon_status = 1;
    }
  }
  else
  {
    if (!strcmp(contig_seq_ptrfct, "atg")
        || !strcmp(contig_seq_ptrfct, "aug"))
    {
      codon_status = 1;
    }
  }

  return codon_status;
}

static void output_close_iteration_xml(const ParseStruct *parsestruct_ptr)
{
  genfile_xprintf(FILEPOINTEROUT, "    </Iteration_hits>\n");
  genfile_xprintf(FILEPOINTEROUT, "  </Iteration>\n");
  genfile_xprintf(FILEPOINTEROUT,
                  "  </MetagenomethreaderOutput_iterations>\n");
}
