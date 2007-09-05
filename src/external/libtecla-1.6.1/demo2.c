/*
 * Copyright (c) 2000, 2001, 2002, 2003, 2004 by Martin C. Shepherd.
 * 
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, and/or sell copies of the Software, and to permit persons
 * to whom the Software is furnished to do so, provided that the above
 * copyright notice(s) and this permission notice appear in all copies of
 * the Software and that both the above copyright notice(s) and this
 * permission notice appear in supporting documentation.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT
 * OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL
 * INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING
 * FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
 * WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 * 
 * Except as contained in this notice, the name of a copyright holder
 * shall not be used in advertising or otherwise to promote the sale, use
 * or other dealings in this Software without prior written authorization
 * of the copyright holder.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <locale.h>

#include <unistd.h>
#include <sys/stat.h>

#include "libtecla.h"

/*
 * If the library is being built with file-system access excluded, this
 * demo program won't have anything to demonstrate.
 */
#ifdef WITHOUT_FILE_SYSTEM
int main(int argc, char *argv[])
{
  fprintf(stderr, "\n"
	  "  This program normally demonstrates tecla's path-lookup\n"
	  "  facility. However libtecla has been installed with\n"
	  "  file-system facilities explicitly excluded, so there is\n"
          "  nothing to demonstrate.\n\n");
  return 1;
}
#else
/*
 * Encapsulate the resources needed by this demo.
 */
typedef struct {
  GetLine *gl;      /* The line editor */
  PathCache *pc;    /* A cache of executables in the user's path */
  PcaPathConf *ppc; /* The configuration argument of pca_path_completions() */
} DemoRes;

/*
 * The following functions allocate and free instances of the above
 * structure.
 */
static DemoRes *new_DemoRes(void);
static DemoRes *del_DemoRes(DemoRes *res);

/*
 * Search backwards for the start of a pathname.
 */
static char *start_of_path(const char *string, int back_from);

/*
 * Find the array indexes of the first character of the first
 * space-delimited word in the specified string, and of the character
 * that follows it.
 */
static int get_word_limits(const char *string, int *wa, int *wb);

/*
 * This is the demonstration completion callback function (defined below).
 */
static CPL_MATCH_FN(demo_cpl_fn);

/* The function which displays the introductory text of the demo */

static void show_demo_introduction(GetLine *gl);

/*.......................................................................
 * This demo takes no arguments. It reads lines of input until the
 * word 'exit' is entered, or C-d is pressed. It replaces the default
 * tab-completion callback with one which when invoked at the start of
 * a line, looks up completions of commands in the user's execution
 * path, and when invoked in other parts of the line, reverts to
 * normal filename completion. Whenever a new line is entered, it
 * extracts the first word on the line, looks it up in the user's
 * execution path to see if it corresponds to a known executable file,
 * and if so, displays the full pathname of the file, along with the
 * remaining arguments.
 */
int main(int argc, char *argv[])
{
  char *line;             /* A line of input */
  DemoRes *res;           /* The resources of the demo */
  int wa,wb;              /* The delimiting indexes of a word in line[] */
  int major,minor,micro;  /* The version number of the library */
/*
 * Allocate the resources needed by this demo.
 */
  res = new_DemoRes();
  if(!res)
    return 1;
/*
 * If the user has the LC_CTYPE or LC_ALL environment variables set,
 * enable display of characters corresponding to the specified locale.
 */
  (void) setlocale(LC_CTYPE, "");
/*
 * Lookup and display the version number of the library.
 */
  libtecla_version(&major, &minor, &micro);
  printf("\n Welcome to the path-search demo of libtecla version %d.%d.%d\n",
	 major, minor, micro);
/*
 * Display some introductory text, left-justifying it within the current
 * width of the terminal and enclosing it in a box of asterixes.
 */
  show_demo_introduction(res->gl);
/*
 * Read lines of input from the user and print them to stdout.
 */
  do {
/*
 * Get a new line from the user.
 */
    line = gl_get_line(res->gl, "$ ", NULL, 0);
    if(!line)
      break;
/*
 * Work out the extent of the first word in the input line, and
 * try to identify this as a command in the path, displaying the
 * full pathname of the match if found.
 */
    if(get_word_limits(line, &wa, &wb) == 0) {
      char *cmd = pca_lookup_file(res->pc, line + wa, wb-wa, 0);
      if(cmd) {
	printf("Command=%s\n", cmd);
	printf("Arguments=%s", line+wb);
      } else {
	printf("Command not found\n");
      };
    };
/*
 * If the user types "exit", quit the program.
 */
    if(strcmp(line, "exit\n")==0)
      break;
  } while(1);
/*
 * Clean up.
 */
  res = del_DemoRes(res);
  return 0;
}

/*.......................................................................
 * This completion callback searches for completions of executables in
 * the user's path when invoked on a word at the start of the path, and
 * performs normal filename completion elsewhere.
 */
static CPL_MATCH_FN(demo_cpl_fn)
{
/*
 * Get the resource object that was passed to gl_customize_completion().
 */
  DemoRes *res = (DemoRes *) data;
/*
 * Find the start of the filename prefix to be completed, searching
 * backwards for the first unescaped space, or the start of the line.
 */
  char *start = start_of_path(line, word_end);
/*
 * Skip spaces preceding the start of the prefix.
 */
  while(start > line && isspace((int)(unsigned char) start[-1]))
    start--;
/*
 * If the filename prefix is at the start of the line, attempt
 * to complete the filename as a command in the path. Otherwise
 * perform normal filename completion.
 */
  return (start == line) ?
    pca_path_completions(cpl, res->ppc, line, word_end) :
    cpl_file_completions(cpl, NULL, line, word_end);
}

/*.......................................................................
 * Search backwards for the potential start of a filename. This
 * looks backwards from the specified index in a given string,
 * stopping at the first unescaped space or the start of the line.
 *
 * Input:
 *  string  const char *  The string to search backwards in.
 *  back_from      int    The index of the first character in string[]
 *                        that follows the pathname.
 * Output:
 *  return        char *  The pointer to the first character of
 *                        the potential pathname, or NULL on error.
 */
static char *start_of_path(const char *string, int back_from)
{
  int i, j;
/*
 * Search backwards from the specified index.
 */
  for(i=back_from-1; i>=0; i--) {
    int c = string[i];
/*
 * Stop on unescaped spaces.
 */
    if(isspace((int)(unsigned char)c)) {
/*
 * The space can't be escaped if we are at the start of the line.
 */
      if(i==0)
        break;
/*
 * Find the extent of the escape characters which precedes the space.
 */
      for(j=i-1; j>=0 && string[j]=='\\'; j--)
	;
/*
 * If there isn't an odd number of escape characters before the space,
 * then the space isn't escaped.
 */
      if((i - 1 - j) % 2 == 0)
	break;
    };
  };
  return (char *)string + i + 1;
}

/*.......................................................................
 * Create a new DemoRes object containing the resources needed by the
 * demo.
 *
 * Output:
 *  return  DemoRes *  The new object, or NULL on error.
 */
static DemoRes *new_DemoRes(void)
{
  DemoRes *res;  /* The object to be returned */
/*
 * Allocate the container.
 */
  res = (DemoRes *)malloc(sizeof(DemoRes));
  if(!res) {
    fprintf(stderr, "new_DemoRes: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_DemoRes().
 */
  res->gl = NULL;
  res->pc = NULL;
  res->ppc = NULL;
/*
 * Create the line editor, specifying a max line length of 500 bytes,
 * and 10000 bytes to allocate to storage of historical input lines.
 */
  res->gl = new_GetLine(500, 10000);
  if(!res->gl)
    return del_DemoRes(res);
/*
 * Enable text attribute formatting directives in prompt strings.
 */
  gl_prompt_style(res->gl, GL_FORMAT_PROMPT);
/*
 * Allocate a cache of the executable files found in the user's path.
 */
  res->pc = new_PathCache();
  if(!res->pc)
    return del_DemoRes(res);
/*
 * Populate the cache with the contents of the user's path.
 */
  if(pca_scan_path(res->pc, getenv("PATH")))
    return del_DemoRes(res);
/*
 * Arrange for susequent calls to pca_lookup_file() and pca_path_completions()
 * to only report files that are executable by the user.
 */
  pca_set_check_fn(res->pc, cpl_check_exe, NULL);
/*
 * Allocate a configuration object for use with pca_path_completions().
 */
  res->ppc = new_PcaPathConf(res->pc);
  if(!res->ppc)
    return del_DemoRes(res);
/*
 * Replace the builtin filename completion callback with one which
 * searches for completions of executables in the user's path when
 * invoked on a word at the start of the path, and completes files
 * elsewhere.
 */
  if(gl_customize_completion(res->gl, res, demo_cpl_fn))
    return del_DemoRes(res);
  return res;
}

/*.......................................................................
 * Delete a DemoRes object.
 *
 * Input:
 *  res     DemoRes *  The object to be deleted.
 * Output:
 *  return  DemoRes *  The deleted object (always NULL).
 */
static DemoRes *del_DemoRes(DemoRes *res)
{
  if(res) {
    res->gl = del_GetLine(res->gl);
    res->pc = del_PathCache(res->pc);
    res->ppc = del_PcaPathConf(res->ppc);
    free(res);
  };
  return NULL;
}

/*.......................................................................
 * Return the limits of the word at the start of a given string, ignoring
 * leading white-space, and interpretting the first unescaped space, tab or
 * the end of the line, as the end of the word.
 *
 * Input:
 *  string   const char *  The string to tokenize.
 * Input/Output:
 *  wa,wb           int *  The indexes of the first character of the word,
 *                         and the character which follows the last
 *                         character of the word, will be assigned to
 *                         *wa and *wb, respectively.
 * Output:
 *  return          int    0 - A word was found.
 *                         1 - No word was found before the end of the
 *                             string.
 */
static int get_word_limits(const char *string, int *wa, int *wb)
{
  int escaped = 0;  /* True if the next character is escaped */
/*
 * Skip leading white-space.
 */
  for(*wa=0; isspace((int)(unsigned char)string[*wa]); (*wa)++)
    ;
/*
 * Find the first unescaped space, stopping early if the end of the
 * string is reached.
 */
  for(*wb = *wa; ; (*wb)++) {
    int c = string[*wb];
    if(c=='\\')
      escaped = !escaped;
    else if((!escaped && isspace((int)(unsigned char)c)) || c=='\0')
      break;
  };
  return *wa == *wb;
}

/*.......................................................................
 * Display introductory text to the user, formatted according to the
 * current terminal width and enclosed in a box of asterixes.
 *
 * Input:
 *  gl      GetLine *   The resource object of gl_get_line().
 */
static void show_demo_introduction(GetLine *gl)
{
  int start;     /* The column in which gl_display_text() left the cursor */
  int i;
/*
 * Break the indtroductory text into an array of strings, so as to
 * avoid overflowing any compiler string limits.
 */
  const char *doc[] = {
    "This program demonstrates the use of the pca_lookup_file() function ",
    "for finding executables in the UNIX PATH. It also demonstrates ",
    "tab completion of the names of executables found in the path. For ",
    "example, if you type:\n\n ta\n\nthen hit the tab key, you will be ",
    "presented with a list of executables such as tar and tail whose names ",
    "start with the string \"ta\". If you decide to add an \"r\" to select ",
    "the tar command, then you type return, the full pathname of the tar ",
    "program will be printed.\n\nThe file demo2.c contains the code ",
    "of this program, and is fully commented to enable its use as ",
    "a working example of how to use the facilities documented in the ",
    "pca_lookup_file man page.\n"};
/*
 * Form the top line of the documentation box by filling the area of
 * the line between a " *" prefix and a "* " suffix with asterixes.
 */
  printf("\n");
  gl_display_text(gl, 0, " *", "* ", '*', 80, 0, "\n");
/*
 * Justify the documentation text within margins of asterixes.
 */
  for(start=0,i=0; i<sizeof(doc)/sizeof(doc[0]) && start >= 0; i++)
    start = gl_display_text(gl, 0, " * ", " * ", ' ', 80, start,doc[i]);
/*
 * Draw the bottom line of the documentation box.
 */
  gl_display_text(gl, 0, " *", "* ", '*', 80, 0, "\n");
  printf("\n");
}

#endif  /* ifndef WITHOUT_FILE_SYSTEM */
