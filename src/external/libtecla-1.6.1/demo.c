/*
 * Copyright (c) 2000, 2001, 2002, 2003, 2004 by Martin C. Shepherd
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

/* The function which displays the introductory text of the demo */

static void show_demo_introduction(GetLine *gl);

/*.......................................................................
 * This program demonstrates how to use gl_get_line() as a line editor to
 * to enable users to enter input. It takes no arguments.
 */
int main(int argc, char *argv[])
{
  char *line;             /* A line of input */
  GetLine *gl;            /* The line editor */
  int major,minor,micro;  /* The version number of the library */
/*
 * Create the line editor, specifying a max line length of 500 bytes,
 * and 10000 bytes to allocate to storage of historical input lines.
 */
  gl = new_GetLine(500, 5000);
  if(!gl)
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
  printf("\n Welcome to the main demo program of libtecla version %d.%d.%d\n",
	 major, minor, micro);
/*
 * Display an introductory banner.
 */
  show_demo_introduction(gl);
/*
 * Load history.
 */
#ifndef WITHOUT_FILE_SYSTEM
  (void) gl_load_history(gl, "~/.demo_history", "#");
#endif
/*
 * Read lines of input from the user and print them to stdout.
 */
  do {
/*
 * Get a new line from the user.
 */
    line = gl_get_line(gl, "$ ", NULL, 0);
    if(!line)
      break;
/*
 * Display what was entered.
 */
    if(printf("You entered: %s", line) < 0 || fflush(stdout))
      break;
/*
 * If the user types "exit", quit the program.
 */
    if(strcmp(line, "exit\n")==0)
      break;
    else if(strcmp(line, "history\n")==0)
      gl_show_history(gl, stdout, "%N  %T   %H\n", 0, -1);
    else if(strcmp(line, "size\n")==0) {
      GlTerminalSize size = gl_terminal_size(gl, 80, 24);
      printf("Terminal size = %d columns x %d lines.\n", size.ncolumn,
	     size.nline);
    } else if(strcmp(line, "clear\n")==0) {
      if(gl_erase_terminal(gl))
	return 1;
    };
  } while(1);
/*
 * Save historical command lines.
 */
#ifndef WITHOUT_FILE_SYSTEM
  (void) gl_save_history(gl, "~/.demo_history", "#", -1);
#endif
/*
 * Clean up.
 */
  gl = del_GetLine(gl);
  return 0;
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
    "This program is a simple shell with which you can experiment with the ",
    "line editing and tab completion facilities provided by the gl_get_line() ",
    "function. The file demo.c also serves as a fully commented example ",
    "of how to use gl_get_line().\n"
  };
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
