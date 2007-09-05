/*
 * Copyright (c) 2002, 2003, 2004 by Martin C. Shepherd
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
#include <setjmp.h>

#ifdef HAVE_SELECT
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#endif

#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>    
#include <sys/types.h>
#include <signal.h>

#include "libtecla.h"

/*
 * The SignalActions object provides a way to temporarily install
 * a signal handler to a given set of signals, and later restore all
 * of the signal handlers that this displaced.
 */
typedef struct {
  int nsignal;               /* The number of signals on the host OS */
  sigset_t mask;             /* The set of signals who's signal handlers */
                             /*  are stored in the following actions[] */
                             /*  array. */
  struct sigaction *actions; /* An array of nsignal actions */
} SignalActions;

static SignalActions *new_SignalActions(void);
static SignalActions *del_SignalActions(SignalActions *si);
static int displace_signal_handlers(SignalActions *si, sigset_t *mask,
				    void (*handler)(int));
static int reinstate_signal_handlers(SignalActions *si);

/* Return resources, restore the terminal to a usable state and exit */

static void cleanup_and_exit(GetLine *gl, SignalActions *si, int status);

/* The function which displays the introductory text of the demo */

static void show_demo_introduction(GetLine *gl);

/* A signal-aware version of select() */

static int demo_sigselect(int n, fd_set *readfds, fd_set *writefds,
			  fd_set *exceptfds, struct timeval *timeout,
			  sigset_t *mask, SignalActions *si);

/*
 * The following variables are accessed from signal handlers.  Note
 * that these variables don't need to be either volatile or
 * sig_atomic_t because:
 *
 * 1. Outside of signal handlers we only access them when signal
 *    delivery is blocked, so we know that no signal handlers can
 *    be accessing them at that time.
 *
 * 2. When the signal handlers that set these variables are installed,
 *    the sa_mask member of the sigaction structure is used to ensure
 *    that only one instance of these signal handlers can be running
 *    at a time, so we also know that there can't be simultaneous
 *    accesses to them by multiple signal handlers.
 */
static GetLine *demo_gl;                /* The line editor object */
static sigjmp_buf demo_setjmp_buffer;   /* The sigsetjmp() buffer */
static int demo_setjmp_signo = -1;      /* The signal that was caught */

/* Signal handlers */

static void demo_signal_handler(int signo);
static void demo_setjmp_handler(int signo);

/*
 * Set the amount of time that gl_get_line() should wait for I/O before
 * returning to let the external event loop continue.
 */
#define DEMO_IO_TIMEOUT 100000000 /* ns => 100ms */

/* The timeout handler */

static GL_TIMEOUT_FN(demo_timeout_fn);

/*.......................................................................
 * This program demonstrates the use of gl_get_line() from an external
 * event loop. It takes no arguments.
 */
int main(int argc, char *argv[])
{
  int major,minor,micro;      /* The version number of the library */
  GetLine *gl=NULL;           /* The resource object of gl_get_line() */
  SignalActions *si=NULL;     /* Temporary storage of displaced signal */
                              /*  handlers. */
  sigset_t all_signal_mask;   /* The set of signals known by gl_get_line() */
/*
 * This program requires select().
 */
#if !defined(HAVE_SELECT)
  fprintf(stderr, "The select() system call isn't available - aborting.\n");
  exit(1);
#else
/*
 * Create the line editor, specifying a maximum line length of 500 bytes,
 * and 10000 bytes to allocate to storage of historical input lines.
 */
  gl = demo_gl = new_GetLine(500, 5000);
  if(!gl)
    cleanup_and_exit(gl, si, 1);

/*
 * Allocate an object in which to temporarily record displaced
 * signal handlers.
 */
  si = new_SignalActions();
  if(!si)
    cleanup_and_exit(gl, si, 1);
/*
 * If the user has the LC_CTYPE or LC_ALL environment variables set,
 * enable display of characters corresponding to the specified locale.
 */
  (void) setlocale(LC_CTYPE, "");
/*
 * Lookup and display the version number of the library.
 */
  libtecla_version(&major, &minor, &micro);
  printf(
    "\n Welcome to the server-mode demo program of libtecla version %d.%d.%d\n",
      major, minor, micro);
/*
 * Display some introductory text, left-justifying it within the current
 * width of the terminal and enclosing it in a box of asterixes.
 */
  show_demo_introduction(gl);
/*
 * Load history.
 */
#ifndef WITHOUT_FILE_SYSTEM
  (void) gl_load_history(gl, "~/.demo_history", "#");
#endif
/*
 * In this demo, rather than having gl_get_line() return immediately
 * when it would otherwise have to wait for I/O, we register a timeout
 * callback which causes gl_get_line() to give up waiting after a short
 * interval.
 */
  gl_inactivity_timeout(gl, demo_timeout_fn, NULL, 0, DEMO_IO_TIMEOUT);
/*
 * Install our signal handlers for process termination, suspension and
 * terminal resize signals. Ignore process continuation signals.
 */
  gl_tty_signals(demo_signal_handler, demo_signal_handler, SIG_DFL,
		 demo_signal_handler);
/*
 * Get a list of all of the signals that gl_get_line() currently catches.
 */
  gl_list_signals(gl, &all_signal_mask);
/*
 * Switch gl_get_line() to non-blocking server mode.
 */
  if(gl_io_mode(gl, GL_SERVER_MODE))
    cleanup_and_exit(gl, si, 1);
/*
 * Instruct gl_get_line() to unblock any signals that it catches
 * while waiting for input. Note that in non-blocking server mode,
 * this is only necessary when using gl_inactivity_timeout() to make
 * gl_get_line() block for a non-zero amount of time.
 */
  gl_catch_blocked(gl);
/*
 * Enter the event loop.
 */
  while(1) {
    int nready;        /* The number of file-descriptors that are */
                       /*  ready for I/O */
    fd_set rfds;       /* The set of file descriptors to watch for */
                       /*  readability */
    fd_set wfds;       /* The set of file descriptors to watch for */
                       /*  writability */
/*
 * Construct the sets of file descriptors to be watched by select(),
 * starting from empty sets.
 */
    FD_ZERO(&rfds);
    FD_ZERO(&wfds);
/*
 * To ensure that no signals are received whos handlers might change
 * the requirements for the contents of the above signal sets, block
 * all of the signals that we are handling.
 */
    sigprocmask(SIG_BLOCK, &all_signal_mask, NULL);
/*
 * Depending on which direction of I/O gl_get_line()s is currently
 * waiting for, add the terminal file descriptor to either the set
 * of file descriptors to watch for readability, or those to watch
 * for writability. Note that at the start of a new line, such as
 * after an error, or the return of a completed line, we need to
 * wait for writability, so that a prompt can be written.
 */
    switch(gl_pending_io(gl)) {
    case GLP_READ:
      FD_SET(STDIN_FILENO, &rfds);
      break;
    default:
      FD_SET(STDIN_FILENO, &wfds);
      break;
    };
/*
 * Wait for I/O to become possible on the selected file descriptors.
 * The following is a signal-aware wrapper around the select() system
 * call. This wrapper guarantees that if any of the signals marked in
 * all_signal_mask arrive after the statement above where we blocked
 * these signals, it will detect this and abort with nready=-1 and
 * errno=EINTR. If instead, we just unblocked the above signals just
 * before calling a normal call to select(), there would be a small
 * window of time between those two statements in which a signal could
 * arrive without aborting select(). This would be a problem, since
 * the functions called by our signal handler may change the type
 * of I/O that gl_get_line() wants us to wait for in select().
 */
    nready = demo_sigselect(STDIN_FILENO + 1, &rfds, &wfds, NULL, NULL,
			    &all_signal_mask, si);
/*
 * We can now unblock our signals again.
 */
    sigprocmask(SIG_UNBLOCK, &all_signal_mask, NULL);
/*
 * Did an I/O error occur?
 */
    if(nready < 0 && errno != EINTR)
      cleanup_and_exit(gl, si, 1);
/*
 * If the terminal file descriptor is now ready for I/O, call
 * gl_get_line() to continue editing the current input line.
 */
    if(FD_ISSET(STDIN_FILENO, &rfds) || FD_ISSET(STDIN_FILENO, &wfds)) {
/*
 * Start or continue editing an input line.
 */
      char *line = gl_get_line(gl, "$ ", NULL, 0);
/*
 * Did the user finish entering a new line?
 */
      if(line) {
/*
 * Before writing messages to the terminal, start a new line and
 * switch back to normal terminal I/O.
 */
	gl_normal_io(gl);
/*
 * Display what was entered.
 */
	if(printf("You entered: %s", line) < 0 || fflush(stdout))
	  break;
/*
 * Implement a few simple commands.
 */
	if(strcmp(line, "exit\n")==0)
	  cleanup_and_exit(gl, si, 0);
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
/*
 * To resume command-line editing, return the terminal to raw,
 * non-blocking I/O mode.
 */
	gl_raw_io(gl);
/*
 * If gl_get_line() returned NULL because of an error or end-of-file,
 * abort the program.
 */
      } else if(gl_return_status(gl) == GLR_ERROR ||
		gl_return_status(gl) == GLR_EOF) {
	cleanup_and_exit(gl, si, 1);
      };
    };
  };
#endif
  return 0;
}

/*.......................................................................
 * This function is called to return resources to the system and restore
 * the terminal to its original state before exiting the process.
 *
 * Input:
 *  gl        GetLine *   The line editor.
 *  si  SignalActions *   The repository for displaced signal handlers.
 *  status        int     The exit code of the process.
 */
static void cleanup_and_exit(GetLine *gl, SignalActions *si, int status)
{
/*
 * Restore the terminal to its original state before exiting the program.
 */
  gl_normal_io(gl);
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
  si = del_SignalActions(si);
/*
 * Exit the process.
 */
  exit(status);
}

/*.......................................................................
 * This is a signal-aware wrapper around the select() system call. It
 * is designed to facilitate reliable signal handling of a given set
 * of signals, without the race conditions that would usually surround
 * the use of select(). See the "RELIABLE SIGNAL HANDLING" section of
 * the gl_get_line(3) man page for further details.
 *
 * Provided that the calling function has blocked the specified set of
 * signals before calling this function, this function guarantees that
 * select() will be aborted by any signal that arrives between the
 * time that the caller blocked the specified signals and this
 * function returns. On return these signals will again be blocked to
 * prevent any signals that arrive after select() returns, from being
 * missed by the caller.
 *
 * Note that this function is written not to be specific to this
 * program, and is thus suitable for use in other programs, whether or
 * not they use gl_get_line().
 *
 * Also note that this function depends on the NSIG preprocessor
 * constant being >= the maximum number of signals available on the
 * host operating system. Under BSD and SysV, this macro is set
 * appropriately in signal.h. On other systems, a reasonably large
 * guess should be substituted. Although nothing terrible will happen
 * if a value that is too small is chosen, signal numbers that exceed
 * the specified value of NSIG will be ignored by this function. A
 * more robust method than depending on nsig would be to use the
 * POSIX sigismember() function to count valid signals, and use this
 * to allocate the array of sigaction structures used to preserve
 * 
 *
 * Input:
 *  n                   int    The number of file descriptors to pay
 *                             attention to at the start of each of the
 *                             following sets of file descriptors.
 *  readfds          fd_set *  The set of file descriptors to check for
 *                             readability, or NULL if not pertinent.
 *  wwritefds        fd_set *  The set of file descriptors to check for
 *                             writability, or NULL if not pertinent.
 *  exceptfds        fd_set *  The set of file descriptors to check for
 *                             the arrival of urgent data, or NULL if
 *                             not pertinent.
 *  timeout  struct timeval *  The maximum time that select() should
 *                             wait, or NULL to wait forever.
 *  mask           sigset_t *  The set of signals to catch.
 *  si       SignalHandlers *  An object in which to preserve temporary
 *                             copies signal handlers.
 * Output:
 *  return              int    > 0  The number of entries in all of the
 *                                  sets of descriptors that are ready
 *                                  for I/O.
 *                               0  Select() timed out.
 *                              -1  Error (see errno).
 */
static int demo_sigselect(int n, fd_set *readfds, fd_set *writefds,
			  fd_set *exceptfds, struct timeval *timeout,
			  sigset_t *mask, SignalActions *si)
{
/*
 * The reason that the the following variables are marked as volatile
 * is to prevent the compiler from placing their values in registers
 * that might not be saved and restored by sigsetjmp().
 */
  volatile sigset_t old_mask;  /* The displaced process signal mask */
  volatile int status;         /* The return value of select() */
/*
 * Make sure that all of the specified signals are blocked. This is
 * redundant if the caller has already blocked signals.
 */
  if(sigprocmask(SIG_BLOCK, mask, (sigset_t *) &old_mask) < 0)
    return -1;
/*
 * Record the fact that no signal has been caught yet.
 */
  demo_setjmp_signo = -1;
/*
 * Now set up the point where our temporary signal handlers will return
 * control if a signal is received.
 */
  if(sigsetjmp(demo_setjmp_buffer, 1) == 0) {
/*
 * Now install the temporary signal handlers that cause the above
 * sigsetjmp() to return non-zero when a signal is detected.
 */
    if(displace_signal_handlers(si, mask, demo_setjmp_handler)) {
      reinstate_signal_handlers(si);
      return 1;
    };
/*
 * Now that we are ready to catch the signals, unblock them.
 */
    sigprocmask(SIG_UNBLOCK, mask, NULL);
/*
 * At last, call select().
 */
    status = select(n, readfds, writefds, exceptfds, timeout);
/*
 * Block the specified signals again.
 */
    sigprocmask(SIG_BLOCK, mask, NULL);
/*
 * Record the fact that no signal was caught.
 */
    demo_setjmp_signo = -1;
  };
/*
 * We can get to this point in one of two ways. Either no signals were
 * caught, and the above block ran to completion (with demo_setjmp_signo=-1),
 * or a signal was caught that caused the above block to be aborted,
 * in which case demo_setjmp_signo will now equal the number of the signal that
 * was caught, and sigsetjmp() will have restored the process signal
 * mask to how it was before it was called (ie. all of the specified
 * signals blocked).
 *
 * First restore the signal handlers to how they were on entry to
 * this function.
 */
  reinstate_signal_handlers(si);
/*
 * Was a signal caught?
 */
  if(demo_setjmp_signo > 0) {
    sigset_t new_mask;
/*
 * Send the signal again, then unblock its delivery, so that the application's
 * signal handler gets invoked.
 */
    raise(demo_setjmp_signo);
    sigemptyset(&new_mask);
    sigaddset(&new_mask, demo_setjmp_signo);
    sigprocmask(SIG_UNBLOCK, &new_mask, NULL);
/*
 * Set the return status to show that a signal was caught.
 */
    errno = EINTR;
    status = -1;
  };
/*
 * Now restore the process signal mask to how it was on entry to this
 * function.
 */
  sigprocmask(SIG_SETMASK, (sigset_t *) &old_mask, NULL);
  return status;
}

/*.......................................................................
 * This is the main signal handler of this demonstration program. If a
 * SIGINT is received by the process, it arranges that the next call
 * to gl_get_line() will abort entry of the current line and start
 * entering a new one. Otherwise it calls the library function which
 * handles terminal resize signals and process suspension and process
 * termination signals. Both of the functions called by this signal
 * handler are designed to be async-signal safe, provided that the
 * rules laid out in the gl_io_mode(3) man page are followed.
 */
static void demo_signal_handler(int signo)
{
  if(signo==SIGINT)
    gl_abandon_line(demo_gl);
  else
    gl_handle_signal(signo, demo_gl, 1);
}

/*.......................................................................
 * The following signal handler is installed while select() is being
 * called from within a block of code protected by sigsetjmp(). It
 * simply records the signal that was caught in setjmp_signo, then
 * causes the sigsetjmp() to return non-zero.
 */
static void demo_setjmp_handler(int signo)
{
  demo_setjmp_signo = signo;
  siglongjmp(demo_setjmp_buffer, 1);
}

/*.......................................................................
 * This optional inactivity timeout function is used in this
 * demonstration to cause gl_get_line() to wait for a small amount of
 * time for I/O, before returning and allowing the event loop to
 * continue. This isn't needed if you want gl_get_line() to return
 * immediately, rather than blocking.
 */
static GL_TIMEOUT_FN(demo_timeout_fn)
{
  return GLTO_CONTINUE;
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
    "To the user this program appears to act identically to the main ",
    "demo program. However whereas the code underlying the main demo ",
    "program uses gl_get_line() in its default configuration, where each ",
    "call blocks the caller until the user has entered a complete input ",
    "line, demo3 uses gl_get_line() in its non-blocking server mode, ",
    "where it must be called repeatedly from an external ",
    "event loop to incrementally accept entry of the input ",
    "line, as and when terminal I/O becomes possible. The well commented ",
    "source code of demo3, which can be found in demo3.c, thus provides ",
    "a working example of how to use gl_get_line() in a manner that ",
    "doesn't block the caller. Documentation of this mode can be found ",
    "in the gl_io_mode(3) man page.\n"
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


/*.......................................................................
 * This is a constructor function for an object who's role is to allow
 * a signal handler to be assigned to potentially all available signals,
 * while preserving a copy of the original signal handlers, for later
 * restration.
 *
 * Output:
 *  return  SignalActions *  The new object, or NULL on error.
 */
static SignalActions *new_SignalActions(void)
{
  SignalActions *si;  /* The object to be returned */
/*
 * Allocate the container.
 */
  si = malloc(sizeof(SignalActions));
  if(!si) {
    fprintf(stderr, "new_SignalActions: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_SignalActions().
 */
  si->nsignal = 0;
  sigemptyset(&si->mask);
  si->actions = NULL;
/*
 * Count the number of signals that are available of the host
 * platform. Note that si->mask has no members set, and that
 * sigismember() is defined to return -1 if the signal number
 * isn't valid.
 */
  for(si->nsignal=1; sigismember(&si->mask, si->nsignal) == 0; si->nsignal++)
    ;
/*
 * Allocate the array of sigaction structures to use to keep a record
 * of displaced signal handlers.
 */
  si->actions = (struct sigaction *) malloc(sizeof(*si->actions) * si->nsignal);
  if(!si->actions) {
    fprintf(stderr, "Insufficient memory for %d sigaction structures.\n",
	    si->nsignal);
    return del_SignalActions(si);
  };
  return si;
}

/*.......................................................................
 * Delete a SignalActions object.
 *
 * Input:
 *  si     SignalActions *  The object to be deleted.
 * Output:
 *  return SignalActions *  The deleted object (always NULL).
 */
static SignalActions *del_SignalActions(SignalActions *si)
{
  if(si) {
    if(si->actions)
      free(si->actions);
    free(si);
  };
  return NULL;
}

/*.......................................................................
 * Replace the signal handlers of all of the signals in 'mask' with
 * the signal handler 'handler'.
 *
 * Input:
 *  si       SignalActions *     The object in which to record the displaced
 *                               signal handlers.
 *  mask          sigset_t *     The set of signals who's signal handlers
 *                               should be displaced.
 *  handler void (*handler)(int) The new signal handler to assign to each
 *                               of the signals marked in 'mask'.
 * Output:
 *  return             int       0 - OK.
 *                               1 - Error.
 */
static int displace_signal_handlers(SignalActions *si, sigset_t *mask,
				    void (*handler)(int))
{
  int signo;                /* A signal number */
  struct sigaction action;  /* The new signal handler */
/*
 * Mark the fact that so far we haven't displaced any signal handlers.
 */
  sigemptyset(&si->mask);
/*
 * Set up the description of the new signal handler. Note that
 * we make sa_mask=mask. This ensures that only one instance of the
 * signal handler will ever be running at one time.
 */
  action.sa_handler = handler;
  memcpy(&action.sa_mask, mask, sizeof(*mask));
  action.sa_flags = 0;
/*
 * Check each of the available signals to see if it is specified in 'mask'.
 * If so, install the new signal handler, record the displaced one in
 * the corresponding element of si->actions[], and make a record in
 * si->mask that this signal handler has been displaced.
 */
  for(signo=1; signo < si->nsignal; signo++) {
    if(sigismember(mask, signo)) {
      if(sigaction(signo, &action, &si->actions[signo]) < 0) {
	fprintf(stderr, "sigaction error (%s)\n", strerror(errno));
	return 1;
      };
      sigaddset(&si->mask, signo);
    };
  };
  return 0;
}

/*.......................................................................
 * Reinstate any signal handlers displaced by displace_signal_handlers().
 *
 * Input:
 *  sig   SignalActions *  The object containing the displaced signal
 *                         handlers.
 * Output:
 *  return          int    0 - OK.
 *                         1 - Error.
 */
static int reinstate_signal_handlers(SignalActions *si)
{
  int signo;                /* A signal number */
/*
 * Check each of the available signals to see if it is specified in
 * si->mask. If so, reinstate the displaced recorded in the
 * corresponding element of si->actions[], and make a record in
 * si->mask that this signal handler has been reinstated.
 */
  for(signo=1; signo < si->nsignal; signo++) {
    if(sigismember(&si->mask, signo)) {
      if(sigaction(signo, &si->actions[signo], NULL) < 0) {
	fprintf(stderr, "sigaction error (%s)\n", strerror(errno));
	return 1;
      };
      sigdelset(&si->mask, signo);
    };
  };
  return 0;
}
