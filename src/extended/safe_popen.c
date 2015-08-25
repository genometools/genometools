/*
  Copyright (c) 2015 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include <string.h>
#include <unistd.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#include "core/ma_api.h"
#include "extended/safe_popen.h"

#ifndef _WIN32
static pid_t safe_fork(void) {
  pid_t childpid;

  if ((childpid = fork()) == (pid_t) -1)
    return (pid_t) -1;

  /* TODO DW Reseed PRNGs in both the parent and the child */

  /* If this is the parent process, there's nothing more to do */
  if (childpid != 0)
    return childpid;

  /* This is the child process */
  /* TODO DW Close all open files.  See
     https://www.safaribooksonline.com/library/view/secure-programming-cookbook/
     0596003943/ch01.html#secureprgckbk-CHP-1-SECT-1 */
  /* TODO DW Permanently drop privileges.  See
     https://www.safaribooksonline.com/library/view/secure-programming-cookbook/
     0596003943/ch01s03.html#secureprgckbk-CHP-1-FNOTE-3
  */

  return 0;
}
#endif

GtSafePipe *gt_safe_popen(const char *path,
                          char *const argv[],
                          char *const envp[],
                          GtError *err) {
#ifndef _WIN32
  int stdin_pipe[2], stdout_pipe[2], had_err = 0;
  GtSafePipe *p = NULL;

  p = gt_malloc(sizeof(*p));
  p->read_fd = p->write_fd = NULL;
  p->child_pid = (pid_t) -1;

  if ((had_err = pipe(stdin_pipe))) {
    gt_error_set(err, "could not open stdin pipe: %s", strerror(errno));
  }
  if (!had_err) {
    if ((had_err = pipe(stdout_pipe))) {
      gt_error_set(err, "could not open stdout pipe: %s", strerror(errno));
    }
    if (!had_err) {
      if (!(p->read_fd = fdopen(stdout_pipe[0], "r"))) {
        gt_error_set(err, "could not open stdout_pipe[0] for reading: %s",
                     strerror(errno));
        had_err = -1;
      }
      if (!had_err) {
        if (!(p->write_fd = fdopen(stdin_pipe[1], "w"))) {
          gt_error_set(err, "could not open stdin_pipe[1] for writing: %s",
                       strerror(errno));
          had_err = -1;
        }
        if (!had_err) {
          if ((p->child_pid = safe_fork()) == (pid_t) -1) {
            gt_error_set(err, "could not fork: %s", strerror(errno));
            had_err = -1;
          }
          if (!had_err) {
            if (!p->child_pid) {
              /* this is the child process */
              (void) close(stdout_pipe[0]);
              (void) close(stdin_pipe[1]);
              if (stdin_pipe[0] != 0) {
                (void) dup2(stdin_pipe[0], 0);
                (void) close(stdin_pipe[0]);
              }
              if (stdout_pipe[1] != 1) {
                (void) dup2(stdout_pipe[1], 1);
                (void) close(stdout_pipe[1]);
              }
              (void) execve(path, argv, envp);
              perror("could not execute external program: ");
              perror(strerror(errno));
              exit(127);
            }
            (void) close(stdout_pipe[1]);
            (void) close(stdin_pipe[0]);
          }
          if (had_err) {
            (void) fclose(p->write_fd);
          }
        }
        if (had_err) {
          (void) fclose(p->read_fd);
        }
      }
      if (had_err) {
        (void) close(stdout_pipe[1]);
        (void) close(stdout_pipe[0]);
      }
    }
    if (had_err) {
      (void) close(stdin_pipe[1]);
      (void) close(stdin_pipe[0]);
    }
  }
  if (had_err) {
    gt_free(p);
    p = NULL;
  }
  return p;
#else
  gt_error_set(err, "Function gt_safe_popen not implemented for windows yet");
  return NULL;
#endif
}

int gt_safe_pclose(GtSafePipe *p) {
#ifndef _WIN32
  int   status;
  pid_t pid = (pid_t) -1;

  if (p->child_pid != (pid_t) -1) {
    do {
      pid = waitpid(p->child_pid, &status, 0);
    } while (pid == (pid_t) -1 && errno == EINTR);
  }
  if (p->read_fd != NULL)
    (void) fclose(p->read_fd);
  if (p->write_fd != NULL)
    (void) fclose(p->write_fd);
  gt_free(p);
  if (pid != (pid_t) -1 && WIFEXITED(status))
    return WEXITSTATUS(status);
  else
    return (pid == (pid_t) -1 ? -1 : 0);
#else
  /* TODO: implement for Windows */
  gt_assert(0);
#endif
}
