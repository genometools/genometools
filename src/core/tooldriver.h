/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TOOLDRIVER_H
#define TOOLDRIVER_H

#include "core/error.h"
#include "core/tool_api.h"

/* The prototype of a tool function. */
typedef int (*GtToolFunc)(int argc, const char **argv, GtError *err);

typedef struct GtLicense GtLicense;

typedef GtLicense* (*GtLicenseConstructor)(const char *argv0,
                                           unsigned int major_version,
                                           unsigned int minor_version);
typedef void       (*GtLicenseDestructor)(GtLicense*);

/* The tool driver module allows one to compile a tool into a separate binary.
   This is mostly useful for stand-alone applications like GenomeThreader.
   The tool driver creates an GtError object, calls <tool>, and reports errors.
*/
int gt_tooldriver(GtToolFunc tool, int argc, char *argv[]);

int gt_tooldriver_with_license(GtToolFunc tool, int argc, char *argv[],
                               GtLicense **license_out,
                               unsigned int major_version,
                               unsigned int minor_version,
                               GtLicenseConstructor, GtLicenseDestructor);

int gt_toolobjdriver(GtToolConstructor, int argc, char *argv[]);

int gt_toolobjdriver_with_license(GtToolConstructor tool_constructor, int argc,
                                  char *argv[], GtLicense **license_out,
                                  unsigned int major_version,
                                  unsigned int minor_version,
                                  GtLicenseConstructor, GtLicenseDestructor);

#endif
