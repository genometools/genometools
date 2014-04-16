$license = <<-END_LICENSE
/*
  Copyright (c) <%=year%> <%=name%> <<%=email%>>
  Copyright (c) <%=year%> Center for Bioinformatics, University of Hamburg

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
END_LICENSE

$repfile = <<-END_REPFILE

#include <stdlib.h>

#include "<%=subdir%>/<%=fkt_pref%>.h"

typedef struct <%=classN%>Class <%=classN%>Class;
typedef struct <%=classN%>Members <%=classN%>Members;

/* TODO: add typedefs for functions here */
/* examples showcasing coding convention:
typedef unsigned int (*<%=classN%>FooFunc)(<%=classN%>*, GtUword*);
typedef int *(*<%=classN%>BarFunc)(<%=classN%>*, GtUword);
*/
typedef void (*<%=classN%>DeleteFunc)(<%=classN%>*);

struct <%=classN%> {
  const <%=classN%>Class *c_class;
  <%=classN%>Members *members;
};

struct <%=classN%>Class {
  size_t size;
  /* XXfunctionNamesXX do not change will be replaced by function names */
};

struct <%=classN%>Members {\
<% if options.refcount %>\n  GtUword refcount;<% end %>
};

const <%=classN%>Class* gt_<%=fkt_pref%>_class_new(size_t size,
<%= " "*(classN.length + fkt_pref.length + 26) %>\
/* XXfunctionPtrXX do not change */);

<%=classN%>*            gt_<%=fkt_pref%>_create(const <%=classN%>Class*);

void*<%= ' ' * (classN.length + 8) %>\
gt_<%=fkt_pref%>_cast(const <%=classN%>Class*,
<%= ' ' * (classN.length + fkt_pref.length + 22) %><%=classN%>*);
END_REPFILE

$hwrapper = <<-END_HWRAPPER
#ifndef <%=guard_macro%>
#define <%=guard_macro%>
<%=content%>
#endif
END_HWRAPPER

$gt_class_c = <<-END_GT_CLASS_C

#include "core/ma.h"

struct <%=klass%>
{
};

<%=creatordecl%>
{
  <%=klass%> *<%=basename%>;
  <%=basename%> = gt_malloc(sizeof (<%=klass%>));
  return <%=basename%>;
}

<%=destructordecl%>
{
  gt_free(<%=basename%>);
}
END_GT_CLASS_C

$gt_class_h = <<-END_GT_CLASS_H

typedef struct <%=klass%> <%=klass%>;

<%=creatordecl%>;
<%=destructordecl%>;
END_GT_CLASS_H

$interface_header = <<-INTERFACE_HEADER

/* The <%=classN%> interface
   TODO: add documentation */
typedef struct <%=classN%> <%=classN%>;
<% if options.refcount %>\

/* Increases the reference count of the <<%=classN%>>. */
<%=classN%>*<%=' '*(max_type_len - 1 - classN.length)%>\
gt_<%=fkt_pref%>_ref(<%=classN%> *<%=fkt_pref%>);\
<% end %>\
<%functions.each do |func, parameters| %>

/* TODO: document and add parameter names */
<%  type =parameters[0] + ' '*(max_type_len - parameters[0].length) %>\
<%  type = type.gsub(/ [*]/, '* ')
    funcname = type + 'gt_' + fkt_pref + '_' +
      func[0..-5].gsub(/(^|.)([A-Z])/) do
        rep = ''
        rep = $1+'_' if $1.length >0
        rep += $2.downcase
      end + '('
%>\
<%=funcname%>\
<%  parameters[1..-1].each_with_index do |parameter, idx|%>\
<%    if idx > 0 %>,
<%=     ' ' * funcname.length%>\
<%    end %>\
<%    if parameter[-1] == '*' %>\
<% puts parameter %>\
<%=     parameter.chop%> *\
<%    else %>\
<%=     parameter%>\
<%    end %>\
<%    if parameter.match /^\#{classN}\\*$/ %>\
<%=     fkt_pref %>\
<%    else %>\
 para<%=idx%> /*TODO: name*/\
<%    end %>\
<%  end %>);\
<%end %>
INTERFACE_HEADER

$ref_func = <<-REF_FUNC
<%=classN%> *gt_<%=fkt_pref%>_ref(<%=classN%> *<%=fkt_pref%>)
{
  gt_assert(<%=fkt_pref%>);
  <%=fkt_pref%>->members->refcount++;
  return <%=fkt_pref%>;
}
REF_FUNC

$create_func = <<-CREATE_FUNC
<%=classN%> *gt_<%=fkt_pref%>_create(const <%=classN%>Class *<%=fkt_pref%>_c)
{
  <%=classN%> *<%=fkt_pref%>;
  gt_assert(<%=fkt_pref%>_c && <%=fkt_pref%>_c->size);
  <%=fkt_pref%> = gt_calloc((size_t) 1, <%=fkt_pref%>_c->size);
  <%=fkt_pref%>->c_class = <%=fkt_pref%>_c;
  <%=fkt_pref%>->members = gt_calloc((size_t) 1, sizeof (<%=classN%>Members));
  return <%=fkt_pref%>;
}
CREATE_FUNC

$cast_func = <<-CAST_FUNC
void *gt_<%=fkt_pref%>_cast(GT_UNUSED const <%=classN%>Class *<%=fkt_pref%>_c,
<%=' ' * (fkt_pref.length + 15)%><%=classN%> *<%=fkt_pref%>)
{
  gt_assert(<%=fkt_pref%>_c && <%=fkt_pref%> &&
            <%=fkt_pref%>->c_class == <%=fkt_pref%>_c);
  return <%=fkt_pref%>;
}
CAST_FUNC

$delete_func = <<-DELETE_FUNC
void gt_<%=fkt_pref%>_delete(<%=classN%> *<%=fkt_pref%>)
{
  if (<%=fkt_pref%> != NULL) {
<% if options.refcount %>\
    if (<%=fkt_pref%>->members->refcount) {
      <%=fkt_pref%>->members->refcount--;
      return;
    }
<%end %>\
    gt_assert(<%=fkt_pref%>->c_class);
    if (<%=fkt_pref%>->c_class->delete_func != NULL)
      <%=fkt_pref%>->c_class->delete_func(<%=fkt_pref%>);
    gt_free(<%=fkt_pref%>->members);
    gt_free(<%=fkt_pref%>);
  }
}
DELETE_FUNC

$interface_func = <<-INTERFACE_FUNC
<%funcname = type + 'gt_' + fkt_pref + '_' +
func[0..-5].gsub(/(^|.)([A-Z])/) do
  rep = ''
  rep = $1+'_' if $1.length >0
  rep += $2.downcase
end + '(' %>
<%=funcname %>\
<%paras.each_with_index do |para,idx| %>\
<%  if idx > 0 %>,
<%=   ' ' * (funcname.length) %>\
<%  end %>\
<%  if para[-1] == '*' %>\
<%=   para.chop %> *\
<%  else %>\
<%=   para %>\
<%  end %>\
<%  if para.match /^\#{classN}\\*$/ %>\
<%=   fkt_pref %>\
<%  else %>\
 para<%=idx%> /*TODO: name*/\
<%  end %>\
<%end %>)
{
  gt_assert(<%=fkt_pref%> != NULL);
  gt_assert(<%=fkt_pref%>->c_class != NULL);
  if (<%=fkt_pref%>->c_class->\
<%=func.gsub(/(^|.)([A-Z])/) do
  rep = ''
  rep = $1+'_' if $1.length >0
  rep += $2.downcase \
end %> != NULL)
    return <%=fkt_pref%>->c_class->\
<%=func.gsub(/(^|.)([A-Z])/) do
  rep = ''
  rep = $1+'_' if $1.length >0
  rep += $2.downcase \
end %>(\
<%paras.each_with_index do |para,idx| %>\
<%  if idx > 0 %>,
<%=   ' ' * (func.length + fkt_pref.length + 23) %>\
<%  end%>\
<%  if para.match /^\#{classN}\\*$/ %>\
<%=   fkt_pref %>\
<%  else %>\
 para<%=idx%> /*TODO: name*/\
<%  end %>\
<%end %>);
  return\
<%if type.match /[*]/ %>\
 NULL;\
<%elsif type.match /void/ %>\
;\
<%else %>\
 0;\
<%end %> /* TODO: check if default return value is sane */
}
INTERFACE_FUNC

$class_new_fkt = <<-CLASS_NEW_FKT
const <%=classN%>Class *gt_<%=fkt_pref%>_class_new(size_t size\
<%functions.each do |func, paras| %>,
<%=  ' ' * (27 + classN.length + fkt_pref.length)%>\
<%=classN%><%=func%> <%=func.gsub(/(^|.)([A-Z])/) do
  rep = ''
  rep = $1+'_' if $1.length > 0
  rep += $2.downcase
end%>\
<%end %>)
{
  <%=classN%>Class *<%=fkt_pref%>_c = gt_class_alloc(sizeof (*<%=fkt_pref%>_c));
  <%=fkt_pref%>_c->size = size;
<%functions.each do |func, paras| %>\
<% name = func.gsub(/(^|.)([A-Z])/) do
  rep = ''
  rep = $1+'_' if $1.length >0
  rep += $2.downcase
end %>\
  <%=fkt_pref%>_c-><%=name%> = <%=name%>;
<%end %>\
  return <%=fkt_pref%>_c;
}
CLASS_NEW_FKT

$interface_file = <<-INTERFACE_CODE
#include "<%=subdir%>/<%=fkt_pref%>_rep.h"

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"

<%=create_func%>
<% if options.refcount %>\
<%=ref_func%>
<% end %>\
<%=cast_func%>\
<%=interface_funcs%>
<%=class_new_fkt%>
<%=delete_func%>
INTERFACE_CODE

$implement_header = <<-IMPLEMENT_HEADER

#include "<%=classbase%>.h"

/* The <<%=iclassN%>> class implements the <<%=classN%>> interface.
   TODO: add documentation */
typedef struct <%=iclassN%> <%=iclassN%>;

/* Return a new <<%=classN%>> object, the implementation beeing of type
   <<%=iclassN%>>.
   TODO: add documentation */
<%=classN%>* <%=ifkt_pref%>_new(void /*TODO: add parameters*/);

/* TODO: add prototypes of implementation <<%=iclassN%>> only functions */
IMPLEMENT_HEADER

$implement_src = <<-IMPLEMENT_SRC

#include "<%=classbase%>_rep.h"
#include "<%=basename%>.h"
#include "core/unused_api.h"

struct <%=iclassN%> {
  <%=classN%> parent_instance;
  /* TODO: add implementation members */
};

#define <%=ifkt_pref%>_cast(cvar) \\
        <%=fkt_pref%>_cast(<%=ifkt_pref%>_class(), cvar)

<%functions.each do |func, parameters| %>\
<%  type = parameters[0] %>\
<%  funcnames << ifkt_pref + '_' +
      func[0..-5].gsub(/(^|.)([A-Z])/) do
        rep = ''
        rep = $1+'_' if $1.length >0
        rep += $2.downcase
      end
    funcname = type + funcnames[-1] + '('
%>\
static <%=funcname%>\
<%  parameters[1..-1].each_with_index do |parameter, idx|%>\
<%    if idx > 0 %>,
<%=     ' ' * (funcname.length + 7)%>\
<%    end %>GT_UNUSED \
<%    if parameter[-1] == '*' %>\
<%=     parameter.chop%> *\
<%    else %>\
<%=     parameter%>\
<%    end %>\
<%    if parameter.match /\#{classN}\\*$/ %>\
<%=     icvar %>\
<%    else %>\
 para<%=idx%> /*TODO: name*/\
<%    end %>\
<%  end %>)
{
<% if parameters[1].match classN %>\
<%=  iclassN%> <%=iclassN.gsub(/^Gt/,'').gsub(/[a-z_]/,'').downcase%>\
= <%=ifkt_pref%>_cast(<%=icvar%>);
<% end %>\
  /* TODO: add functionality */
  return\
<%  if type.match /[*]/ %>\
 NULL;\
<%  elsif type.match /void/ %>\
;\
<%  else %>\
 0;\
<%  end %> /* TODO: check if default return value is sane */
}

<%end %>\
/* TODO: implement <<%=iclassN%>> only functions */

/* map static local methods to interface */
const <%=classN%>Class* <%=ifkt_pref%>_class(void)
{
  static const <%=classN%>Class *this_c = NULL;
  if (this_c == NULL) {
    this_c = <%=fkt_pref%>_class_new(sizeof (<%=iclassN%>),
<%funcnames.each_with_index do |name,idx| %>\
<%  if idx > 0 %>,
<%  end %>\
<%= ' ' * (fkt_pref.length + 24) %><%=name%>\
<%end %>);
  }
  return this_c;
}

<%=classN%>* <%=ifkt_pref%>_new(void /* TODO: add parameters */)
{
  <%=classN%> *<%=cvar%>;
  <%=iclassN%> *<%=icvar%>;
  <%=cvar%> = <%=fkt_pref%>_create(<%=ifkt_pref%>_class());
  <%=icvar%> = <%=ifkt_pref%>_cast(<%=cvar%>);
  /* TODO: initialise */
  return <%=cvar%>;
}
IMPLEMENT_SRC
