#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"

#include "esa-maxpairs.pr"

typedef struct
{
  unsigned int userdefinedleastlength;
  Str *indexname;
} Maxpairsoptions;

static OPrval parse_options(Maxpairsoptions *maxpairsoptions,
                            int *parsed_args,
                            int argc, const char **argv,Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[-l minlength] -db indexname",
                         "Enumerate maximal paris of minimum length.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  option = option_new_uint_min("l","Specify minimum length",
                               &maxpairsoptions->userdefinedleastlength,
                               (unsigned int) 20,
                               (unsigned int) 1,env);
  option_parser_add_option(op, option, env);

  option = option_new_string("ii",
                             "Specify input index",
                             maxpairsoptions->indexname, NULL, env);
  option_parser_add_option(op, option, env);
  option_is_mandatory(option);

  oprval = option_parser_parse(op, parsed_args, argc, argv,
                               versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_maxpairs(int argc, const char **argv, Env *env)
{
  bool haserr = false;
  int parsed_args;
  Maxpairsoptions maxpairsoptions;
  OPrval oprval;

  env_error_check(env);

  maxpairsoptions.indexname = str_new(env);
  oprval = parse_options(&maxpairsoptions,&parsed_args, argc, argv, env);
  if(oprval == OPTIONPARSER_OK)
  {
    assert(parsed_args == argc);
    if(callenummaxpairs(maxpairsoptions.indexname,
                        (uint32_t) maxpairsoptions.userdefinedleastlength,
                        env) != 0)
    {
      haserr = true;
    }
  }
  str_delete(maxpairsoptions.indexname,env);
  if(oprval == OPTIONPARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if(oprval == OPTIONPARSER_ERROR)
  {
    return -1;
  }
  return haserr ? -1 : 0;
}
