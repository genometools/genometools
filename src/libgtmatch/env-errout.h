#ifndef ENV_ERROUT_H
#define ENV_ERROUT_H

#define ENVERROUT\
        if (env_error_is_set(env))\
        {\
          fprintf(stderr, "%s: %s\n",argv[0],env_error_get(env));\
        } else\
        {\
          fprintf(stderr, "%s: error expected\n",argv[0]);\
        }

#endif
