/* $Id: b64gen.c 453 2005-01-07 00:07:51Z dvd $ */
/* prints definitions for trailing characters in BASE64 pattern (see xsd.c) */
int main() {
  char b64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  int i;
  printf("#define PAT_BASE64_2 \"");
  for(i=0;i!=64;++i)
    if((i&0x30)==i) printf("%c",b64[i]);
  printf("\"\n");
  printf("#define PAT_BASE64_1 \"");
  for(i=0;i!=64;++i)
    if((i&0x3C)==i) printf("%c",b64[i]);
  printf("\"\n");
}
