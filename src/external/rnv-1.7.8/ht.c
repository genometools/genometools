/* $Id: ht.c,v 1.14 2004/01/23 20:26:45 dvd Exp $ */

#include <stdlib.h> /*NULL*/
#include <assert.h> /*assert*/
#include "m.h"
#include "ht.h"

#define LOAD_FACTOR 2

void ht_init(struct hashtable *ht,int len,int (*hash)(int),int (*equal)(int,int)) {
  assert(len>0);
  ht->tablen=1; len*=LOAD_FACTOR;
  while(ht->tablen<len) ht->tablen<<=1;
  ht->limit=ht->tablen/LOAD_FACTOR;
  ht->table=(int*)m_alloc(ht->tablen<<1,sizeof(int)); /* the second half is hash values */
  ht->hash=hash; ht->equal=equal;
  ht_clear(ht);
}

void ht_clear(struct hashtable *ht) {
  int i;
  ht->used=0; for(i=0;i!=ht->tablen;++i) ht->table[i]=-1;
}

void ht_dispose(struct hashtable *ht) {
  m_free(ht->table); ht->table=NULL;
}

#define first(ht,hv) (hv&(ht->tablen-1))
#define next(ht,i) (i==0?ht->tablen-1:i-1)

int ht_get(struct hashtable *ht,int i) {
  int hv=ht->hash(i),j;
  for(j=first(ht,hv);;j=next(ht,j)) {
    int tj=ht->table[j];
    if(tj==-1) break;
    if(ht->equal(i,tj)) return tj;
  }
  return -1;
}

void ht_put(struct hashtable *ht,int i) {
  int hv=ht->hash(i),j;
  if(ht->used==ht->limit) {
    int tablen=ht->tablen; int *table=ht->table;
    ht->tablen<<=1; ht->limit<<=1;
    ht->table=(int*)m_alloc(ht->tablen<<1,sizeof(int));
    for(j=0;j!=ht->tablen;++j) ht->table[j]=-1;
    for(j=0;j!=tablen;++j) {
      if(table[j]!=-1) {
	int hvj=table[j|tablen]; int k;
	for(k=first(ht,hvj);ht->table[k]!=-1;k=next(ht,k));
	ht->table[k]=table[j]; ht->table[k|ht->tablen]=hvj;
      }
    }
    m_free(table);
  }
  for(j=first(ht,hv);ht->table[j]!=-1;j=next(ht,j)) assert(!ht->equal(i,ht->table[j]));
  ht->table[j]=i;
  ht->table[ht->tablen|j]=hv;
  ++ht->used;
}

static int del(struct hashtable *ht,int i,int eq) {
  if(ht->used!=0) {
    int hv=ht->hash(i),j;
    for(j=first(ht,hv);;j=next(ht,j)) {
      int tj=ht->table[j];
      if(tj==-1) break;
      if(eq?i==tj:ht->equal(i,tj)) {
	do {
	  int k=j,j0;
	  ht->table[j]=-1;
	  for(;;) {
	    j=next(ht,j);
	    if(ht->table[j]==-1) break;
	    j0=first(ht,ht->table[j|ht->tablen]);
	    if((k<=j0||j0<j)&&(j0<j||j<=k)&&(j<=k||k<=j0)) break;
	  }
	  ht->table[k]=ht->table[j]; ht->table[k|ht->tablen]=ht->table[j|ht->tablen];
	} while(ht->table[j]!=-1);
	--ht->used;
	return tj;
      }
    }
  }
  return -1;
}
int ht_del(struct hashtable *ht,int i) {return del(ht,i,0);}
int ht_deli(struct hashtable *ht,int i) {return del(ht,i,1);}
