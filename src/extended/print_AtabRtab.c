#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/error.h"
#include "core/types_api.h"
#include "affinealign_linear.h"

typedef enum {
  R,
  D,
  I,
  X //unknown
} Edge;

typedef struct {
  GtUword Rvalue, Dvalue, Ivalue;
  
  Edge Redge,
       Dedge,
       Iedge;
}Atabentry;

typedef struct {
  GtUword idx;
  Edge edge;
}Rnode;

typedef struct {
  Rnode R,D,I;
}Rtabentry;

/*static char* edge_to_char(Edge edge)
{
  switch (edge){
    case 0: return "R";
    case 1: return "D";
    case 2: return "I";
    default: return "X";
    }
}*/

/*void print(Atabentry *Atabcolumn, Rtabentry *Rtabcolumn,const GtUword ulen, const GtUword colindex)
{
    FILE *data_A, *data_R;
    data_A =fopen("data_A","a");
    data_R =fopen("data_R","a");
    Atabentry *a;
    Rtabentry *r;
    fprintf(data_A, "****************************************\n"
                     GT_WU"\n"
                     "******************************************\n", 
                     colindex);
    fprintf(data_R, "****************************************\n"
                     GT_WU"\n"
                     "******************************************\n", 
                     colindex);
    for(a = Atabcolumn;a <= Atabcolumn+ulen ; a++)
    {   
      fprintf(data_A,"Atab[%d].R: ("GT_WU", %s) I: ("GT_WU", %s) D: ("GT_WU", %s)\n",
                    (int)(a-Atabcolumn),
                    a->Rvalue, edge_to_char(a->Redge),
                    a->Ivalue, edge_to_char(a->Iedge),
                    a->Dvalue, edge_to_char(a->Dedge));
    }
    for(r = Rtabcolumn;r <= Rtabcolumn+ulen ; r++)
    {   
      fprintf(data_R,"Rtab[%d].R:("GT_WU", %s) I: ("GT_WU", %s) D: ("GT_WU", %s)\n",
                  (int)(r-Rtabcolumn),
                  r->R.idx, edge_to_char(r->R.edge),
                  r->I.idx, edge_to_char(r->I.edge),
                  r->D.idx, edge_to_char(r->D.edge));
    }
    fclose(data_A);
    fclose(data_R);
}*/
